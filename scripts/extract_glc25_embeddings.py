import ee
import pandas as pd
import math
import time
import os
import json

def initialize_gee_safe():
    """
    Robust Earth Engine initialization with manual credentials fallback.
    """
    project_id = '927064910362'
    
    try:
        ee.Initialize(project=project_id)
        print(f"Success: GEE Initialized (Project: {project_id})")
        return True
    except Exception:
        try:
            from google.oauth2.credentials import Credentials
            creds_path = os.path.expanduser('~/.config/earthengine/credentials')
            if not os.path.exists(creds_path):
                creds_path = os.path.expanduser('~/.google/earthengine/credentials')
            
            if os.path.exists(creds_path):
                with open(creds_path, 'r') as f:
                    c = json.load(f)
                    client_id = c.get('client_id', '517222506229-vsmmajv00ul0bs7p89v5m89qs8eb9359.apps.googleusercontent.com')
                    client_secret = c.get('client_secret')
                    refresh_token = c.get('refresh_token')
                    
                    if refresh_token:
                        creds = Credentials(
                            token=None,
                            refresh_token=refresh_token,
                            client_id=client_id,
                            client_secret=client_secret,
                            token_uri='https://oauth2.googleapis.com/token',
                            scopes=c.get('scopes', ['https://www.googleapis.com/auth/earthengine'])
                        )
                        ee.Initialize(credentials=creds, project=project_id)
                        print(f"Success: GEE Initialized via manual credentials (Project: {project_id})")
                        return True
        except Exception as e:
            print(f"Error: Manual credentials load failed: {e}")
            
    print("\nGEE Initialization failed. Please run: earthengine authenticate --force")
    return False

def extract_glc25_embeddings(csv_path, output_name_prefix, drive_folder="GLC25_embeddings"):
    """
    Robustly extract Alpha Earth Embeddings for GLC25 data from GEE.
    Uses Spatial Sorting to ensure 25k chunks remain small AND spatially grouped.
    """
    if not initialize_gee_safe():
        return

    # 2. Load and Deduplicate Data
    print(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)
    
    rename_dict = {}
    if 'lon' in df.columns and 'longitude' not in df.columns: rename_dict['lon'] = 'longitude'
    if 'lat' in df.columns and 'latitude' not in df.columns: rename_dict['lat'] = 'latitude'
    if rename_dict:
        df = df.rename(columns=rename_dict)
    
    if not all(c in df.columns for c in ['longitude', 'latitude', 'year']):
        print(f"Error: Missing columns. Found: {df.columns.tolist()}")
        return

    # Deduplicate: Essential for 5.1M row datasets to avoid redundant calls
    df = df[['longitude', 'latitude', 'year']].drop_duplicates().reset_index(drop=True)
    print(f"Processing {len(df)} unique coordinate triplets.")

    # 3. Process Year-by-Year
    years = sorted(df['year'].unique())
    CHUNK_SIZE = 25000  # Smaller chunks = zero "Computed value too large" errors

    for yr in years:
        yr_df = df[df['year'] == yr].copy()
        
        # SPATIAL SORTING: Slices of the sorted list will form longitudinal/latitudinal strips.
        # This keeps tile-fetching localized per task while guaranteeing fixed task size.
        print(f"\n--- Year {int(yr)} ({len(yr_df)} points) ---")
        yr_df = yr_df.sort_values(['longitude', 'latitude']).reset_index(drop=True)
        
        n_points = len(yr_df)
        n_chunks = math.ceil(n_points / CHUNK_SIZE)

        # Embedding image setup
        emb_cols = [f"A{i:02d}" for i in range(64)]
        img = ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL") \
                .filter(ee.Filter.calendarRange(int(yr), int(yr), "year")) \
                .mosaic() \
                .select(emb_cols)

        for i in range(n_chunks):
            start = i * CHUNK_SIZE
            end = min((i + 1) * CHUNK_SIZE, n_points)
            sub_df = yr_df.iloc[start:end]
            
            task_name = f"{output_name_prefix}_{int(yr)}_batch{i:03d}"
            print(f"  -> {i+1}/{n_chunks}: Queuing '{task_name}' ({len(sub_df)} points)...")

            coords = sub_df[['longitude', 'latitude']].values.tolist()
            
            features = []
            SUB_CHUNK_LIMIT = 5000
            for j in range(0, len(coords), SUB_CHUNK_LIMIT):
                batch_coords = coords[j:j+SUB_CHUNK_LIMIT]
                features.append(ee.Feature(ee.Geometry.MultiPoint(batch_coords), {'year': int(yr)}))

            # Unpack for sampling
            fc = ee.FeatureCollection(features).map(lambda f: 
                ee.FeatureCollection(
                    ee.List(f.geometry().coordinates()).map(lambda c: 
                        ee.Feature(ee.Geometry.Point(c), {
                            'year': f.get('year'),
                            'longitude': ee.List(c).get(0),
                            'latitude': ee.List(c).get(1)
                        })
                    )
                )
            ).flatten()

            # Sampling results
            sampled = img.sampleRegions(
                collection=fc,
                properties=['year', 'longitude', 'latitude'],
                scale=10,
                tileScale=16,
                geometries=False
            ).filter(ee.Filter.notNull(['A00']))

            task = ee.batch.Export.table.toDrive(
                collection=sampled,
                description=task_name,
                folder=drive_folder,
                fileFormat='CSV'
            )
            
            try:
                task.start()
                time.sleep(1) # Faster queuing
            except Exception as e:
                print(f"    ! Queuing failed: {e}. Waiting 10s...")
                time.sleep(10)
                task.start()

    print("\nSuccess: All clustered and sorted tasks queued! Monitor at: https://code.earthengine.google.com/tasks")

if __name__ == "__main__":
    csv_path = "/home/james-longo/Projects/AlphaSDM_benchmarks/GeoLifeCLEF 2025/GLC25_PO_metadata_train.csv"
    extract_glc25_embeddings(csv_path, "GLC25_PO_Final_Batch", drive_folder="AlphaSDM_GLC25_FixedSize_Extraction")
