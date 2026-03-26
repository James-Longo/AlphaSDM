import ee
import pandas as pd
import math
import time
import os
import json
import argparse

def initialize_gee_safe(project_id='927064910362'):
    """
    Robust Earth Engine initialization with manual credentials fallback.
    """
    try:
        # Check for local configuration override (e.g. from gcloud or earthengine-api)
        creds_path = os.path.expanduser('~/.config/earthengine/credentials')
        if os.path.exists(creds_path):
            with open(creds_path, 'r') as f:
                c = json.load(f)
                project_id = c.get('project', project_id)

        ee.Initialize(project=project_id)
        print(f"Success: GEE Initialized (Project: {project_id})")
        return True
    except Exception as e:
        print(f"Error: GEE Initialization failed: {e}")
        return False

def extract_embeddings(csv_path, output_prefix, drive_folder, lon_col='longitude', lat_col='latitude', year_col='year', chunk_size=25000):
    """
    Extract Alpha Earth satellite embeddings for any list of coordinates.
    """
    if not initialize_gee_safe():
        return

    # 1. Load Data
    print(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)
    
    # 2. Flexible Column Detection
    rename_dict = {}
    cols = df.columns.tolist()
    
    # Lon case-insensitive check
    if lon_col not in cols:
        for c in ['lon', 'long', 'x']:
            if c in [x.lower() for x in cols]:
                rename_dict[[x for x in cols if x.lower() == c][0]] = lon_col
                break
    
    # Lat case-insensitive check
    if lat_col not in cols:
        for c in ['lat', 'y']:
            if c in [x.lower() for x in cols]:
                rename_dict[[x for x in cols if x.lower() == c][0]] = lat_col
                break

    # Year case-insensitive check
    if year_col not in cols:
        for c in ['year', 'time', 'date']:
            if c in [x.lower() for x in cols]:
                rename_dict[[x for x in cols if x.lower() == c][0]] = year_col
                break

    if rename_dict:
        print(f"Auto-detected columns: {rename_dict}")
        df = df.rename(columns=rename_dict)
    
    if not all(c in df.columns for c in [lon_col, lat_col, year_col]):
        print(f"Error: Missing columns. Required: {lon_col}, {lat_col}, {year_col}. Found: {df.columns.tolist()}")
        return

    # 3. Deduplicate (Essential for 10m global datasets)
    initial_len = len(df)
    df = df[[lon_col, lat_col, year_col]].drop_duplicates().reset_index(drop=True)
    if initial_len != len(df):
        print(f"Dropped {initial_len - len(df)} duplicate coordinate triplets.")
    
    print(f"Processing {len(df)} unique points.")

    # 4. Year-by-Year Export
    years = sorted(df[year_col].unique())
    
    for yr in years:
        yr_df = df[df[year_col] == yr].copy()
        if len(yr_df) == 0: continue
        print(f"\n--- Year {int(yr)} ({len(yr_df)} points) ---")
        
        # Spatial Sorting: Ensures better tile-fetching locality on GEE side
        yr_df = yr_df.sort_values([lon_col, lat_col]).reset_index(drop=True)
        
        n_points = len(yr_df)
        n_chunks = math.ceil(n_points / chunk_size)

        emb_cols = [f"A{i:02d}" for i in range(64)]
        img = ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL") \
                .filter(ee.Filter.calendarRange(int(yr), int(yr), "year")) \
                .mosaic() \
                .select(emb_cols)

        for i in range(n_chunks):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, n_points)
            sub_df = yr_df.iloc[start:end]
            
            task_name = f"{output_prefix}_{int(yr)}_batch{i:03d}"
            print(f"  -> {i+1}/{n_chunks}: Queuing '{task_name}' ({len(sub_df)} points)...")

            coords = sub_df[[lon_col, lat_col]].values.tolist()
            
            # Use MultiPoint geometries for efficient upload
            features = []
            SUB_CHUNK_LIMIT = 5000
            for j in range(0, len(coords), SUB_CHUNK_LIMIT):
                batch_coords = coords[j:j+SUB_CHUNK_LIMIT]
                features.append(ee.Feature(ee.Geometry.MultiPoint(batch_coords), {year_col: int(yr)}))

            # Unpack and Sample server-side
            fc = ee.FeatureCollection(features).map(lambda f: 
                ee.FeatureCollection(
                    ee.List(f.geometry().coordinates()).map(lambda c: 
                        ee.Feature(ee.Geometry.Point(c), {
                            year_col: f.get(year_col),
                            lon_col: ee.List(c).get(0),
                            lat_col: ee.List(c).get(1)
                        })
                    )
                )
            ).flatten()

            sampled = img.sampleRegions(
                collection=fc,
                properties=[year_col, lon_col, lat_col],
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
                time.sleep(0.5)
            except Exception as e:
                print(f"    ! Queuing failed: {e}. Waiting 10s...")
                time.sleep(10)
                task.start()

    print("\nSuccess: All tasks queued. Monitor progress at: https://code.earthengine.google.com/tasks")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generalized Alpha Earth Embedding Extractor (10m Resolution).")
    parser.add_argument("--csv", required=True, help="Path to input coordinate CSV")
    parser.add_argument("--prefix", required=True, help="Description prefix for GEE tasks/output files")
    parser.add_argument("--folder", required=True, help="Target Google Drive folder name")
    parser.add_argument("--lon", default="longitude", help="Longitude column name (auto-detects 'lon', 'long', 'x')")
    parser.add_argument("--lat", default="latitude", help="Latitude column name (auto-detects 'lat', 'y')")
    parser.add_argument("--year", default="year", help="Year column name (auto-detects 'time', 'date')")
    parser.add_argument("--chunk", type=int, default=25000, help="Points per GEE export task (default: 25k)")
    
    args = parser.parse_args()
    
    extract_embeddings(
        csv_path=args.csv,
        output_prefix=args.prefix,
        drive_folder=args.folder,
        lon_col=args.lon,
        lat_col=args.lat,
        year_col=args.year,
        chunk_size=args.chunk
    )
