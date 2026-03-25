import ee
import pandas as pd
import math
import time
import os

def extract_glc25_embeddings(csv_path, output_name_prefix, drive_folder="GLC25_embeddings"):
    """
    Robustly extract Alpha Earth Embeddings for GLC25 data from GEE.
    Splits tasks by year and by small chunks to prevent GEE OOM/Complexity errors.
    """
    # 1. Initialize GEE
    try:
        ee.Initialize()
    except Exception as e:
        print("GEE Initialization failed. Ensure you've authenticated with 'earthengine authenticate'.")
        return

    # 2. Load Data
    print(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)
    
    # Handle column names from the GLC2025 CSV
    rename_dict = {}
    if 'lon' in df.columns and 'longitude' not in df.columns: rename_dict['lon'] = 'longitude'
    if 'lat' in df.columns and 'latitude' not in df.columns: rename_dict['lat'] = 'latitude'
    if rename_dict:
        df = df.rename(columns=rename_dict)
    
    # Required columns
    if not all(c in df.columns for c in ['longitude', 'latitude', 'year']):
        print(f"Error: Required columns (longitude, latitude, year) not found.")
        print(f"Available columns: {df.columns.tolist()}")
        return

    # 3. Process Year-by-Year
    years = sorted(df['year'].unique())
    print(f"Detected years: {years}")

    CHUNK_SIZE = 100000  # 100k points per task is a safe sweet spot for 64 bands

    for yr in years:
        yr_df = df[df['year'] == yr].copy().reset_index(drop=True)
        n_points = len(yr_df)
        n_chunks = math.ceil(n_points / CHUNK_SIZE)
        print(f"\n--- Processing Year {int(yr)} ({n_points} points, {n_chunks} tasks) ---")

        # Get the embedding image for this year
        emb_cols = [f"A{i:02d}" for i in range(64)]
        img = ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL") \
                .filter(ee.Filter.calendarRange(int(yr), int(yr), "year")) \
                .mosaic() \
                .select(emb_cols)

        for i in range(n_chunks):
            start = i * CHUNK_SIZE
            end = min((i + 1) * CHUNK_SIZE, n_points)
            sub_df = yr_df.iloc[start:end]
            
            task_name = f"{output_name_prefix}_{int(yr)}_batch{i:02d}"
            print(f"  -> Task {i+1}/{n_chunks}: Queuing '{task_name}' ({len(sub_df)} points)...")

            # Efficient MultiPoint upload for the chunk
            # Instead of ee.Feature(Point), we upload as a MultiPoint list to minimize tree size
            # Then we zip them on the server side (like in R/gee_logic.R but in Python)
            
            coords = sub_df[['longitude', 'latitude']].values.tolist()
            years_list = [int(yr)] * len(sub_df) # though it's all one year here
            
            # We can use sampleRegions on a FeatureCollection created from points
            # To stay efficient, we split the sub_df into further internal multipoint features if needed
            # but for 100k it should be fine to do direct if we can.
            
            features = []
            # Internal sub-chunking of geometries for GeoJSON limit (50 MB)
            SUB_CHUNK_LIMIT = 5000
            for j in range(0, len(coords), SUB_CHUNK_LIMIT):
                batch_coords = coords[j:j+SUB_CHUNK_LIMIT]
                # MultiPoint is much more compact than thousands of Feature(Point)
                geom = ee.Geometry.MultiPoint(batch_coords)
                features.append(ee.Feature(geom, {'year': int(yr), 'batch_idx': i}))

            fc = ee.FeatureCollection(features).map(lambda f: 
                ee.FeatureCollection(
                    ee.List(f.geometry().coordinates()).map(lambda c: 
                        ee.Feature(ee.Geometry.Point(c), {'year': f.get('year')})
                    )
                )
            ).flatten()

            # Apply sampling
            sampled = img.sampleRegions(
                collection=fc,
                properties=['year'],
                scale=10,
                tileScale=16,
                geometries=False
            ).filter(ee.Filter.notNull(['A00']))

            # Start export
            task = ee.batch.Export.table.toDrive(
                collection=sampled,
                description=task_name,
                folder=drive_folder,
                fileFormat='CSV'
            )
            
            try:
                task.start()
                time.sleep(1) # Breathe
            except Exception as e:
                print(f"    ! Error queuing task: {e}")
                time.sleep(10) # Wait and try again if it's a rate limit
                task.start()

    print("\nAll tasks queued! You can monitor them at https://code.earthengine.google.com/tasks")

if __name__ == "__main__":
    import sys
    csv_path = "/home/james-longo/Projects/autoSDM_benchmarks/GeoLifeCLEF 2025/GLC25_PO_metadata_train.csv"
    extract_glc25_embeddings(csv_path, "GLC25_PO_Embeddings", drive_folder="AlphaSDM_GLC25_Extraction")
