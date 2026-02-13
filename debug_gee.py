import ee
try:
    project_id = "your-project-id"
    print(f"Attempting to initialize GEE with project: {project_id}")
    ee.Initialize(project=project_id)
    print("Success! GEE initialized.")
    # Try a simple operation to confirm access
    print(f"Default root: {ee.data.getAssetRoot()}")
except Exception as e:
    print(f"Failed to initialize: {e}")
