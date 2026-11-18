import glob, os

# --- Inspect the actual files on disk ---
for root in ["/data1/lorraine-20419/benchmark_cmip6/historical",
             "/data1/lorraine-20419/ocean_data/ERSSTv5",
             "/data1/lorraine-20419/ocean_data/CPC"]:
    print("="*80)
    print(root)
    for p in sorted(glob.glob(os.path.join(root, "**", "*"), recursive=True)):
        if os.path.isfile(p):
            print("  ", p)
