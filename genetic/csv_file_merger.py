import os
import glob
import pandas as pd

def merge_csv_files(download_folder):
    # Merge downloaded CSV files
    csv_files = glob.glob(os.path.join(download_folder, '*.csv'))
    if not csv_files:
        print("No CSV files were downloaded.")
        return None

    merged_df = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)
    merged_output_path = os.path.join(download_folder, 'merged_output.csv')
    merged_df.to_csv(merged_output_path, index=False)
    print("Files have been merged into:", merged_output_path)
    
    # Clean up the downloaded CSV files to avoid redundancy
    for f in csv_files:
        os.remove(f)
    
    return merged_df
