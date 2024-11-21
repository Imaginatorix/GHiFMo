import os
import time
import glob
import pandas as pd

def automated_admet(
    smiles_list,
    chromedriver_path=r".\chromedriver-win64\chromedriver.exe",
    download_folder=r".\downloads",
    batch_size=1000):

    merged_output_path = os.path.join(download_folder, 'merged_output.csv')
    # Load the merged file to ensure data is accessible for column extraction
    merged_df = pd.read_csv(merged_output_path)

    # Define columns to extract
    columns_to_extract = ['Lipinski', 'PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
                          'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth']
    
    # Check if all columns exist in the merged dataset
    missing_columns = [col for col in columns_to_extract if not col in merged_df.columns]
    if missing_columns:
        print(f"The following columns are missing in the dataset: {missing_columns}")
        return None

    # Extract the relevant columns
    extracted_data = merged_df[columns_to_extract]
    print("Extracted Data:", extracted_data)
    return extracted_data
