import os
import time
import glob
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

def automated_admet(
    smiles_file_path=r"C:\A. Personal Files\ReSearch\Final\download\smiles.xlsx",
    chromedriver_path=r"C:\A. Personal Files\ReSearch\Final\chromedriver-win64\chromedriver.exe",
    download_folder=r"C:\A. Personal Files\ReSearch\Final\download",
    batch_size=15):

    try:
        df = pd.read_excel(smiles_file_path)
        smiles_list = df['SMILES'].tolist()
    except FileNotFoundError:
        print("SMILES file not found. Please check the file path.")
        return None
    except KeyError:
        print("The specified 'SMILES' column was not found in the file.")
        return None

    options = webdriver.ChromeOptions()
    options.add_experimental_option("detach", True)
    prefs = {"download.default_directory": download_folder}
    options.add_experimental_option("prefs", prefs)

    service = Service(chromedriver_path)
    driver = webdriver.Chrome(service=service, options=options)
    
    try:
        driver.get("https://admet.ai.greenstonebio.com/")

        for i in range(0, len(smiles_list), batch_size):
            smiles_batch = smiles_list[i:i + batch_size]
            smiles_string = "\n".join(smiles_batch)

            search = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "text-smiles-input"))
            )
            search.clear()
            search.send_keys(smiles_string)
            time.sleep(3)

            submit_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.ID, "predict-button"))
            )
            submit_button.click()

            download_button = WebDriverWait(driver, 10).until(
                # EC.element_to_be_clickable((By.CSS_SELECTOR, "button.btn.btn-primary"))
                EC.element_to_be_clickable((By.XPATH, '//button[text()="Download Results"]'))
            )
            download_button.click()

            time.sleep(5)  # Wait for download to complete (improve logic here if necessary)
            driver.get("https://admet.ai.greenstonebio.com")

    finally:
        driver.quit()

    # Get all CSV files in the download folder
    csv_files = glob.glob(os.path.join(download_folder, '*.csv'))

    # Check if no CSV files are found
    if not csv_files:
        print("No CSV files were downloaded.")
        return None

    # Check if there is only one CSV file
    if len(csv_files) == 1:
        print("Only one CSV file found. Skipping merging.")
        merged_df = pd.read_csv(csv_files[0])  # Directly load the single CSV file
    else:
        # If multiple files exist, merge them
        merged_df = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)
        merged_output_path = os.path.join(download_folder, 'merged_output.csv')
        merged_df.to_csv(merged_output_path, index=False)
        print("Files have been merged into:", merged_output_path)

        # Remove the original CSV files after merging
        for f in csv_files:
            os.remove(f)

    # Define the columns to extract (for ADMET-Ai only)
    columns_to_extract = [
        'Lipinski', 'PPBR_AZ', 'VDss_Lombardo', 'CYP3A4_Veith', 'CYP3A4_Substrate_CarbonMangels',
        'CYP2D6_Veith', 'CYP2D6_Substrate_CarbonMangels', 'Clearance_Hepatocyte_AZ', 'Clearance_Microsome_AZ',
        'Half_Life_Obach', 'DILI', 'hERG'
    ]

    # Check if all columns exist in the dataset
    missing_columns = [col for col in columns_to_extract if col not in merged_df.columns]
    if missing_columns:
        print(f"The following columns are missing in the dataset: {missing_columns}")
        return None

    # Extract the relevant data
    extracted_data = merged_df[columns_to_extract].values
    return extracted_data
