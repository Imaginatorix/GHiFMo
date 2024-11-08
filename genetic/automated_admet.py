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
    batch_size = 15):
    # Load SMILES data from the specified Excel file
    try:
        df = pd.read_excel(smiles_file_path)
        smiles_list = df['SMILES'].tolist()
    except FileNotFoundError:
        print("SMILES file not found. Please check the file path.")
        return None
    except KeyError:
        print("The specified 'SMILES' column was not found in the file.")
        return None

    # Set up Selenium options
    options = webdriver.ChromeOptions()
    options.add_experimental_option("detach", True)
    prefs = {"download.default_directory": download_folder}
    options.add_experimental_option("prefs", prefs)

    # Initialize the WebDriver
    service = Service(chromedriver_path)
    driver = webdriver.Chrome(service=service, options=options)
    
    try:
        # Open ADMET Lab website
        driver.get("https://admetlab3.scbdd.com/server/screening")
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "profile-tab"))).click()

        # Process SMILES in batches
        for i in range(0, len(smiles_list), batch_size):
            smiles_batch = smiles_list[i:i + batch_size]
            smiles_string = "\n".join(smiles_batch)

            # Enter SMILES data into the input field
            search = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "exampleFormControlTextarea1"))
            )
            search.clear()
            search.send_keys(smiles_string)
            time.sleep(3)

            # Submit the batch for processing
            submit_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "button.btn.btn-success"))
            )
            submit_button.click()

            # Wait for results to load and download the results
            WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "text-center")))
            download_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CLASS_NAME, "btn-outline-success"))
            )
            download_button.click()

            # Wait for download to complete
            time.sleep(5)
            driver.get("https://admetlab3.scbdd.com/server/screening")
            WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "profile-tab"))).click()

    finally:
        # Close the driver after all downloads are complete
        driver.quit()

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

    # Define columns to extract
    columns_to_extract = ['Lipinski', 'PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
                          'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth']
    
    # Check if all columns exist in the merged dataset
    missing_columns = [col for col in columns_to_extract if col not in merged_df.columns]
    if missing_columns:
        print(f"The following columns are missing in the dataset: {missing_columns}")
        return None

    # Extract the relevant columns
    extracted_data = merged_df[columns_to_extract].values
    return extracted_data

automated_admet()
