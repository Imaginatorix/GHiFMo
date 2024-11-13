import os
import time
import glob
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

def automated_admet(
    smiles_list,
    chromedriver_path=r".\chromedriver-win64\chromedriver.exe",
    download_folder=r".\downloads",
    batch_size=1000):

    # Set up Selenium options
    options = webdriver.ChromeOptions()
    options.add_experimental_option("detach", True)
    prefs = {"download.default_directory": os.path.abspath(download_folder)}
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
                EC.element_to_be_clickable((By.XPATH, '//button[@aria-label="Submit"]'))
            )
            submit_button.click()

            download_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "button.btn.btn-outline-success"))
            )
            download_button.click()

            # Wait for download to complete
            time.sleep(5)
            driver.get("https://admetlab3.scbdd.com/server/screening")
            WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "profile-tab"))).click()

    finally:
        # Close the driver after all downloads are complete
        driver.quit()

    # Find all CSV files in the download folder
    csv_files = glob.glob(os.path.join(download_folder, '*.csv'))
    csv_files = [f for f in csv_files if not f.endswith("merged_output.csv")]
    if not csv_files:
        print("No CSV files were found.")
        return None
    
    # Merge all CSV files
    merged_df = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)
    merged_output_path = os.path.join(download_folder, 'merged_output.csv')
    merged_df.to_csv(merged_output_path, index=False)
    print("Files have been merged into:", merged_output_path)

    # Clean up the downloaded CSV files to avoid redundancy
    for f in csv_files:
        os.remove(f)

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
