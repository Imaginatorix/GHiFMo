import os
import time
import glob
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


# def automated_admet(
#     smiles_file_path=r"C:\A. Personal Files\ReSearch\Final\download\smiles.xlsx",
#     chromedriver_path=r"C:\A. Personal Files\ReSearch\Final\chromedriver-win64\chromedriver.exe",
#     download_folder=r"C:\A. Personal Files\ReSearch\Final\download",
#     batch_size = 15):

# Refactored automated_admet function
def automated_admet(smiles_file_path, chromedriver_path, download_folder, batch_size):
    # Load SMILES data from the specified Excel file
    try:
        df = pd.read_excel(smiles_file_path)
        smiles_list = df['SMILES'].tolist()
    except (FileNotFoundError, KeyError) as e:
        print(f"Error loading SMILES file: {e}")
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
        # Open ADMET Lab website and process SMILES in batches
        driver.get("https://admetlab3.scbdd.com/server/screening")
        for i in range(0, len(smiles_list), batch_size):
            smiles_batch = smiles_list[i:i + batch_size]
            smiles_string = "\n".join(smiles_batch)

            # Interact with webpage elements and download results
            search = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "exampleFormControlTextarea1"))
            )
            search.clear()
            search.send_keys(smiles_string)
            time.sleep(3)

            submit_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "button.btn.btn-success"))
            )
            submit_button.click()

            # Download the results
            download_button = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.CLASS_NAME, "btn-outline-success"))
            )
            download_button.click()
            time.sleep(5)  # Ensure download completes before navigating away

    finally:
        driver.quit()

    print("Data extraction and download completed.")

automated_admet()
