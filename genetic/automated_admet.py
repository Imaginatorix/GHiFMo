import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.keys import Keys
import time
import os

# Set Chrome options
options = webdriver.ChromeOptions()
options.add_experimental_option("detach", True)

# Initialize WebDriver
PATH = r"C:\A. Personal Files\ReSearch\selenium\chromedriver-win64\chromedriver.exe"
service = Service(PATH)
driver = webdriver.Chrome(service=service)

# Open Google
driver.get("https://www.google.com")

# Locate the search input field using the 'name' attribute, which is more reliable
search = driver.find_element(By.NAME, "q")

# Send search query to the input field
search.send_keys("Stephen Curry")
search.send_keys(Keys.RETURN)  # Press 'Enter' to submit the search

# Wait for a while before navigating back
time.sleep(10)
driver.back()

time.sleep(60)

# Close the driver (optional)
driver.quit()
