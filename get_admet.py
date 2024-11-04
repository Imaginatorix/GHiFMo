import pandas as pd

# File location to .xlsx file
def get_admet(file_path='C:\A. Personal Files\ReSearch\A. Admet\selenium\merged_excel.xlsx'):
    # Load the Excel file
    try:
        df = pd.read_excel(file_path)
    except FileNotFoundError:
        return "File not found. Please check the file path and name."
    
    # Columns to extract  
    columns_to_extract = ['smiles' , 'Lipinski','PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
                          'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth']
    # Check if all columns exist in the dataset
    missing_columns = [col for col in columns_to_extract if col not in df.columns]
    if missing_columns:
        return f"The following columns are missing in the dataset: {missing_columns}"
    
    # Extract the relevant columns
    extracted_data = df[columns_to_extract]
    # No output yet. Not stored in .csv nor .xlsx file. 
   
    return extracted_data

# Call the function to test (file 'merged_files.xlsx' should be present in the directory)
# get_sa() # Uncomment this line to run it in a local environment

extracted_data = get_admet()
print(extracted_data)

# extracted_columns = get_admet()
# print(extracted_columns)

# When conversion is needed, just uncomment:
# extracted_columns.to_excel('extracted_data.xlsx', index=False)

