import pandas as pd

def get_admet(file_path='merged_files.xlsx'):
    # Load the Excel file
    try:
        df = pd.read_excel(file_path)
    except FileNotFoundError:
        return "File not found. Please check the file path and name."
    
    # Define the columns to extract
    columns_to_extract = ['Caco2', 'HIA', 'PPB', 'logVDss', 'CYP3A4-inh', 'CYP3A4-sub', 
                          'CYP2D6-inh', 'CYP2D6-sub', 'cl-plasma', 't0.5', 'DILI', 'hERG', 'Synth']
    
    # Check if all columns exist in the dataset
    missing_columns = [col for col in columns_to_extract if col not in df.columns]
    if missing_columns:
        return f"The following columns are missing in the dataset: {missing_columns}"
    
    # Extract the relevant columns
    extracted_data = df[columns_to_extract]
    # No output yet. Not stored in .csv nor .xlsx file. 
    # When conversion is needed, just uncomment:
    #extracted_columns.to_excel('extracted_data.xlsx', index=False)
    return extracted_data
get_sa()
