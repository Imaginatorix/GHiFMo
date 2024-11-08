import pandas as pd

def extract_columns(merged_df, columns_to_extract):
    # Check if all columns exist in the merged dataset
    missing_columns = [col for col in columns_to_extract if col not in merged_df.columns]
    if missing_columns:
        print(f"The following columns are missing in the dataset: {missing_columns}")
        return None

    # Extract the relevant columns
    extracted_data = merged_df[columns_to_extract].values
    return extracted_data
