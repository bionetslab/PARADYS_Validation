import os
import pandas as pd

# Directory containing CSV files
os.chdir('../')

cohort = 'COAD_Decoy_075'
directory = f'results/{cohort}/'

# List of all contained CSV files.
csv_files = [file for file in os.listdir(directory) if file.endswith('.csv')]

# Resulting output DataFrame.
merged_df = pd.DataFrame()


for csv_file in csv_files:
    file_path = os.path.join(directory, csv_file)
    df = pd.read_csv(file_path)
    
    # Add a new column named 'patient' with the file name 
    patient_name = os.path.splitext(os.path.splitext(csv_file)[0])[0]
    df['patient'] = patient_name[:12]
    
    # Renam 'drivers' column to 'driver' column.
    df.rename(columns={'drivers':'driver'}, inplace=True)

    # Merge into common DataFrame.
    merged_df = pd.concat([merged_df, df], ignore_index=True)

# Save merged DataFrame.
merged_df.to_csv(f'{directory}/merged_results.csv', sep='\t')
