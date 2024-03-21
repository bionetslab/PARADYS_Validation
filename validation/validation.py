import pandas as pd


def fdr_patient(truth_df,false_df):
    '''
    Calculates FDR for every patient

    Parameters
    ----------
    truth_df : DataFrame
        Original results, columns are patient and driver
    false_df : DataFrame
        Decoy mutations results

    Returns
    -------
    fdr_df : DataFrame
        Columns are patient and FDR

    '''
    
    # Create an empty dictionary to store FDR values for each patient
    fdr_per_patient = {}
    
    # Get the list of unique patients
    patients = truth_df['patient'].unique()
    
    for patient in patients:
        # Filter truth_df and false_df for the current patient
        truth_patient_df = truth_df[truth_df['patient'] == patient]
        false_patient_df = false_df[false_df['patient'] == patient]
    
        # Remove duplicate driver entries.
        truth_patient_df = truth_patient_df.drop_duplicates(subset=['driver'])
        false_patient_df = false_patient_df.drop_duplicates(subset=['driver'])

        # Calculate True Positives (TP)
        tp = len(truth_patient_df[(truth_patient_df['driver'].isin(false_patient_df['driver']))])
    
        # Calculate False Positives (FP)
        fp = len(false_patient_df) - tp
    
        # Calculate False Discovery Rate (FDR) for the current patient
        if (tp + fp) > 0:
            fdr = fp / (tp + fp)
        else:
            fdr = 0  # Handle the case where there are no positive predictions
    
        # Store the FDR value for the current patient in the dictionary
        fdr_per_patient[patient] = fdr
    
    # Print FDR values for each patient
    for patient, fdr in fdr_per_patient.items():
        print(f"Patient {patient}: FDR = {fdr}")
    fdr_df = pd.DataFrame(list(fdr_per_patient.items()), columns=['patient', 'FDR'])
    return fdr_df


if __name__ == "__main__":

    truth_df_path = '../fdrs/PRAD/'
    false_df_path = '../fdrs/PRAD_075/'

    # Load truth and false DataFrames.
    truth_df = pd.read_csv(truth_df_path+'merged_results.csv', sep='\t')
    false_df = pd.read_csv(false_df_path+'merged_results.csv', sep='\t')

    fdr_df = fdr_patient(truth_df, false_df)
    fdr_df.to_csv(false_df_path+'patient_fdrs.csv', sep='\t')

