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

def rank_cluster_drivers_by_pagerank(input_dir : str, clusters_file : str, cluster : float, k : int):
    """Aggregates PageRank scores of driver genes for given list of clusters and returns 
    top drivers for each cluster.

    Args:
        input_dir (str): Directory storing patient-specific PageRank results.
        clusters_file (str): File storing cluster relations for patients.
        cluster (float): Identifier of cluster to be analyzed.
        k (int): Number of highest ranked drivers to return.
    """
    # Extract patients belonging to given cluster.
    clusters_df = pd.read_csv(clusters_file, sep=',')
    subset_df = clusters_df[clusters_df['4_cluster_label']==cluster]

    patients = set(subset_df['Patient'])
    
    # Iterate over PageRank score files of all patients and aggregate PageRank scores.
    driver_scores = {}
    for pat in patients:
        pat_file = input_dir+f'{pat}_scores.csv'
        pat_df = pd.read_csv(pat_file, sep='\t')
        # Iterate over all driver entries.
        for driver, score in zip(pat_df['driver'], pat_df['score']):
            if not driver in driver_scores:
                driver_scores[driver]=score
            else:
                driver_scores[driver] = driver_scores[driver]+score
    
    sorted_pairs = sorted(driver_scores.items(), key=lambda x: x[1], reverse=True)
    sorted_drivers = [x[0] for x in sorted_pairs]
    return sorted_drivers[:k]
        
        
if __name__ == "__main__":

    """ truth_df_path = '../fdrs/PRAD/'
    false_df_path = '../fdrs/PRAD_075/'

    # Load truth and false DataFrames.
    truth_df = pd.read_csv(truth_df_path+'merged_results.csv', sep='\t')
    false_df = pd.read_csv(false_df_path+'merged_results.csv', sep='\t')

    fdr_df = fdr_patient(truth_df, false_df)
    fdr_df.to_csv(false_df_path+'patient_fdrs.csv', sep='\t') """
    
    # Input for cluster ranking according to PageRank results.
    input_dir = "pagerank_scores/"
    clusters_file = "../results/clustering/clusters/clusters_PARADYS_DysRegNet_PRAD.csv"
    cluster = 4.0
    k = 20
    driver_scores = rank_cluster_drivers_by_pagerank(input_dir, clusters_file, cluster, k)
    [print(x) for x in driver_scores]

