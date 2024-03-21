import random
import os
import json
import sklearn
import sklearn.cluster
import json
import numpy as np
import pandas as pd
from scipy import stats
from os.path import exists
# import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

def matrix_results_binary(d,t, patients, genes_M, dic_p_m):
    '''
    Creates binary matrix of patients and drivers
    
    Returns
    -------
    b : dataframe
        Binary matrix of patients and drivers
    patients_f : list
        patients that have at least one driver gene
    dic_patients_f : dict
        dictionary of patients and indices
    
    '''
    df = pd.DataFrame(index=patients, columns=genes_M, data=np.zeros([len(patients),len(genes_M)]))
    patients_f=list()
    indices=list()
    a=0
    for p in patients:
        data_p=data[data['patient']==p]
      
        
        sigs=data_p['driver']
        for i in sigs:
            indices.append(a)
            a+=1
        genes_s0 = set(sigs)
        p0 = p#p.replace(".", "-")
        if (p0) in dic_p_m.keys():
            genes_s = [value for value in dic_p_m[p0] if value in genes_s0]
        else:
            genes_s=[]
        for i in list(genes_s0):
            df.loc[p,i]=1
    d['id']=indices
    
    b = df.loc[(df).any(1), (df!=0).any(0)]   
    #b.to_excel('results '+cancer+'/matrix kappa 3 binary.xlsx', index=True)
    sumas=list()
    drivers=list()
    for col in b:
        a=b[col].sum()
        sumas.append(a)
        drivers.append(col)
    indices=np.argsort(np.array(sumas))[::-1][:t]
    top_drivers=[drivers[i] for i in indices]
    drivers_to_drop = list(set(drivers)-set(top_drivers))
    b_tops = b.drop(drivers_to_drop,axis=1)
    #b_tops.to_excel('results '+cancer+'/matrix kappa 3 binary top 20.xlsx', index=True)
    return b_tops, d, top_drivers



def matrix_results_binary_cross_validation(d, r, t):
    """
    Creates a binary matrix of patients and drivers for the r dataset

    Parameters
    ----------
    r : int
        index of cross-validation dataset

    Returns
    -------
    b : dataframe
        bnary matrix of patients and driver
    d : dataframe
        dictionary of patients and indices

    """
    patients_f = patients
    df = pd.DataFrame(index=patients_f, columns=genes_M, data=np.zeros([len(patients_f),len(genes_M)]))
    to_delete = set(random.sample(range(len(patients_f)), int(len(patients_f)*0.2)))
    patients_g = [x for i,x in enumerate(patients_f) if not i in to_delete]
    a=0
    indices=list()
    for p in patients:
        if p in patients_g:
           
            data_p=data[data['patient']==p]

            sigs=data_p['driver']
            for i in sigs:
                indices.append(a)
                a+=1
            genes_s0 = set(sigs)
            p0 = p#p.replace(".", "-")
            if (p0) in dic_p_m.keys():
                genes_s = [value for value in dic_p_m[p0] if value in genes_s0]
            else:
                genes_s=[]
            for i in list(genes_s0):
                df.loc[p,i]=1
                df.loc[p,i]=1
        else:
            indices.append(str(-1))
       
    d[r+1]=indices
    b = df.loc[(df).any(1), (df!=0).any(0)] 
    sumas=list()
    drivers=list()
    for col in b:
        a=b[col].sum()
        sumas.append(a)
        drivers.append(col)
    indices=np.argsort(np.array(sumas))[::-1][:t]
    top_drivers=[ drivers[i] for i in indices]
    drivers_to_drop = list(set(drivers)-set(top_drivers))
    b_tops = b.drop(drivers_to_drop,axis=1)
    return b_tops, d




def consensus_matrix(df,t, patients_f, CV):
    """
    Creates consensus matrix

    Parameters
    ----------
    df : DataFrame
        all data

    Returns
    -------
    CM : Numpy array
        Consensus Matrix
    df : DataFrame
        Updated data

    """
    CM = np.zeros([len(patients_f),len(patients_f)])
    for r in range(0,CV):
        print(r)
        b_r, df = matrix_results_binary_cross_validation(df, r, t)
        X = b_r.to_numpy()[:,1:]
        model = NMF(n_components=N, init='random', random_state=0)
        W = model.fit_transform(X)
        H = model.components_
        shape = W.shape
        #M = np.zeros([len(patients_f),len(patients_f)])
        
        for i in range(0,len(W)):
            for j in range(0,len(W)):
                max_i=W[i][0]
                max_j=W[j][0]
                for k in range(0, shape[1]):
                    if W[i][k]>=max_i:
                        max_i=W[i][k]
                        cluster_i = k
                    if W[j][k]>=max_j:
                        max_j=W[j][k]
                        cluster_j = k
                if cluster_i == cluster_j:
                    index_i = df.loc[df[r+1]==i].index.values[0]
                    id_i = df.loc[index_i,'id']
                    index_j = df.loc[df[r+1]==j].index.values[0]
                    id_j = df.loc[index_j,'id']
                    CM[id_i,id_j]+=1
    cols = list(range(1,CV+1))
    df['total']=df[cols].apply(lambda x: (x!='-1').sum(), axis=1)#df[cols].gt(str(-1)).sum(axis=1)
    for i in range(0,len(patients_f)):
        for j in range(0, len(patients_f)):
            if CM[i,j]!=0:
                sets=0
                index_i= df.loc[df['id']==i].index.values[0]
                index_j= df.loc[df['id']==j].index.values[0]
                for k in range(1,CV+1):
                    if df.loc[index_i,k] != '-1' and df.loc[index_j,k]!='-1':
                        sets+=1
                CM[i,j]/= sets
    return CM, df



def cluster_similarity_matrix(cons, N):
    """
    Aggregated Hierarchical Clustering

    Parameters
    ----------
    cons : Numpy Array
        Consensus Matrix

    Returns
    -------
    clusters: Numpy Array
        Array indicating the cluster of every patient

    """
    clusters = sklearn.cluster.AgglomerativeClustering(n_clusters=N).fit_predict(cons)
    return clusters



def top_t(t, path_results):
    b, df, top_drivers = matrix_results_binary(d,t)
    patients_f = df[df['id']!='-1']['patient']
    a, df = consensus_matrix(df, t, patients_f)
    clusters=cluster_similarity_matrix(a)
    to_drop = df.loc[df['id']=='-1'].index.values
    df_clean = df.drop(to_drop)
    df_clean['cluster']=clusters
    df_clean.to_csv(path_results)
    return df_clean


""" 
# Directory containing CSV files
directory = 'results BRCA/kappa = 3'

# Get a list of all CSV files in the directory
csv_files = [file for file in os.listdir(directory) if file.endswith('.csv')]

# Initialize an empty DataFrame to store the merged data
merged_df = pd.DataFrame()

# Loop through each CSV file and concatenate it to the merged DataFrame
for csv_file in csv_files:
    file_path = os.path.join(directory, csv_file)
    df = pd.read_csv(file_path)
    
    # Add a new column named 'patient' with the file name 

    patient_name = os.path.splitext(os.path.splitext(csv_file)[0])[0]
    df['patient'] = patient_name[:-10]
    # Concatenate the DataFrame to the merged DataFrame
    merged_df = pd.concat([merged_df, df], ignore_index=True)
     """
    
os.chdir('../')
cancer = 'PRAD'
# Read merged driver-patient results. Needs to contain columns 'patient' and 'driver'.
print("Reading results data...")
merged_df_path = f'fdrs/{cancer}/merged_results.csv'
data=pd.read_csv(merged_df_path, sep='\t')

# Loading per-patient mutation data.
with open(f'filtered_data/{cancer}/dic_patients_mutations.json') as i:
    dic_p_m = json.load(i)

# Output path.
path_to_save = f'clustering/{cancer}/'

# Number of top drivers to consider for clustering.
t = 30

# Name of considered cohort.
cancer = 'PRAD'

CV=100

for N in [2,3,4]:
    print(f"Starting clustering with N={N}...")
    patients=set(data['patient'])
    patients=list(patients)
    drivers=set(data['driver'])
    drivers=list(drivers)
    genes_M = drivers 
    # data={'patients':patients}
    d = pd.DataFrame(data)
    print("Creating binary matrix from results...")
    b, df, top_drivers = matrix_results_binary(d,t, patients, drivers, dic_p_m)
    patients_f = df[df['id']!='-1']['patient']
    print("Building consensus matrix...")
    a, df = consensus_matrix(df, t, patients_f, CV)
    print("Clustering similarity matrix...")
    clusters=cluster_similarity_matrix(a, N)
    to_drop = df.loc[df['id']=='-1'].index.values
    df_clean = df.drop(to_drop)
    df_clean2 = pd.DataFrame()
    data['cluster']=clusters
    data.to_csv(path_to_save+'clusters'+cancer+'_top'+str(t)+'N='+str(N)+',CV='+str(CV)+'.csv', sep='\t')
