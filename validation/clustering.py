import random
import os
import json
import sklearn
import sklearn.cluster
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF




CV =100
output_path = 'clusters'
input_path = 'results'

'''
FUNCTIONS:
    matrix_results_binary: 
matrix_results_binary_cross_validation: Constructs a binary matrix representing the presence or absence of specific genes in patients based on driver predictions excluding a subset of patients for each iteration.
    consensus_matrix: Constructs a consensus matrix using Non-negative Matrix Factorization (NMF) to identify patterns across multiple iterations of cross-validation.
cluster_similarity_matrix: Performs hierarchical clustering on the consensus matrix to determine patient clusters.
'''




def matrix_results_binary_cross_validation(d, r, patients, drivers, input_data):
    '''
    Creates a binary matrix of patients and drivers for the r dataset (excluding a subset of 20% of the patients)

    Parameters
    ----------
    d : dataframe
        DESCRIPTION.
    r : int
        iteration id
    patients : list
    drivers : list
        driver genes
    input_data : dataframe
        Driver predictions

    Returns
    -------
    b_tops : dataframe
        binary matrix with patients and top t drivers, filled with zeros and one if a gene is driver for a given patient
    d : dataframe
        stores for every patient if it was selected or not (-1) for each iteration  and the number of times it was selected

    '''
    patients_f = patients
    df = pd.DataFrame(index=patients, columns=drivers, data=np.zeros([len(patients),len(drivers)]))
    to_delete = set(random.sample(range(len(patients_f)), int(len(patients_f)*0.2)))
    patients_g = [x for i,x in enumerate(patients_f) if not i in to_delete]
    a=0
    indices=list()
    for p in patients:
        if p in patients_g:
            indices.append(a)
            a+=1
            data_p=input_data[input_data['patient']==p]

            sigs=data_p['driver']
            genes_s= set(sigs)
            p0 = p
            for i in genes_s:
                df.loc[p,i]=1
        else:
            indices.append(str(-1))
       
    d[r+1]=indices
    b = df.loc[df.any(axis=1), (df != 0).any(axis=0)]

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




def consensus_matrix(df, patients,drivers, input_data):
    """
    Creates consensus matrix

    Parameters
    ----------
    df : DataFrame
        column patients
    patients : list
    drivers : list
        driver genes
    input_data : dataframe
        Driver predictions

    Returns
    -------
    CM : Numpy array
        Consensus Matrix
    df : DataFrame
        Updated data

    """
    CM = np.zeros([len(patients),len(patients)])
    for r in range(0,CV):
        b_r, df = matrix_results_binary_cross_validation(df, r, patients, drivers, input_data)
        X = b_r.to_numpy()[:,1:]
        model = NMF(n_components=N, init='random', random_state=0)
        W = model.fit_transform(X)
        H = model.components_
        shape = W.shape

        
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
    for i in range(0,len(patients)):
        for j in range(0, len(patients)):
            if CM[i,j]!=0:
                sets=0
                index_i= df.loc[df['id']==i].index.values[0]
                index_j= df.loc[df['id']==j].index.values[0]
                for k in range(1,CV+1):
                    if df.loc[index_i,k] != '-1' and df.loc[index_j,k]!='-1':
                        sets+=1
                CM[i,j]/= sets
    return CM, df



def cluster_similarity_matrix(cons,N):
    """
    Aggregated Hierarchical Clustering

    Parameters
    ----------
    cons : Numpy Array
        Consensus Matrix
    
    N : int
        Number of clusters

    Returns
    -------
    clusters: Numpy Array
        Array indicating the cluster of every patient

    """
    clusters = sklearn.cluster.AgglomerativeClustering(n_clusters=N).fit_predict(cons)
    return clusters



cancer = 'PRAD'
cluster_df = pd.DataFrame()
for N in [2,3,4]:
    print(N)
    t=30
    input_data = pd.read_csv(input_path+'/PRAD/merged_results_PRAD_DysRegNet_NoDir.csv',sep='\t')
    patients=set(input_data['patient'])
    patients=list(patients)
    drivers=set(input_data['driver'])
    drivers=list(drivers)


    with open('filtered data/PRAD/dic_patients_mutations.json') as i:
        dic_p_m = json.load(i)
    patients_id = [patient[:12] for patient in patients]
    
    data = {'patients': patients_id,'id':list(range(len(patients_id)))}
    d = pd.DataFrame(data)
    a, df = consensus_matrix(d, t, patients, drivers, input_data)
    clusters=cluster_similarity_matrix(a,N)
    clusters2 = [str(N) + '.' + str(x) for x in clusters]
    cluster_df['cluster '+str(N)] = clusters2
    
cluster_df['patient']=patients_id
cluster_df.set_index('patient', inplace=True)
#cluster_df.to_csv(output_path+'/'+cancer+'_top '+str(t)+'_DysRegNet_NoDir.csv')