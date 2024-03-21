# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:37:26 2024

@author: sangoibe
"""

import os
import json
import sklearn.preprocessing
from collections import defaultdict
import pandas as pd
import mygene

os.chdir("../")

# We assume input data to be stored in directory ../filtere_data/<CANCER_TYPE>/.

# Data format specifications (only concerning mutation preprocessing)
effect = 'effect'
gene = 'gene'
sample = 'Sample_ID'
#sample = 'sample'

# Path to unprocessed mutation data.
path = 'data PRAD/TCGA-PRAD.muse_snv.tsv'
#path = r'data/M_matrix.csv'

# Cancer cohort to be analyzed.
cancer='PRAD'

# Path of dysregulation network feather files.
network_path = f'data/{cancer}/{cancer}_nets_DysRegNet_no_dir.fea'

# Output path for filtered data.
output_path = f'filtered_data/{cancer}'

# Read in mutation data.
# data=pd.read_csv(path, sep='\t')
# Commented: data=pd.read_csv(path, sep='\t')
# cnvs=pd.read_csv(r'data/TCGA-PRAD_CNVS.tsv', sep='\t')
# Commented: cnvs=pd.read_csv(r'data PRAD/TCGA-PRAD_CNVS.tsv', sep='\t')

# Load dysregulation network.
print("Loading network feather file...")
networks = pd.read_feather(network_path)

# Make sure that output directory exists.
if not os.path.exists(output_path):
    os.makedirs(output_path)

## Uncomment if networks aren't binarized
# networks = binarize(networks)

def create_dictionary_genes_ids(cnvs):
    '''
    Transltaion of gene symbols to gene ids
    Stores a dictionary in the outputpath 

    Parameters
    ----------
    cnvs : dataFrame
        CNV original files with column named Gene Symbols

    Returns
    -------
    dic_genes_ids : dict
        Keys are gene ids and values are gene symbols

    '''
    mg = mygene.MyGeneInfo()
    geneList = cnvs['Gene Symbol']
    dic_genes_ids=dict()
    # with open(f'{output_path}/dic_ensembl_genes.json') as i:
    #     dic_genes_ids=json.load(i)
    for a in geneList:
        d = a.partition('.')[0]
        b = mg.querymany(d, scopes='ensembl.gene', fields='symbol', species='human', returnall=True)
        
        c=b['out']
        if 'symbol' in c[0].keys():
        #print(c[0]['symbol'])
            dic_genes_ids[a]=c[0]['symbol']
    a_file = open(f'{output_path}/dic_ensembl_genes.json', "w")
    json.dump(dic_genes_ids, a_file)
    a_file.close()
    return dic_genes_ids





def binarize(networks):
    '''
    Creates a binarized edge matrix given threshold for edge weights.
    To be used only in case that original dysregulation networks are not binarized

            Parameters:
                    networks (DataFrame): original networks

            Returns:
                    binarized_data (DataFrame): binarized edge matrix
    '''
    samples = list()
    edges = list()
    d=networks.iloc[: , 1:].to_numpy()
    #d=sklearn.preprocessing.binarize(d,threshold=0.1)
    d=networks.iloc[: , :-1].to_numpy()
    for col in networks.columns:
        edges.append(col)
    for index, row in networks.iterrows():
        samples.append(networks.loc[index, 'patient id'])
    edges.pop(0)
    binarized_data=pd.DataFrame(d, index=samples, columns=edges)
    return binarized_data


def set_up_mutations(effect, gene, sample, cancer):
    '''
    Aggregates cnvs and snvs and stores different files to be used later in the analysis

    Parameters
    ----------
    effect : String
        Name of the column of the original mutation file that specifies the effect of every mutation
    gene : String
        Name of the column of the original mutation file that specifies the gene ID
    sample : String
        Name of the column of the original mutation file that specifies the sample ID
    cancer : String
        Name of the cohort

    Returns
    -------
    None.
    It stores:
        dic_patients_mutations --> dictionary, keys are patients and values are a list of mutated genes for the given patient
        dic_mutations_patients --> dictionary, keys are genes and patients are a list of patients that have the given gene mutated
        genes_M --> list of genes that were mutated at least for one patient
        patients_M --> list of patients that were present in the original mutation data

    '''

    type_to_filter_out = 'synonymous_variant'
    # read the data
    genes = []
    patients = []
    for index, row in data.iterrows():   
        if data.loc[index,effect]!=type_to_filter_out:
            genes.append(data.loc[index,gene])
            patients.append(data.loc[index,sample][:-1])
    
    # add also cnvs
    if cancer == 'BRCA':
        dic_genes_ids = create_dictionary_genes_ids(cnvs)
        patients_cnvs = list(cnvs.columns)[1:-1]
        for index, row in cnvs.iterrows():
            for p in patients_cnvs:
                if cnvs.loc[index,p]!=0:
                    if cnvs.loc[index,'Gene Symbol'] in dic_genes_ids.keys():
                        genes.append(dic_genes_ids[cnvs.loc[index,'Gene Symbol']])
                        patients.append(p[:-1])

                        
    else:
        patients_cnvs = list(cnvs.columns)[1:]
        for index, row in cnvs.iterrows():
            for p in patients_cnvs:
                if cnvs.loc[index,p]!=0:
                    genes.append(cnvs.loc[index,'Gene Symbol'])
                    patients.append(p)
                    
    # create dictionaries               
    dic_patients_mutations = defaultdict(list)      
    dic_mutations_patients = defaultdict(list) 
    for i,j in zip(genes,patients):
        dic_mutations_patients[i].append(j)
        dic_patients_mutations[j].append(i)
    
    print(len(patients))
    print(len(genes))
    patients=list(set(patients))
    genes=list(set(genes))
    
    a_file = open(f'{output_path}/dic_patients_mutations.json', "w")
    json.dump(dic_patients_mutations, a_file)
    a_file.close()
    a_file = open(f'{output_path}/genes_M.json', "w")
    json.dump(genes, a_file)
    a_file.close()
    a_file = open(f'{output_path}/dic_mutations_patients.json', "w")
    json.dump(dic_mutations_patients, a_file)
    a_file.close()
    a_file = open(f'{output_path}/patients_M.json', "w")
    json.dump(patients, a_file)
    a_file.close() 


def set_up_networks(binarized_networks):
    '''
    Creates preprocessed files of the networks to be used in the analysis

    Parameters
    ----------
    binarized_networks : DataFrame
        Indices are patients and columns are edges (strings of pairs of genes separated by coma and inside ())
        Cell values are 0 or 1 if the given edge was present in the patient network

    Returns
    -------
    None.
    Stores:
        genes_networks --> all genes that were present in the networks (at least once)
        patients --> all patients that have network and mutation data (intersection with patients_M)
        dic_patients_edges --> dictionary, keys are patients and values are edges present in the network (as '(gene1,gene2)')
        dic_patients_graphs --> dictionary, keys are patients and values are dictionaries that have tfs as keys and genes connected to them as values
        edges_counts --> dictionary, keys are edges and values are the number of patients whose networks contain the given edge (to be used only for analysing frequent edges instead of all)
    '''
    # check that we have mutation info from every patient of the network original data, if not remove the patient
    # Load patient data
    with open(f'{output_path}/patients_M.json') as j:
        patients = set(json.load(j))

    
    # patients_nets = set(binarized_networks.index)
    # comment previous lines and uncomment next 4 lines--> uncomment if patients are a column and not indices
    patients_nets=list(binarized_networks['patient id'])
    patients_nets=set(patients_nets)
    patients_common = patients.intersection(patients_nets)
    print("Lenght of patient intersection: ", len(patients_common))
    binarized_networks.set_index('patient id', inplace=True)

    
    columns = binarized_networks.columns # [:-1] --> uncomment if patients are a column and not indices
    

    genes_networks = []
    edges = []
    all_dic_edges = []
    
    for p in patients_common:
        edges_p = [col for col in columns if binarized_networks.loc[p, col] != 0]
        edges.append(edges_p)
        
        sources, targets = [], []
        for i in edges_p:
            gene0, gene1 = i.split(',')
            sources.append(gene0.strip("() '"))
            targets.append(gene1.strip("() '"))
        
        genes_networks.extend(sources)
        genes_networks.extend(targets)
        
        dict_edges_p = defaultdict(list)
        for source, target in zip(sources, targets):
            dict_edges_p[source].append(target)
        all_dic_edges.append(dict_edges_p)
            
    
   
    
    genes_networks = list(set(genes_networks))
    a_file = open(f'{output_path}/genes_networks.json', "w")
    json.dump(genes_networks, a_file)
    a_file.close()
    a_file = open(f'{output_path}/patients.json', "w")
    json.dump(list(patients_common), a_file)
    a_file.close()
    zip_iterator = zip(patients_common, edges)
    dic_patients_edges = dict(zip_iterator)
    b_file = open(f'{output_path}/dic_patients_edges.json', "w")
    json.dump(dic_patients_edges, b_file)
    b_file.close()
    zip_iterator = zip(patients_common, all_dic_edges)
    dic_patients_graphs = dict(zip_iterator)
    b_file = open(f'{output_path}/dic_patients_graphs.json', "w")
    json.dump(dic_patients_graphs, b_file)
    b_file.close()
    
  
    counts=[]
    edges_totals = binarized_networks.columns
    for col in binarized_networks.columns:
        counts.append(binarized_networks[col].sum())
    zip_iterator = zip(edges_totals, counts)
    edges_counts = dict(zip_iterator)
    b_file = open(f'{output_path}/edges_counts.json', "w")
    json.dump(edges_counts, b_file)
    b_file.close()

#print("Preprocessing mutations...")
#set_up_mutations(effect,gene, sample, cancer=cancer)
print("Preprocessing networks...")
set_up_networks(networks)

