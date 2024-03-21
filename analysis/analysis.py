import json
import numpy as np
import pandas as pd
from scipy import stats
from os.path import exists
import multiprocessing
from multiprocessing import Pool
import itertools
import os
from scipy.stats import chisquare, chi2
import networkx as nx

cancer = 'PRAD_Decoy_075' 
kappa = 3
scores = False

os.chdir("../")



print("Loading preprocessed data...")

with open('filtered_data/'+cancer+'/patients.json') as i:
    patients = json.load(i)

print("Number of patients:", len(patients))

with open('filtered_data/'+cancer+'/genes_M.json') as i:
    genes_M = json.load(i)

with open('filtered_data/'+cancer+'/dic_patients_edges.json') as i:
    dic_patients_edges = json.load(i)

with open('filtered_data/'+cancer+'/dic_mutations_patients.json') as i:
    dic_mutations_patients = json.load(i)

with open('filtered_data/'+cancer+'/dic_patients_mutations.json') as i:
    dic_patients_mutations = json.load(i)   
    


def calculate_page_rank(data, d=0.85, directed=True):
    '''
    Creates a linegraph and applies pagerank algorithm to get driver scores

    Parameters
    ----------
    data : Dataframe
        Results data, columns are patient, drivers, dysregulations and p-values
    d : float, optional
        Dumping factor for the pagerank algorithm. The default is 0.85.
    directed : Bool, optional
        True if original dysregulation network is directed, false otherwise. The default is True.

    Returns
    -------
    data : DataFrame
        
    Same as input, with additional column with scores and ordered in descending score

    '''
    # Extract unique drivers and edges
    total_drivers = set(data['drivers'])
    target_edges = set(data['dysregulations'])
    
    # Create a directed or undirected graph
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    # Add vertices
    G.add_nodes_from(total_drivers.union(target_edges))

    # Assign attributes to vertices (genes and edges)
    vertex_weights = {}
    vertex_names = {}
    vertex_colors = {}
    edge_weights = {}

    for vertex in total_drivers.union(target_edges):
        vertex_weights[vertex] = data[data['drivers'] == vertex]['dysregulations'].count() / len(data['patient'])
        vertex_names[vertex] = vertex
        vertex_colors[vertex] = 0.5 if vertex in target_edges else 1

    nx.set_node_attributes(G, vertex_weights, 'weight')

    # Preprocess edges to remove unwanted characters and create a dictionary for quick access
    processed_edges = {edge.replace("(", "").replace(")", "").replace("'", "").replace(" ", ""): edge for edge in target_edges}
    
    # Add edges between connected edges with weight 1
    for edge in target_edges:
        edge1 = edge.partition(', ')[2].replace("(", "").replace(")", "").replace("'", "").replace(" ", "")
        for edge2 in target_edges:
            edge0 = edge2.partition(', ')[0].replace("(", "").replace("'", "").replace(" ", "")
            if edge0 == edge1:
                G.add_edge(edge, processed_edges[edge2], weight=1)
    
    


    # Add edges between genes and edges with weight -log(p-value)
    for driver in total_drivers:
        data_driver = data[data['drivers'] == driver]
        for edge in data_driver['dysregulations']:
            G.add_edge(edge, driver, weight=-np.log(data_driver[data_driver['dysregulations'] == edge]['p-values'].values[0]))

    # Apply PageRank algorithm and store scores
    pr = nx.pagerank(G, alpha=d, personalization=vertex_weights, weight='weight')
    # Add scores as a new column in the original data DataFrame
    data['PageRank_Score'] = data['drivers'].map(pr)

    # Sort the DataFrame based on the scores
    data = data.sort_values(by='PageRank_Score', ascending=False)

    return data


 

def get_relations(kappa, cancer, patients, dic_patients_mutations):
    '''
    Returns and stores all connections between mutations and dysregulations

    Parameters
    ----------
    kappa : int
        Maximum path length between edge and mutated gene
    cancer : str
        Cancer type
    patients : list
        List of patient IDs
    dic_patients_mutations : dict
        Dictionary containing mutated genes for each patient

    Returns
    -------
    data : DataFrame
        For every mutated gene for every patient, list of dysregulations connected to it in a max distance kappa
        Columns are genes (mutations), patients (patient ids) and connected dysregulations (genes connected to the mutation)

    '''
    with open('filtered_data/'+cancer+'/dic_patients_graphs.json') as i:
        dic_patients_graphs = json.load(i)

    gs=[] # mutations
    ps=[] # patients
    ds=[] # connected genes
    for p in patients:
        genes_p = dic_patients_mutations.get(p, [])
        edges = dic_patients_graphs.get(p, {})

        for g in genes_p:
            if g in edges:
                connected_dysregulations = set()
                genes_stack = [g]

                for _ in range(kappa):
                    next_genes = set()

                    for gene in genes_stack:
                        if gene in edges:
                            next_genes.update(set(edges[gene]) - connected_dysregulations)

                    connected_dysregulations.update(next_genes)
                    genes_stack = list(next_genes)

                if connected_dysregulations:
                    for dysregulation in connected_dysregulations:
                        gs.append(g)
                        ps.append(p)
                        ds.append(dysregulation)
    data = pd.DataFrame({
        'connected dysregulations': ds,
        'genes': gs,
        'patients': ps
    })
    data.to_csv('filtered_data/'+cancer+'/gene_sets_kappa' + str(kappa) + '.csv', index=False)
    return data




def mutation_condition (gene,patient):
    '''
    Mutation condition

    Parameters
    ----------
    gene : string or list
    patient : string

    Returns
    -------
    bool
        True if gene mutated for that given patient, False otherwise

    '''
    patient = patient
    if type(gene)== list:
        gene = gene[0]
    if gene in genes_M:
        if patient in dic_mutations_patients[gene]: 
            return True
        else:
            return False
    else:
        return False
        



def dysregulation_condition(edge, patient):
    '''
    Dysregulation condition

    Parameters
    ----------
    edge : string
    patient : string

    Returns
    -------
    bool
        True if edge is present in the patient's dysregulatory network, False otherwise

    '''
    if edge in dic_patients_edges[patient]:
        return True
    else:
        return False


def obs(p, edge_string, gene):
    '''
    Create observation matrix for Chi2 test

    Parameters
    ----------
    p : string
        Patient.
    edge_string : string
        Edge.
    gene : string

    Returns
    -------
    o : DataFrame
        2x2 Observation matrix.

    '''
    o = np.zeros((2, 2), dtype=int)
    dis = dysregulation_condition(edge_string, p)
    mut = mutation_condition(gene,p)
    if dis and mut:
        o[0, 0] += 1
    elif dis:
        o[0, 1] += 1
    elif mut:
        o[1, 0] += 1
    else:
        o[1, 1] += 1
    return o


def get_o(edge_string, gene):
    o=np.array([[0,0],[0,0]])
    obs_results = np.array([obs(p, edge_string, gene) for p in patients])
    return np.sum(obs_results, axis=0)

    # for result in map(obs, patients, itertools.repeat(edge_string, len(patients)), itertools.repeat(gene, len(patients))):
    #     o=o+result
    # return o


def chi(edge, gene, edge_string):
    '''
    Calculate chi2 statistic for a given gene and edge

    Parameters
    ----------
    edge : list --> first value is source, second value is target
        
    gene : string --> name of the mutated gene to analyse wrt the given edge
        
    
    edge_string : string --> edge in the string format --> (gene1,gene2)
        

    Returns
    -------
    a : float
        Chi-square value for that gene and edge

    '''
    o = list(get_o(edge_string, gene))
    row_totals = np.array([np.sum(o, axis=1)])
    col_totals = np.array([np.sum(o, axis=0)])
    n = np.sum(o)
    e = np.dot(row_totals.T, col_totals) / n
    chisq, p = chisquare(o, e)
    chisq = np.sum(chisq)
    p = 1 - chi2.cdf(chisq, 1)   
    return p



# only if not all edges are goint to be analysed:
# print("Opening edge counts...")
# with open('filtered_data/'+cancer+'/edges_counts.json') as h:
#     edges_counts = json.load(h)




print("Checking if gene sets data already exists...")

if exists('filtered_data/'+cancer+'/gene_sets_kappa'+str(kappa)+'.csv'): 
    gene_sets = pd.read_csv('filtered_data/'+cancer+'/gene_sets_kappa'+str(kappa)+'.csv')
else:
	gene_sets = get_relations(kappa,cancer, patients, dic_patients_mutations)
    
    
    
'''
iterates through every patient and computes patient-specific results
'''
for p in patients:
    print(p + ' calculating...')
    
    data = pd.DataFrame()  # Results dataframe that will be stored in the results folder
    genes_sigs = []  # Drivers (chi2 p_val < 0.05)
    chis_sig = []  # P_vals that were < 0.05 (p_vals of drivers)
    edge_strings=[] # dysregulations associated to drivers
    analysed_genes = set()  # To analyze every gene only once
    gene_sets_p = gene_sets[gene_sets['patients'] == p] # acces for the given patient the mutations and genes connected to mutations
    genes = set(gene_sets_p['connected dysregulations']) # get all mutations and genes connected to mutations, the rest of genes are not of interest for analysis

    for edge in dic_patients_edges[p]:
        edge_string = edge
        edge = [e.strip("() '") for e in edge.split(', ')]
        edge0, edge1 = edge
        
        if edge0 in genes and edge0 not in analysed_genes: # make sure the edge is not analysed and the tf is mutated or connected to mutation, otherwise not of interest
            analysed_genes.add(edge0)
            target_genes_set = set(gene_sets_p[gene_sets_p['connected dysregulations'] == edge0]['genes']) # get candidate drivers (mutations that are connected to the given dysregulation)
            visited = set()
            for gene in target_genes_set: # iterate every candidate driver
                key = str([edge, gene])
                if key not in visited: # make sure the combination edge and gene is not already analysed
                    visited.add(key)
                    chi_g_e = chi(edge, gene, edge_string)
                    if chi_g_e < 0.05: # if p_val < 0.05 then the mutation is considered a driver
                        edge_strings.append(edge_string)
                        chis_sig.append(chi_g_e)
                        genes_sigs.append(gene)
                   
              
            
   
    data = pd.DataFrame({'drivers': genes_sigs, 'p-values': chis_sig, 'dysregulations': edge_strings, 'patient': p})
    print('Number of drivers:', len(genes_sigs))
    if genes_sigs:
        if scores:
            df = calculate_page_rank(data, d=0.85)  # Concatenate all data at once
            df.to_csv(f'results/{cancer}/{p}_scores.csv', index=False)
        else:
            data.to_csv(f'results/{cancer}/{p}.csv', index=False)

    print(p + ' results computed')


print("FINISHED")
