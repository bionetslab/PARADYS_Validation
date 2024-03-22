import json
import numpy as np
import pandas as pd
import os
from scipy.stats import chisquare, chi2
import networkx as nx

def calculate_page_rank(data, d=0.85, directed=True):
    '''
    Creates a linegraph and applies pagerank algorithm to get driver scores.

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
        Same as input, with additional column with scores and ordered in descending score.

    '''
    # Extract unique drivers and edges.
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

    # Preprocess edges to remove unwanted characters and create a dictionary for quick access.
    processed_edges = {edge.replace("(", "").replace(")", "").replace("'", "").replace(" ", ""): edge for edge in target_edges}
    
    # Add edges between connected edges with weight 1.
    for edge in target_edges:
        edge1 = edge.partition(', ')[2].replace("(", "").replace(")", "").replace("'", "").replace(" ", "")
        for edge2 in target_edges:
            edge0 = edge2.partition(', ')[0].replace("(", "").replace("'", "").replace(" ", "")
            if edge0 == edge1:
                G.add_edge(edge, processed_edges[edge2], weight=1)
    
    # Add edges between genes and edges with weight -log(p-value).
    for driver in total_drivers:
        data_driver = data[data['drivers'] == driver]
        for edge in data_driver['dysregulations']:
            G.add_edge(edge, driver, weight=-np.log(data_driver[data_driver['dysregulations'] == edge]['p-values'].values[0]))

    # Apply PageRank algorithm and store scores.
    pr = nx.pagerank(G, alpha=d, personalization=vertex_weights, weight='weight')
    # Add scores as a new column in the original data DataFrame.
    data['PageRank_Score'] = data['drivers'].map(pr)

    # Sort the DataFrame based on the scores
    data = data.sort_values(by='PageRank_Score', ascending=False)

    return data


def get_relations(kappa, cancer, patients, dic_patients_mutations, 
                  dic_patients_graphs, output_dir, save_file=False):
    '''
    Computes and stores all connections between mutations and dysregulations.

    Parameters
    ----------
    kappa : int
        Maximum shortest path length between edge and mutated gene.
    cancer : str
        Cancer type.
    patients : list
        List of patient IDs.
    dic_patients_mutations : dict
        Dictionary containing mutated genes for each patient.
    dic_patients_graphs: dict
        Dictionary with keys as patients and values are dictionaries that have TFs as keys and
        genes connected to them as values.
    output_dir : str
        Name of output directory for storing kappa neighborhood relation file.
    save_file : bool
        Whether or not to save the computed resulting neighborhood relation csv file (large!).

    Returns
    -------
    data : DataFrame
        For every mutated gene for every patient, list of dysregulations connected to it 
        within a max shortest-path distance kappa.
        Columns are genes ('genes'), patients ('patients') and 
        genes connected to the mutation ('connected dysregulations').

    '''
    found_mutations=[] 
    found_patients=[] 
    connected_genes=[] 
    
    for p in patients:
        genes_p = dic_patients_mutations.get(p, [])
        edges = dic_patients_graphs.get(p, {})
        # Process all mutated genes.
        for g in genes_p:
            if g in edges:
                connected_dysregulations = set()
                genes_stack = [g]
                # Compute kappa neighborhood.
                for _ in range(kappa):
                    next_genes = set()

                    for gene in genes_stack:
                        if gene in edges:
                            next_genes.update(set(edges[gene]) - connected_dysregulations)

                    connected_dysregulations.update(next_genes)
                    genes_stack = list(next_genes)

                if connected_dysregulations:
                    for dysregulation in connected_dysregulations:
                        found_mutations.append(g)
                        found_patients.append(p)
                        connected_genes.append(dysregulation)
    
    data = pd.DataFrame({
        'connected dysregulations': connected_genes,
        'genes': found_mutations,
        'patients': found_patients
    })
    
    # Store neighborhood file if desired.
    if save_file:
        data.to_csv(output_dir+'gene_sets_kappa'+str(kappa)+'.csv', index=False)
    
    return data


def mutation_condition (gene,patient):
    '''
    Check if gene is mutated for given patient.

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
    Check if dysregulation edge is present in patient's network.

    Parameters
    ----------
    edge : string
    patient : string

    Returns
    -------
    bool
        True if edge is present in the patient's dysregulatory network, False otherwise.

    '''
    if edge in dic_patients_edges[patient]:
        return True
    else:
        return False


def obs(p, edge_string, gene):
    '''
    Create per-patient observation matrix for Chi2 test.

    Parameters
    ----------
    p : string
        Patient ID.
    edge_string : string
        Edge in string format, e.g. (gene1, gene2).
    gene : string
        Gene label.
        
    Returns
    -------
    o : DataFrame
        2x2 observation matrix.

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


def get_o(edge_string, gene, patients):
    """Computes aggregated observation matrix over all patients.

    Args:
        edge_string (str): Edge in the string format, e.g. (gene1, gene2).
        gene (str): Name of the mutated gene to analyse wrt the given edge.
        patients (list): Patients in considered cohort.
    Returns:
        np.array: 2x2 observation matrix.
    """
    o=np.array([[0,0],[0,0]])
    obs_results = np.array([obs(p, edge_string, gene) for p in patients])
    return np.sum(obs_results, axis=0)


def chi(gene, edge_string, patients):
    '''
    Calculate chi2 test statistic for given gene and edge.

    Parameters
    ----------
    gene : string 
        Name of the mutated gene to analyse wrt the given edge.
    edge_string : string 
        Edge in the string format, e.g. (gene1,gene2).
    patients: list
        List of all patients in the cohort.

    Returns
    -------
    a : float
        Chi-square value for that gene and edge

    '''
    # Compute observation matrix.
    o = list(get_o(edge_string, gene, patients))
    row_totals = np.array([np.sum(o, axis=1)])
    col_totals = np.array([np.sum(o, axis=0)])
    n = np.sum(o)
    # Compute expectance matrix.
    e = np.dot(row_totals.T, col_totals) / n
    chisq, p = chisquare(o, e)
    chisq = np.sum(chisq)
    p = 1 - chi2.cdf(chisq, 1)   
    return p


if __name__ == "__main__":
    
    os.chdir("../")
    
    ### Input parameters.
    
    # Name of cancer cohort to be analyzed.
    cancer = 'PRAD_Decoy_075' 
    # Size of neighbor to consider for putative driver identification.
    kappa = 3
    # Whether to compute driver scores based on PageRank algorithm.
    scores = False
    # Input directory containing preprocessed files.
    input_dir = f'filtered_data/{cancer}/'
    # Output directory for storing result files.
    output_dir = f'results/{cancer}/'

    ### Loading of preprocessed data.
    
    print("Loading preprocessed data...")

    with open(input_dir+'patients.json') as i:
        patients = json.load(i)

    with open(input_dir+'genes_M.json') as i:
        genes_M = json.load(i)

    with open(input_dir+'dic_patients_edges.json') as i:
        dic_patients_edges = json.load(i)

    with open(input_dir+'dic_mutations_patients.json') as i:
        dic_mutations_patients = json.load(i)

    with open(input_dir+'dic_patients_mutations.json') as i:
        dic_patients_mutations = json.load(i)   
        
    with open(input_dir+'dic_patients_graphs.json') as i:
        dic_patients_graphs = json.load(i)

    # Checking if gene dysregulation neighborhood information has already been computed and stored.
    # Speeds up computation if this information can simply be loaded.
    gene_sets_file = 'gene_sets_kappa'+str(kappa)+'.csv'
    if os.exists(input_dir+gene_sets_file): 
        gene_sets = pd.read_csv(input_dir+gene_sets_file)
    else:
        save_file = True
        gene_sets = get_relations(kappa,cancer, patients, dic_patients_mutations, 
                                  dic_patients_graphs, input_dir, save_file)

    
    # Process all given patients.
    for p in patients:

        # Results dataframe that will be stored in the results folder.
        data = pd.DataFrame()
        # List of significant drivers.  
        genes_sigs = []
        # Pvalues of significant drivers.  
        chis_sig = []  
        # Dysregulations associated to drivers.
        edge_strings=[] 
        # To analyze every gene only once.
        analysed_genes = set()  
        
        # For the given patient, access mutations and genes connected to mutations.
        gene_sets_p = gene_sets[gene_sets['patients'] == p]
        # Sample all mutations and genes connected to mutations, the remaining genes are not of interest for analysis.
        genes = set(gene_sets_p['connected dysregulations'])

        # Process all patient-specific edges.
        for edge in dic_patients_edges[p]:
            edge_string = edge
            edge = [e.strip("() '") for e in edge.split(', ')]
            edge0, edge1 = edge

            # Only consider edges that are not analysed yet and whose tf is mutated or connected to mutation.
            if edge0 in genes and edge0 not in analysed_genes: 
                analysed_genes.add(edge0)
                # Extract candidate drivers (i.e. mutations that are connected to the given dysregulation).
                target_genes_set = set(gene_sets_p[gene_sets_p['connected dysregulations'] == edge0]['genes']) 
                visited = set()
                # Process every candidate driver.
                for gene in target_genes_set: 
                    key = str([edge, gene])
                    # Make sure the combination of edge and gene is not already analysed.
                    if key not in visited: 
                        visited.add(key)
                        chi_g_e = chi(gene, edge_string, patients)
                        # Significance thresholding.
                        if chi_g_e < 0.05: 
                            edge_strings.append(edge_string)
                            chis_sig.append(chi_g_e)
                            genes_sigs.append(gene)

        # Create driver results dataframe.
        data = pd.DataFrame({'drivers': genes_sigs, 'p-values': chis_sig, 'dysregulations': edge_strings, 'patient': p})
        
        # If at least one driver has been detected, check if scoring is desired.
        if genes_sigs:
            if scores:
                df = calculate_page_rank(data, d=0.85) 
                df.to_csv(output_dir+f'{p}_scores.csv', index=False)
            else:
                data.to_csv(output_dir+f'{p}.csv', index=False)

    print("FINISHED")
