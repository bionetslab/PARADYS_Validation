import pandas as pd
import numpy as np
import graph_tool.all as gt
import os
import json

def compute_edge_patient_dict(networks : dict) -> dict:
    """Computes inverse dictionary of networks dictionary.

    Args:
        networks (dict): Keys are patients, values are sets of edges.

    Returns:
        dict: Keys are edges, values are sets of patients.
    """
    edge_patient_dict = dict()
    for pat in networks.keys():
        patient_edges = networks[pat]
        for edge in patient_edges:
            if not edge in edge_patient_dict.keys():
                edge_patient_dict[edge] = set()
            edge_patient_dict[edge].add(pat)
    
    return edge_patient_dict

def run_pagerank(patient_data : pd.DataFrame, edge_patient_dict : dict, 
                 number_patients : int) -> dict:
    """Run PageRank ranking algorithm of detected driver scores.

    Args:
        patient_data (pd.DataFrame): DataFrame storing patient-specific drivers and dysregulations.
        edge_patient_dict (dict): Keys are edges, values are sets of patients that possess this edge.
        number_patients (int): Number of patients in the whole cohort.

    Returns:
        dict: Keys are drivers, values are scores for each driver.
    """
    drivers = set(patient_data['drivers'])
    drivers_array = np.array(patient_data['drivers'])
    edges = set(patient_data['dysregulations'])
    edges = {eval(x) for x in edges}
    edges_array = np.array(patient_data['dysregulations'])
    pvalue_array = np.array(patient_data['score'])
    nodes = list(drivers.union(edges))
    # Init graph object.
    graph = gt.Graph(directed=True)
    
    # Add nodes to graph.
    vertices = graph.add_vertex(len(nodes))
    
    # Add labels to nodes.
    vertex_names = graph.new_vertex_property("string")
    vertex_weights = graph.new_vertex_property("double")
    edge_weights = graph.new_edge_property("double")
    
    # Set names and weights on vertices.
    for v, node in zip(vertices, nodes):
        vertex_names[v] = str(node)
        if node in drivers:
            vertex_weights[v]=0
        else: # "Node" is actually a dysregulation edge.
            if not node in edge_patient_dict.keys():
                vertex_weights[v] = 1/number_patients
            else:
                vertex_weights[v] = len(edge_patient_dict[node])/number_patients
    
    # Normalize vertex weights to unit-sum (otherwise centrality values diverge here).
    vertex_weights_scaled = vertex_weights.get_array()/np.sum(vertex_weights.get_array())
    vertex_weights.a = vertex_weights_scaled
    
    # Add edges between connected dysregulation edges with weight 1.
    for source in edges:
        for target in edges:
            if source[0]==target[1]:
                e = graph.add_edge(nodes.index(source), nodes.index(target))
                edge_weights[e]=1
    
    # Add edges between genes and edges with weight 1-log(pval+1).
    for driver in drivers:
        # Extract corresponding dysregulation edges to driver.
        driver_indices = np.where(drivers_array == driver)[0]
        driver_edges = edges_array[driver_indices]
        driver_pvalues = pvalue_array[driver_indices]
        for (edge, pvalue) in zip(driver_edges.tolist(), driver_pvalues.tolist()):
            e = graph.add_edge(nodes.index(eval(edge)), nodes.index(driver))
            edge_weights[e]=1-np.log10(pvalue+1)
    
    # Apply PageRank algorithm using vertex and edge weights.
    pagerank_prop = graph.new_vertex_property("float")
    pagerank = gt.pagerank(graph, damping=0.85, pers=vertex_weights,
                                      weight=edge_weights, prop=pagerank_prop)
    pagerank_res = pagerank.get_array()
    driver_names = {str(x) for x in drivers}
    sorted_drivers = [vertex_names[i] for i in np.argsort(pagerank_res) if vertex_names[i] in driver_names]
    sorted_drivers.reverse()
    sorted_scores = [pagerank_res[i] for i in np.argsort(pagerank_res) if vertex_names[i] in driver_names]
    sorted_scores.reverse()
    
    out_dict = {'driver': sorted_drivers , 'score': sorted_scores}
    return out_dict 

if __name__=="__main__":
    
    # Read drivers file.
    drivers_path = "../results/drivers/PRAD_DysRegNet_TPM.csv"
    out_dir = "pagerank_scores/"
    drivers_df = pd.read_csv(drivers_path, sep=';')
    
    # Load patient-edge data (from dysregulation network).
    patient_edge_file = '../data/patient_edges_PRAD.json'
    with open(patient_edge_file) as f:
        dic_patient_edges = json.load(f)
    # Transform edges strings into tuples.
    for key, values in dic_patient_edges.items():
        dic_patient_edges[key] = [eval(x) for x in values]
    
    # Compute edge-patient information for easier lookup how many patients have one specific edge.
    edge_patient_dict = compute_edge_patient_dict(dic_patient_edges)
    all_patients = set(dic_patient_edges.keys())
    number_patients = len(all_patients)
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # Extract all patients.
    patients = set(drivers_df['patient'])
    
    # Process each patient indivually with PageRank analysis.
    for pat in patients:
        
        # Extract corresponding drivers and dysregulations.
        patient_df = drivers_df[drivers_df['patient']==pat]
        print(patient_df)
        
        pagerank_dict = run_pagerank(patient_df, edge_patient_dict, number_patients)
        pagerank_df = pd.DataFrame(pagerank_dict)
        
        # Store pagerank results in individual file.
        file_name = f"{pat}_scores.csv"
        pagerank_df.to_csv(out_dir+file_name, sep='\t')
        
        