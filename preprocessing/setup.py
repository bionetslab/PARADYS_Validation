import os
import json
from collections import defaultdict
import pandas as pd
import mygene
import sklearn

def create_dictionary_genes_ids(cnvs : pd.DataFrame, output_path : str) -> dict:
    '''
    Transltaion of gene symbols to gene ids. Stores dictionary in the output path. 

    Parameters
    ----------
    cnvs : pd.DataFrame
        CNV dataframe with column named 'Gene Symbols'.
    output_path : str
        Name of the output path where to write files to.
    
    Returns
    -------
    dic_genes_ids : dict
        Keys are gene ids and values are gene symbols.

    '''

    mg = mygene.MyGeneInfo()
    geneList = cnvs['Gene Symbol'].to_list()
    dic_genes_ids=dict()

    # Reuse previously computed translation dict if already avaiable.
    # with open(f'{output_path}/dic_ensembl_genes.json') as i:
    #     dic_genes_ids=json.load(i)
    
    # Process all genes that appear in CNV DataFrame.
    for gene in geneList:
        d = gene.partition('.')[0]
        b = mg.querymany(d, scopes='ensembl.gene', fields='symbol', species='human', returnall=True)
        c=b['out']
        if 'symbol' in c[0].keys():
        #print(c[0]['symbol'])
            dic_genes_ids[gene]=c[0]['symbol']
    
    out_file = open(f'{output_path}/dic_ensembl_genes.json', "w")
    json.dump(dic_genes_ids, out_file)
    out_file.close()
    
    return dic_genes_ids


def binarize(networks : pd.DataFrame) -> pd.DataFrame:
    '''
    Creates a binarized edge matrix given threshold for edge weights.
    To be used only in case that original dysregulation networks are weighted.

    Parameters
    ----------
        networks (DataFrame): Original (weighted) networks.

    Returns
    -------
        binarized_data (DataFrame): DataFrame with binarized (thresholded) edge weights.
    '''

    # Remove first column in networks dataframe.
    d=networks.iloc[: , 1:].to_numpy()
    
    # Threshold edge weights.
    d=sklearn.preprocessing.binarize(d,threshold=0.1)
    #d=networks.iloc[: , :-1].to_numpy()
    
    edges = networks.columns.to_list()
    samples = networks['patient id'].to_list()

    # First edge (column) are patients, drop this.
    edges.pop(0)
    
    # Reconstruct binarized DataFrame based on thresholded edge weight matrix.
    binarized_data=pd.DataFrame(d, index=samples, columns=edges)
    
    return binarized_data


def set_up_mutations(snv_data, cnv_data, effect, gene, sample, cancer, output_path):
    '''
    Aggregates CNVs and SNVs and stores patient-mutation dictionaries as input for analysis.

    Parameters
    ----------
    snv_data : pd.DataFrame
        DataFrame storing SNV mutation information.
    cnv_data : pd.DataFrame
        DataFrame storing CNV mutation information.
    effect : String
        Name of the column of the original mutation file that specifies the effect of every mutation
    gene : String
        Name of the column of the original mutation file that specifies the gene ID
    sample : String
        Name of the column of the original mutation file that specifies the sample ID
    cancer : String
        Name of the cohort
    output_path: str
        Name of output path where files shall be stored.

    Returns
    -------
    None.
    It writes the following files:
        dic_patients_mutations --> dictionary, keys are patients and values are a list of mutated genes for the given patient
        dic_mutations_patients --> dictionary, keys are genes and patients are a list of patients that have the given gene mutated
        genes_M --> list of genes that were mutated at least for one patient
        patients_M --> list of patients that were present in the original mutation data
    '''
    
    # Effect entry of variant that is supposed to be skipped, i.e. not considered as variant.
    type_to_filter_out = 'synonymous_variant'
    
    genes = []
    patients = []

    # Process SNV data.
    for index, _ in snv_data.iterrows():
        # Check if variant is supposed to be ignored.   
        if snv_data.loc[index,effect]!=type_to_filter_out:
            # If not, append correspondig gene to mutated genes list.
            genes.append(snv_data.loc[index,gene])
            patients.append(snv_data.loc[index,sample][:-1])
    
    # Process CNV data.
    if cancer == 'BRCA':
        # BRAC cohort needs additional translation of gene symbols.
        dic_genes_ids = create_dictionary_genes_ids(cnv_data, output_path)
        patients_cnvs = list(cnv_data.columns)[1:-1]
        for index, _ in cnv_data.iterrows():
            for p in patients_cnvs:
                # Check if gene contains copy number variation, i.e. entry unequal zero.
                if cnv_data.loc[index,p]!=0:
                    if cnv_data.loc[index,'Gene Symbol'] in dic_genes_ids.keys():
                        # Translate gene ID based on myGene dictionary.
                        genes.append(dic_genes_ids[cnv_data.loc[index,'Gene Symbol']])
                        patients.append(p[:-1])

    # Process CNVs for PRAD/COAD cohort. 
    else:
        patients_cnvs = list(cnv_data.columns)[1:]
        for index, _ in cnv_data.iterrows():
            for p in patients_cnvs:
                if cnv_data.loc[index,p]!=0:
                    genes.append(cnv_data.loc[index,'Gene Symbol'])
                    patients.append(p)
                    
    # Create output dictionaries storing patient-mutation information.            
    dic_patients_mutations = defaultdict(list)
    dic_mutations_patients = defaultdict(list) 
    for i,j in zip(genes,patients):
        dic_mutations_patients[i].append(j)
        dic_patients_mutations[j].append(i)
    
    # Remove duplicated from patients and genes list.
    patients=list(set(patients))
    genes=list(set(genes))
    
    # Write output files.
    pat_muts_file = open(f'{output_path}/dic_patients_mutations.json', "w")
    json.dump(dic_patients_mutations, pat_muts_file)
    pat_muts_file.close()
    
    genes_file = open(f'{output_path}/genes_M.json', "w")
    json.dump(genes, genes_file)
    genes_file.close()
    
    mut_pats_file = open(f'{output_path}/dic_mutations_patients.json', "w")
    json.dump(dic_mutations_patients, mut_pats_file)
    mut_pats_file.close()

    patients_file = open(f'{output_path}/patients_M.json', "w")
    json.dump(patients, patients_file)
    patients_file.close() 


def set_up_networks(binarized_networks : pd.DataFrame, patients_are_index : bool):
    '''
    Creates preprocessed files of patient-dysregulation data to be used in the analysis.

    Parameters
    ----------
    binarized_networks : DataFrame
        Indices are patients and columns are edges (pairs of strings of genes, separated by coma).
        Cell entries are 0 or 1 if the given edge was present in the patient network.
    patients_are_index : bool
        True if patient names are index of the input DataFrame. False means that they are stored in
        the column 'patient id'.

    Returns
    -------
    None.
    It writes the following files:
        genes_networks --> all genes that were present in the networks (at least once)
        patients --> all patients that have network and mutation data (intersection with patients_M)
        dic_patients_edges --> dictionary, keys are patients and values are edges present in the network (as '(gene1,gene2)')
        dic_patients_graphs --> dictionary, keys are patients and values are dictionaries that have tfs as keys and genes connected to them as values
        edges_counts --> dictionary, keys are edges and values are the number of patients whose networks contain the given edge (to be used only for analysing frequent edges instead of all)
    '''

    # Make sure that we have mutation info from every patient of the network original data.
    # If not remove the patient.
    with open(f'{output_path}/patients_M.json') as j:
        patients = set(json.load(j))

    # Set patients as index if necessary.
    if not patients_are_index:
        binarized_networks.set_index('patient id', inplace=True)
        binarized_networks.drop(columns=['patient id'], inplace=True)
    
    patients_nets = set(binarized_networks.index.to_list())
    # Intersect patients from dysregulation network and mutation data.
    patients_common = patients.intersection(patients_nets)

    columns = binarized_networks.columns 
    genes_networks = []
    edges = []
    all_dic_edges = []
    
    # Iterate over all common patients.
    for p in patients_common:

        # Extract and collect per-patient dysregulation networks.
        edges_p = [col for col in columns if binarized_networks.loc[p, col] != 0]
        edges.append(edges_p)
        
        sources, targets = [], []
        # Extract source and target genes of all per-patient edges.
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
            
    # Write output files.
    genes_networks = list(set(genes_networks))
    gene_nets_file = open(f'{output_path}/genes_networks.json', "w")
    json.dump(genes_networks, gene_nets_file)
    gene_nets_file.close()

    patients_file = open(f'{output_path}/patients.json', "w")
    json.dump(list(patients_common), patients_file)
    patients_file.close()
    
    zip_iterator = zip(patients_common, edges)
    dic_patients_edges = dict(zip_iterator)
    pat_edges_file = open(f'{output_path}/dic_patients_edges.json', "w")
    json.dump(dic_patients_edges, pat_edges_file)
    pat_edges_file.close()

    zip_iterator = zip(patients_common, all_dic_edges)
    dic_patients_graphs = dict(zip_iterator)
    pat_graphs_file = open(f'{output_path}/dic_patients_graphs.json', "w")
    json.dump(dic_patients_graphs, pat_graphs_file)
    pat_graphs_file.close()
    
    counts=[]
    edges_totals = binarized_networks.columns
    # Compute number of times a dysregulation edge is contained in per-patient networks.
    for col in binarized_networks.columns:
        counts.append(binarized_networks[col].sum())
    zip_iterator = zip(edges_totals, counts)
    edges_counts = dict(zip_iterator)
    edges_file = open(f'{output_path}/edges_counts.json', "w")
    json.dump(edges_counts, edges_file)
    edges_file.close()


if __name__ == '__main__':

    os.chdir("../")

    # Cancer cohort to be analyzed. # In our case, selection from {PRAD, COAD, BRCA}.
    cancer='PRAD'
    
    ## Output path specification.
    output_path = f'filtered_data/{cancer}'
    
    # Make sure that output directory exists.
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    ##  Data format specification concerning mutation preprocessing.

    # Name of the column of the mutation data file that specifies the effect of a given mutation.
    effect = 'effect'
    # Name of the column of the mutation data file that specifies the name of a given gene.
    gene = 'gene'
    # Name of the column of the mutation data file that specifies the sample label.
    sample = 'Sample_ID'

    # Path to mutation data.
    snv_path = f'data/{cancer}/TCGA-PRAD.muse_snv.tsv'
    cnv_path = f'data/{cancer}/TCGA-PRAD_CNVS.tsv'

    # Read in mutation data.
    snvs=pd.read_csv(snv_path, sep='\t')
    cnvs=pd.read_csv(cnv_path, sep='\t')

    ##  Dysregulation network input specification.
    
    # Set this to False if the input dysregulation networks are weighted (i.e. not thresholded yet).
    is_network_binarized = True 
    network_path = f'data/{cancer}/SSN_nets.fea'

    # Load dysregulation network.
    print("Loading network feather file...")
    networks = pd.read_feather(network_path)

    # Threshold dsyregulation networks if necessary.
    if not is_network_binarized:
        networks = binarize(networks)

    # Preprocess mutations.
    print("Preprocessing mutations...")
    set_up_mutations(snvs, cnvs, effect, gene, sample, cancer, output_path)
    
    # Preprocess dysregulation networks.
    print("Preprocessing networks...")
    set_up_networks(networks)

