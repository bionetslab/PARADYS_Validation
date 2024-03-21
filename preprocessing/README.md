# Preprocessing workflow

In order to preprocess input data for the application of PARADYS analysis, you have to run 
the script 'setup.py' with paths and filenames according to your input data storage. Below you 
can find the necessary changes you have to make in order to apply 'setup.py' to the desired dataset for
replication.

## Mutation input

Depending where the SNV and CNV dataframe files are stored, you have to adapt the corresponding 
paths in `snv_path` and `cnv_path`. We assume files to be tab-separated. Furthermore, the variables
`effect`, `gene`, and `sample` need to be set to the corresponding names of the columns in the
dataframes where variant effect, gene name, and sample ID are stored, respectively.

## Dysregulation network input

Depending on where your input dysregulation files (we assume them to be given in the .fea format)
are located, you need to update the variable `network_path` accordingly. Directly above in the 
variable `is_network_binarized` you need to specify whether the given input networks are already
binarized, i.e. entries in the dysregulation dataframe only take values 0 (edge not present) or 
1 (edge present).

## Output files

All output files are written to the `output_path`. This includes the following data for the 
mutational input:

- dic_patients_mutations.json: dictionary, keys are patients and values are a list of mutated genes for the given patient

- dic_mutations_patients.json: dictionary, keys are genes and patients are a list of patients that have the given gene mutated
- genes_M.json: list of genes that were mutated at least for one patient
- patients_M.json: list of patients that were present in the original mutation data



The output files of the dysregulation networks preprocessing include:

- genes_networks.json: all genes that were present in the networks (at least once)
- patients.json: all patients that have network and mutation data (intersection with patients_M)
- dic_patients_edges.json: dictionary, keys are patients and values are edges present in the network (as '(gene1,gene2)')
- dic_patients_graphs.json: dictionary, keys are patients and values are dictionaries that have tfs as keys and genes connected to them as values
- edges_counts.json: dictionary, keys are edges and values are the number of patients whose networks contain the given edge (to be used only for analysing frequent edges instead of all)
