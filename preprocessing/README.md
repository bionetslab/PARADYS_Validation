# Preprocessing output files

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
