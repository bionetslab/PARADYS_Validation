# PARADYS

In this repository, we present the code used to generate validation results presented in our paper "Patient specific ranking of genes driving dysregulations in cancer" for reproduction purposes.

# Usage

## Input preprocessing

In the `data/` directory, you can find our used dysregulation network files, as well as mutation and phenotype files. For usage, please first decompress the desired network and mutation file.

For preprocessing the desired input data, you need to navigate into `preprocessing/` and run the python script `setup.py`. You need to consider making the following changes in the script to run on your desired input:

### Mutation input

Depending where the SNV and CNV dataframe files are stored, you have to adapt the corresponding paths in `snv_path` and `cnv_path`. We assume files to be tab-separated. Furthermore, the variables `effect`, `gene`, and `sample` need to be set to the corresponding names of the columns in the dataframes where variant effect, gene name, and sample ID are stored, respectively.

### Dsyregulation network input

Depending on where your input dysregulation files (we assume them to be given in the .fea format) are located, you need to update the variable `network_path` accordingly. Directly above in the variable `is_network_binarized` you need to specify whether the given input networks are already binarized, i.e. entries in the dysregulation dataframe only take values 0 (edge not present) or 1 (edge present).

For a detailed description of the written output files, see `preprocessing/README.md`.

## Analysis

After preprocessing, you can detect patient-specific drivers for a given cohort by using the script `analysis/analysis.py`. You need to adapt the following parts in the script to fit your desired input data:

- kappa: Specifies the size of the considered shortest-path neighborhood in the dysregulation networks.
- input_dir: Name of directory that contains preprocessed input data. The list of necessary input data can be found in the `preprocessing/` directory.
- output_dir: Name of the output directory where to store patient-specific result files.

The resulting per-patient driver files are stored in `output_dir`. In such a file, the column `'drivers'` stores the detected driver gene, `'p-values'` the associated P-value for the edge contained in the column `'dysregulations'`, and `'patient'` the ID of the corresponding patient.

## Validation

In the directory `validation/` you can find scripts for all analyses presented in our paper. This includes clustering, PageRank driver scoring, plotting, and FDR calculation. 

## Results

The `results/` directory contains the data on which our validation results and plots are based upon.
