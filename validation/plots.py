import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd
import numpy as np
import seaborn as sns
import scipy as sc
import os
    
def get_patients_of_cluster(cluster_file : str, cluster_identifier : float, cluster_level : int) -> list:
    """Extract sample IDs of given cluster and returns them as list.
    """
    # Open cluster information file and extract corresponding samples.
    cluster_info = pd.read_csv(cluster_file)
    column_identifier = f'{cluster_level}_cluster_label'
    cluster_info.set_index('Patient', inplace=True)
    cluster_samples = cluster_info.index[cluster_info[column_identifier]==cluster_identifier].to_list()
    return cluster_samples
    
def get_phenotypes_of_cluster(cluster_file : str, cluster_identifier : float, cluster_level : int,
                              phenotype_file: str, phenotype : str) -> list:
    """Extract list of phenotype labels of all samples in given cluster (e.g. gleason/PSA).
    """
    samples = get_patients_of_cluster(cluster_file, cluster_identifier, cluster_level)
    # Extract corresponding column from phenotype data file.
    data = pd.read_csv(phenotype_file)
    # Extract cluster samples from phenotype file.
    data_cluster = data[data["_PATIENT"].isin(samples)]
    return data_cluster[phenotype].to_list()


def get_survival_times_of_cluster(cluster_file : str, cluster_identifier : float, cluster_level : int,
                                  survival_file : str) -> list:
    """Extract list of survival times of patients in given cluster.

    Args:
        cluster_file (str): Path of file where cluster-sample relation is stored.
        cluster_identifier (float): Float describing cluster ID.
        cluster_level (int): Level of cluster in hierarchy (i.e. 4 for cluster 4.1)
        survival_file (str): Path of file storing survival time.

    Returns:
        list: Survival times of patients in cluster.
    """
    samples = get_patients_of_cluster(cluster_file, cluster_identifier, cluster_level)
    # Extract corresponding column from phenotype data file.
    data = pd.read_csv(survival_file, sep='\t')
    # Extract cluster samples from phenotype file.
    data_cluster = data[data["_PATIENT"].isin(samples)]
    return data_cluster['OS.time'].to_list()

def compute_mwu_pval(arr1 : np.ndarray, arr2 : np.ndarray, mode : str) -> float:
    """Compute MWU pvalue of two given feature vectors.

        Args:
            mode (str): Alternative hypothesis of MWU test. Either 'two-sided', 'greater' or 'less'.
    """
    _, p = sc.stats.mannwhitneyu(arr1, arr2, alternative=mode)
    return p

def plot_gleason_os_xcell(cluster_file_path : str, phenotype_file : str, survival_file : str,
                          xcell_results_path : str):
    """Plots Figure containing gleason scores and OS time per 4.x cluster as well as xcell analysis 
    results.
    """
    # Create one large mosaic plot.
    fig, axs = plt.subplot_mosaic([['A', 'A', 'B', 'B'], ['C', 'C', 'D', 'D']])
    fig.set_size_inches(10, 10)
    axs_list = list(axs.values())
    
    # Add gleason and OS subplot to mosaic plot.
    cluster_level = 4
    label_fontsize = 15
    
    clusters = [4.0, 4.1, 4.2, 4.3]
    phenotype = "gleason_score"
    gleason_data = []
    os_data = []
    for clust in clusters:
        gleason = get_phenotypes_of_cluster(cluster_file_path, clust, cluster_level, 
                                        phenotype_file, phenotype)
        gleason_data.append(np.array(gleason))
        os = get_survival_times_of_cluster(cluster_file_path, clust, cluster_level, survival_file)
        os_data.append(np.array(os))    
    
    axs_list[0].boxplot(gleason_data, patch_artist=True, medianprops = dict(color="firebrick",linewidth=2.3), boxprops=dict(facecolor='mistyrose'))
    axs_list[0].set_ylabel('Gleason score', fontsize=label_fontsize)
    axs_list[0].set_xticks([y + 1 for y in range(len(gleason_data))], labels=['4.0', '4.1', '4.2', '4.3'])
    axs_list[0].set_xlabel('Cluster', fontsize=label_fontsize)
    axs_list[0].grid(True, axis='y')
    axs_list[0].yaxis.set_tick_params(labelsize=label_fontsize)
    axs_list[0].xaxis.set_tick_params(labelsize=label_fontsize)
    
    # Add OS time subplot.
    axs_list[1].boxplot(os_data, patch_artist=True, medianprops=dict(color='firebrick', linewidth=2.3), boxprops=dict(facecolor='mistyrose'))
    axs_list[1].set_ylabel('Survival time [days]', fontsize=label_fontsize)
    axs_list[1].set_xticks([y + 1 for y in range(len(os_data))], labels=['4.0', '4.1', '4.2', '4.3'])
    axs_list[1].set_xlabel('Cluster', fontsize=label_fontsize)
    axs_list[1].grid(True, axis='y')
    axs_list[1].yaxis.set_tick_params(labelsize=label_fontsize)
    axs_list[1].xaxis.set_tick_params(labelsize=label_fontsize)
    
    # Add mean xcell scores to third mosaic plot.
    clusters = ['4.0', '4.1', '4.2', '4.3']
    immune_related_cells = [
    "aDC",
    "DC",
    "iDC",
    "pDC"]
    index = sorted(immune_related_cells)
    columns = clusters
    means_df = pd.DataFrame(index=index, columns=columns)
    
    for clust in clusters:
        mean_file1 = xcell_results_path+f'cluster{clust}/xcell_score_means.csv'
        mean_df1 = pd.read_csv(mean_file1, sep='\t')
        # Subset only immune cells.
        mean_df1 = mean_df1[mean_df1['Unnamed: 0'].isin(immune_related_cells)]
        mean_df1.sort_values(by=['Unnamed: 0'], inplace=True)
        means_df[clust] = mean_df1['mean'].to_list()
    
    sns.set(font_scale=1.3)
    s = sns.heatmap(means_df, annot=True, fmt='.1f', ax=axs_list[2], cbar_kws={'label': 'Mean xCell score'})
    s.set_xlabel('Cluster', fontsize=label_fontsize)
    s.set_ylabel('Cell type', fontsize=label_fontsize)
    s.set_xticklabels(s.get_xmajorticklabels(), fontsize=label_fontsize)
    s.set_yticklabels(s.get_ymajorticklabels(), fontsize=label_fontsize)
    
    # Add xcell pvalues to fourth mosaic plot.
    clusters = ['4.1', '4.2', '4.3']
    pvals = pd.DataFrame(index=immune_related_cells, columns=clusters, dtype=float)
    
    for clust in clusters:
        score_file1 = xcell_results_path+'cluster4.0/xcell_scores.tsv'
        score_file2 = xcell_results_path+f'cluster{clust}/xcell_scores.tsv'
        score_df1 = pd.read_csv(score_file1, sep='\t')
        score_df2 = pd.read_csv(score_file2, sep='\t')
        score_df1.set_index('Unnamed: 0', inplace=True)
        score_df2.set_index('Unnamed: 0', inplace=True)
        # Iterate over all immune cells and compute MWU pvalue for 4.0-4.x.
        for cell in immune_related_cells:
            scores1 = score_df1.loc[cell].to_list()
            scores2 = score_df2.loc[cell].to_list()
            pvalue = compute_mwu_pval(np.array(scores1), np.array(scores2), mode='greater')
            pvals.loc[cell, clust] = pvalue
    
    pvalues_signif = pvals.copy()
    # Log transform pvalues.
    pvalues_signif = pvalues_signif.transform(lambda x: -1*np.log10(x))
    print(pvalues_signif)
    print(pvalues_signif.index)
    print(pvalues_signif.columns)
    print(type(pvalues_signif.iloc[0,0]))
    sns.set(font_scale=1.5)
    s = sns.heatmap(pvalues_signif, annot=True, fmt='.2f', vmin=0, vmax=17, ax=axs_list[3], cbar_kws={'label': '-log10(p-value)'}) #annot_kws={"size": label_fontsize-1}
    s.set_xlabel('Target cluster', fontsize=label_fontsize) 
    #s.set_ylabel('Cell type')
    s.set_xticklabels(s.get_xmajorticklabels(), fontsize=label_fontsize)
    s.set_yticklabels(s.get_ymajorticklabels(), fontsize=label_fontsize)
    
    
    # Add mosaic labels above plots.
    for label, ax in axs.items():
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='x-large', va='bottom', fontfamily='sans-serif', fontweight='bold')
    
    fig.tight_layout()
    plt.savefig("gleason_os_xcell.pdf", format='pdf')
    
    
def plot_fdrs():
    """Plot FDR results from all used tools and input networks.
    """
    # Create one large mosaic plot.
    fig, axs = plt.subplot_mosaic([['A', 'B'], ['C', 'D']])
    fig.set_size_inches(10, 10)
    axs_list = list(axs.values())
    label_fontsize = 17
    
    # Axis label sizes.
    axis_label_sizes = "large"
    
    ### Add PARADYS DysRegNet FDRS in top left plot.
    file_path = '../results/fdrs/PARADYS_DysRegNet/fdrs_dysregnet_all_cohorts.csv'
    df = pd.read_csv(file_path, sep=';')
    df.rename(columns={"cohort" : "Cohort"}, inplace=True)
    df.sort_values(by=['Cohort'], inplace=True, ascending=False)
    sns.set_theme(style='darkgrid')
    sns.set(font_scale=3.0)
    sns.violinplot(x = df['decoy_rate'],
                y = df['FDR'], 
                hue = df['Cohort'],
                cut=0,
                ax=axs_list[0],
                linewidth=1.4)
    # axs.set_title("FDRs of PARADYS using DysRegNet networks")
    axs_list[0].legend(prop={'size':18})
    axs_list[0].set_ylabel("FDR", fontsize=label_fontsize)
    axs_list[0].set_xlabel("Rate of decoy mutations", fontsize=label_fontsize)
    axs_list[0].set_axisbelow(True)
    axs_list[0].grid(axis='y')
    axs_list[0].yaxis.set_tick_params(labelsize=label_fontsize)
    axs_list[0].xaxis.set_tick_params(labelsize=label_fontsize)
    
    ### Add PRODIGY FDRS on bottom right plot.
    file_path = '../results/fdrs/PRODIGY/fdrs_prodigy.csv'
    df = pd.read_csv(file_path, sep=',')
    rate_dict = {25 : '0.25', 5 : '0.50', 75 : '0.75'}
    df.replace({"Rate": rate_dict}, inplace=True)
    df.sort_values(by=['Cohort'], inplace=True, ascending=False)
    
    sns.set_theme(style='darkgrid')
    sns.set(font_scale=3.0)
    sns.violinplot(x = df['Rate'],
                y = df['FDR'], 
                hue = df['Cohort'],
                ax=axs_list[3],
                cut=0, linewidth=1.4)
    # axs.set_title("FDRs of PARADYS using DysRegNet networks")
    #axs_list[1].legend(loc='upper right', title='Cohort', title_fontsize='small')
    axs_list[3].get_legend().remove()
    axs_list[3].set_ylabel("FDR", fontsize=label_fontsize)
    axs_list[3].set_xlabel("Rate of decoy mutations", fontsize=label_fontsize)
    axs_list[3].set_axisbelow(True)
    axs_list[3].grid(axis='y')
    axs_list[3].yaxis.set_tick_params(labelsize=label_fontsize)
    axs_list[3].xaxis.set_tick_params(labelsize=label_fontsize)
    
    ### Add PARADYS with SSN networks plots on bottom left.
    cohorts = ['prad', 'coad']
    file_path = '../results/fdrs/PARADYS_SSN/'
    # Merge FDR files for different cohorts into one dataframe.
    combined_df = pd.DataFrame(columns=['patient', 'FDR', 'Rate', 'Cohort'])
    for cohort in cohorts:
        file = file_path+f'fdrs_{cohort}.csv'
        df = pd.read_csv(file)
        df['Cohort']=cohort.upper()
        df.drop(columns=['Unnamed: 0'], inplace=True)
        combined_df = pd.concat([combined_df, df], ignore_index=True)
        
    combined_df.sort_values(by=['Cohort'], inplace=True, ascending=False)
    
    sns.set_theme(style='darkgrid')
    sns.set(font_scale=3.0)
    sns.violinplot(x = combined_df['Rate'],
                y = combined_df['FDR'], 
                hue = combined_df['Cohort'],
                cut=0,
                ax=axs_list[2],
                linewidth=1.4)
    axs_list[2].get_legend().remove()
    axs_list[2].set_ylabel("FDR", fontsize=label_fontsize)
    axs_list[2].set_xlabel("Rate of decoy mutations", fontsize=label_fontsize)
    axs_list[2].set_axisbelow(True)
    axs_list[2].grid(axis='y')
    axs_list[2].yaxis.set_tick_params(labelsize=label_fontsize)
    axs_list[2].xaxis.set_tick_params(labelsize=label_fontsize)
    
    ### Add PARADYS with DysRegNet_NoDir networks FDRS on top right.
    cohorts = ['prad', 'coad', 'brca']
    file_path = '../results/fdrs/PARADYS_DysRegNet_NoDir/'
    # Merge FDR files for different cohorts into one dataframe.
    combined_df = pd.DataFrame(columns=['patient', 'FDR', 'Rate', 'Cohort'])
    for cohort in cohorts:
        file = file_path+f'fdrs_{cohort}.csv'
        df = pd.read_csv(file)
        df['Cohort']=cohort.upper()
        combined_df = pd.concat([combined_df, df], ignore_index=True)
        
    combined_df.sort_values(by=['Cohort'], inplace=True, ascending=False)
    
    sns.set_theme(style='darkgrid')
    sns.set(font_scale=3.0)
    sns.violinplot(x = combined_df['Rate'],
                y = combined_df['FDR'], 
                hue = combined_df['Cohort'],
                cut=0,
                ax=axs_list[1],
                linewidth=1.4)
    axs_list[1].get_legend().remove()
    axs_list[1].set_ylabel("FDR", fontsize=label_fontsize)
    axs_list[1].set_xlabel("Rate of decoy mutations", fontsize=label_fontsize)
    axs_list[1].set_axisbelow(True)
    axs_list[1].grid(axis='y')
    axs_list[1].yaxis.set_tick_params(labelsize=label_fontsize)
    axs_list[1].xaxis.set_tick_params(labelsize=label_fontsize)
    
    # Add mosaic labels above plots.
    for label, ax in axs.items():
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='x-small', va='bottom', fontfamily='sans-serif', fontweight='bold')
    
    fig.tight_layout()
    plt.savefig("fdrs.pdf", format='pdf')
    
    
def get_top_k_drivers(k : int, drivers_file : str, clusters : list, output_path : str):
    """Extract the top k most occuring drivers in a given list of clusters.

    Args:
        k (int): Number of top frequent drivers to return.
        drivers_file (str): Name of input file storing driver frequencies per cluster.
        clusters (list): List of cluster ID strings to consider (in the format '3.1').
        output_path (str): Path to output directory.
    """
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    drivers_df = pd.read_csv(drivers_file, sep=';')
    
    for clust in clusters:
        column_name = f'Driver_{clust[0]}_cluster_label_{clust}'
        cluster_df = drivers_df.sort_values(by=[column_name], ascending=False)
        top_drivers = list(cluster_df['driver'])[:k]
        file_name = f'top_{k}_cluster_{clust}.txt'
        with open(output_path+file_name, 'w') as f:
            for gene in top_drivers:
                f.write(f'{gene}\n')
    

def plot_digest_pvalues(digest_path : str, tools : list, clusters : list):
    """Plots pvalues of DIGEST analysis.

    Args:
        digest_path (str): Path to digest pvalue files.
        tools (list): List of tools to consider in digest_path.
        clusters (list): List of cluster IDs to consider.
    """
    # Iterate over all tools and combine cluster pvalues into one dataframe.
    dataframes = dict()
    for tool in tools:
        tool_path = digest_path+tool+"/"
        # Init tool's dataframe.
        tool_df = pd.DataFrame(index=['GO.BP', 'GO.CC', 'GO.MF', 'KEGG'], columns=clusters)
        for clust in clusters:
            cluster_file = tool_path+f'{tool}_pvalues_cluster{clust}.csv'
            cluster_df = pd.read_csv(cluster_file)
            cluster_df.rename(columns={'Unnamed: 0':'Database'}, inplace=True)
            for _, row in cluster_df.iterrows():
                pval_log = -1.0 * np.log10(row['JI-based'])
                tool_df.loc[row['Database'], clust] = np.float64(pval_log)
        # Append tool dataframe to collection of all dataframes.
        dataframes[tool]=tool_df
        
    # Init mosaic subplot.
    fig, axs = plt.subplot_mosaic([['A', 'A', 'B', 'B'], ['.', 'C', 'C', '.']])
    fig.set_size_inches(7, 7)
    axs_list = list(axs.values())
    
    # Plot heatmap for DysRegNet pvalues with -log transformed pvalues.
    pvals_dysregnet = dataframes['dysregnet']
    pvals_dysregnet = pvals_dysregnet[pvals_dysregnet.columns].astype(float)
    s = sns.heatmap(pvals_dysregnet, annot=True, fmt='.2f', vmin=1.3, vmax=3, ax=axs_list[0], cmap='magma', cbar_kws={'label': '-log10(p-value)'})
    s.set_xlabel('Patient Cluster (DysRegNet)')
    
    hitndrive_df = dataframes['hitndrive']
    hitndrive_df = hitndrive_df[hitndrive_df.columns].astype(float)
    s = sns.heatmap(hitndrive_df, annot=True, fmt='.2f', vmin=1.3, vmax=3, ax=axs_list[1], cmap='magma', cbar_kws={'label': '-log10(p-value)'})
    s.set_xlabel('Patient Cluster (HITNDRIVE)')
    
    prodigy_df = dataframes['prodigy']
    prodigy_df = prodigy_df[prodigy_df.columns].astype(float)
    s = sns.heatmap(prodigy_df, annot=True, fmt='.2f', vmin=1.3, vmax=3, ax=axs_list[2], cmap='magma', cbar_kws={'label': '-log10(p-value)'})
    s.set_xlabel('Patient Cluster (PRODIGY)')
    
    # Add mosaic labels above plots.
    for label, ax in axs.items():
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='x-large', va='bottom', fontfamily='sans-serif', fontweight='bold')
    
    fig.tight_layout()
    plt.savefig("digest_pvalues.pdf", format='pdf')
    
def change_cluster_labels(label : str) -> str:
    """Changes cluster labeling from 2.x_2.y to 2.x vs 2.y.

    Args:
        label (str): Current label to be changed.

    Returns:
        str: Update label.
    """
    label = label.replace('_', '')
    return label[:3] + "-" + label[3:] 
    

def plot_pvalue_distributions(pvalue_file : str):
    """Creates dot plot for distributions of pvalues of MWU on PSA value and Gleason score for
    different competitors.

    Args:
        pvalue_file (str): Name of file storing pvalues for cohorts, tools and phenotypes.
    """
    label_fontsize = 58
    
    pval_df = pd.read_csv(pvalue_file, sep='\t')
    methods_list = ['PARADYS_DysRegNet', 'PNC', 'HITnDrive', 'PRODIGY']
    ordered_clusters = ['2.0_2.1', '3.0_3.1', '3.0_3.2', '3.1_3.2', '4.0_4.1', '4.0_4.2', '4.0_4.3', '4.1_4.2', '4.1_4.3', '4.2_4.3']
    
    # Extract parts that correspond to PRAD cohort, Gleason score.
    psa_df = pval_df[(pval_df['phenotype']=='psa_value') & (pval_df['cohort']=='PRAD') & (pval_df['method'].isin(methods_list)) & (pval_df['cluster'].isin(ordered_clusters))]
    gleason_df = pval_df[(pval_df['phenotype']=='gleason_score') & (pval_df['cohort']=='PRAD') & (pval_df['method'].isin(methods_list)) & (pval_df['cluster'].isin(ordered_clusters))]
    # Transform naming of cluster from e.g. 2.0_2.1 to 2.0 vs. 2.1.
    psa_df['cluster']= psa_df['cluster'].transform(change_cluster_labels)
    gleason_df['cluster'] = gleason_df['cluster'].transform(change_cluster_labels)
    cluster_combis = set(psa_df['cluster'])
    
    # Transform PSA dataframe such that methods are used as columns with pvalue entries.
    psa_transformed = pd.DataFrame(index=sorted(list(cluster_combis)), columns=methods_list)
    gleason_transformed =pd.DataFrame(index=sorted(list(cluster_combis)), columns=methods_list)
    
    for _, row in psa_df.iterrows():
        cluster = row['cluster']
        pval = row['p_value']
        method = row['method']
        # Insert pvalue into new dataframe.
        psa_transformed.loc[cluster, method] = pval
    
    psa_transformed['cluster']=psa_transformed.index
    
    # Transform Gleason score dataframe.
    for _, row in gleason_df.iterrows():
        cluster = row['cluster']
        pval = row['p_value']
        method = row['method']
        gleason_transformed.loc[cluster, method] = pval
    
    gleason_transformed['cluster'] = gleason_transformed.index
        
    psa_transformed['Phenotype']='PSA value'
    gleason_transformed['Phenotype']='Gleason score'
    combined_df = pd.concat([psa_transformed, gleason_transformed], ignore_index=True)
        
    
    # Create dot plots using pvalues as dots.
    # sns.set_theme(style="whitegrid")

    # Make the pairgrid.
    g = sns.PairGrid(combined_df, hue="Phenotype",
                     x_vars=combined_df.columns[:4], y_vars=['cluster'], hue_order=['PSA value', 'Gleason score'],
                     height=20, aspect=.25)
    #g.fig.set_size_inches(70,40)
    g.tick_params(axis='y', labelsize=label_fontsize)
    
    # Draw a dot plot using the stripplot function
    g.map(sns.stripplot, size=35, orient="h", jitter=True,
          linewidth=0.5, edgecolor="w")
    
    # Use the same x axis limits on all columns and add better labels
    g.set(xlim=(-0.1, 1.1), ylabel="")
    g.add_legend(title='', fontsize=label_fontsize)

    # Use semantically meaningful titles for the columns
    titles = ["PARADYS (DysRegNet)", "PNC", "HITnDrive",
              "PRODIGY"]

    for ax, title in zip(g.axes.flat, titles):

        # Make the grid horizontal instead of vertical
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)
        ax.set_xlabel("p-value", fontsize=label_fontsize)
        ax.xaxis.set_tick_params(labelsize=label_fontsize)
        ax.axvline(x = 0.05, color = 'black', linestyle='--', linewidth=4.0)
        
        # Set a different title for each axes
        ax.set_title(title, {'fontsize' : label_fontsize})

    
    sns.despine(left=True, bottom=True)
    plt.savefig('phenotype_pvals.pdf', format='pdf')


if __name__=='__main__':
    
    # Input paths for Gleason/OS/Xcell plots.
    cluster_file_path = "../results/clustering/clusters/clusters_PARADYS_DysRegNet_PRAD.csv"
    phenotype_file = "../data/PRAD_phenotypes.csv"
    survival_file = "../data/PRAD_survival.txt"
    xcell_results_path = "../results/clustering/xcell/"
    #plot_gleason_os_xcell(cluster_file_path, phenotype_file, survival_file, xcell_results_path)
    
    # Input paths for FDR plots.
    #plot_fdrs()
    
    # Input files for top K drivers extraction.
    k = 20
    drivers_file = '../results/clustering/driver_frequencies/cluster_drivers_PRODIGY_PRAD.csv'
    clusters = ['4.0', '4.1', '4.2', '4.3']
    output_path = './top_drivers_prodigy/'
    # get_top_k_drivers(k, drivers_file, clusters, output_path) 
    
    # Plot pvalues of DIGEST on top 20 drivers of PARADYS, PRODIGY, HITNDRIVE.
    digest_path = "../results/clustering/digest/"
    tools = ['dysregnet', 'hitndrive', 'prodigy']
    clusters = ['4.0', '4.1', '4.2', '4.3']
    # plot_digest_pvalues(digest_path, tools, clusters)
    
    # Plot pvalues of PSA value and Gleason score for several tools.
    pvals_file = '../results/clustering/phenotype_pvalues/p_values_phenotypes.csv'
    #plot_pvalue_distributions(pvals_file)

