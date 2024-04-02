import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd
import numpy as np
import seaborn as sns
import scipy as sc
    
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
    axs_list[0].set_ylabel('Gleason score')
    axs_list[0].set_xticks([y + 1 for y in range(len(gleason_data))], labels=['4.0', '4.1', '4.2', '4.3'])
    axs_list[0].set_xlabel('Cluster')
    axs_list[0].grid(True, axis='y')
    
    # Add OS time subplot.
    axs_list[1].boxplot(os_data, patch_artist=True, medianprops=dict(color='firebrick', linewidth=2.3), boxprops=dict(facecolor='mistyrose'))
    axs_list[1].set_ylabel('Survival time [days]')
    axs_list[1].set_xticks([y + 1 for y in range(len(os_data))], labels=['4.0', '4.1', '4.2', '4.3'])
    axs_list[1].set_xlabel('Cluster')
    axs_list[1].grid(True, axis='y')
    
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
    
    s = sns.heatmap(means_df, annot=True, fmt='.1f', ax=axs_list[2], cbar_kws={'label': 'Mean xCell score'})
    s.set_xlabel('Cluster')
    s.set_ylabel('Cell type')
    
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
    s = sns.heatmap(pvalues_signif, annot=True, fmt='.2f', vmin=0, vmax=17, ax=axs_list[3], cbar_kws={'label': '-log10(p-value)'})
    s.set_xlabel('Target cluster')
    #s.set_ylabel('Cell type')
    
    
    # Add mosaic labels above plots.
    for label, ax in axs.items():
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='xx-large', va='bottom', fontfamily='sans-serif', fontweight='bold')
    
    fig.tight_layout()
    plt.savefig("gleason_os_xcell.pdf", format='pdf')
    
    
def plot_fdrs():
    """Plot FDR results from all used tools and input networks.
    """
    # Create one large mosaic plot.
    fig, axs = plt.subplot_mosaic([['A', 'B'], ['C', 'D']])
    fig.set_size_inches(12, 12)
    axs_list = list(axs.values())
    
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
    axs_list[0].set_ylabel("FDR")
    axs_list[0].set_xlabel("Rate of decoy mutations")
    axs_list[0].set_axisbelow(True)
    axs_list[0].grid(axis='y')
    
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
    axs_list[3].set_ylabel("FDR")
    axs_list[3].set_xlabel("Rate of decoy mutations")
    axs_list[3].set_axisbelow(True)
    axs_list[3].grid(axis='y')
    
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
    axs_list[2].set_ylabel("FDR")
    axs_list[2].set_xlabel("Rate of decoy mutations")
    axs_list[2].set_axisbelow(True)
    axs_list[2].grid(axis='y')
    
    ### Add PARADYS with DysRegNet_NoDir networks FDRS on top right.
    cohorts = ['prad', 'coad']
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
    axs_list[1].set_ylabel("FDR")
    axs_list[1].set_xlabel("Rate of decoy mutations")
    axs_list[1].set_axisbelow(True)
    axs_list[1].grid(axis='y')
    
    # Add mosaic labels above plots.
    for label, ax in axs.items():
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
                fontsize='x-small', va='bottom', fontfamily='sans-serif', fontweight='bold')
    
    fig.tight_layout()
    plt.savefig("fdrs.pdf", format='pdf')
    

def plot_fdrs_paradys():
    """Plots violinplot of FDR of PARADYS colored by cohort and categorized by % simulated.
    """
    # Plot FDR files in multiple boxplots.
    file_path = '../results/fdrs/PARADYS_DysRegNet/fdrs_dysregnet_all_cohorts.csv'
    df = pd.read_csv(file_path, sep=';')
    df.rename(columns={"cohort" : "Cohort"}, inplace=True)
    #plt.ylim(0,0.04)
    sns.set_theme(style='darkgrid')
    sns.set(font_scale=3.0)
    axs = sns.violinplot(x = df['decoy_rate'],
                y = df['FDR'], 
                hue = df['Cohort'],
                cut=0)
    # axs.set_title("FDRs of PARADYS using DysRegNet networks")
    axs.set_ylabel("FDR")
    axs.set_xlabel("Rate of decoy mutations")
    plt.show()
    
def plot_fdrs_prodigy_new():
    # Plot FDR files in multiple boxplots.
    file_path = '../results/fdrs/PRODIGY/fdrs_prodigy.csv'
    df = pd.read_csv(file_path, sep=',')
    rate_dict = {25 : '0.25', 5 : '0.50', 75 : '0.75'}
    df.replace({"Rate": rate_dict}, inplace=True)
    
    sns.set_theme(style='darkgrid')
    sns.set(font_scale=3.0)
    axs = sns.violinplot(x = df['Rate'],
                y = df['FDR'], 
                hue = df['Cohort'],
                cut=0, linewidth=2.5)
    # axs.set_title("FDRs of PARADYS using DysRegNet networks")
    axs.legend(loc='upper right')
    axs.set_ylabel("FDR")
    axs.set_xlabel("Rate of decoy mutations")
    plt.show()
    
    
def plots_fdrs_prodigy():
    """Combines PRODIGY FDR results and plots them together in violinplot.
    """
    cohorts = ['BRCA', 'COAD', 'PRAD']
    rates = ['025', '05', '075']
    rates_float = [0.25, 0.5, 0.75]
    
    df = pd.DataFrame(columns=["FDR", "Rate", "Cohort"])
     
    for type in cohorts:
        for rate, rate_float in zip(rates, rates_float):
            file_path = f'FPRs/False_mutations_PRODIGY/FDRs/{type}_{rate}_1_fdr_patients.csv'
            sub_df = pd.read_csv(file_path, sep=',')
            sub_df.drop(columns=['patient'], inplace=True)
            sub_df['Rate'] = rate_float
            sub_df['Cohort'] = type
            df = pd.concat([df, sub_df], axis=0)
    
    print(df)
    sns.set_style('whitegrid')
    sns.set(font_scale=3.0)
    axs = sns.violinplot(x = df['Rate'],
                y = df['FDR'], 
                hue = df['Cohort'],
                cut=0)
    # axs.set_title("FDRs of PARADYS using DysRegNet networks")
    
    axs.set_ylabel("FDR")
    axs.set_xlabel("Rate of decoy mutations")
    
    plt.show()

    

if __name__=='__main__':
    
    # Input paths for Gleason/OS/Xcell plots.
    cluster_file_path = "../results/clustering/clusters/clusters_PARADYS_DysRegNet_PRAD.csv"
    phenotype_file = "../data/PRAD_phenotypes.csv"
    survival_file = "../data/PRAD_survival.txt"
    xcell_results_path = "../results/clustering/xcell/"
    plot_gleason_os_xcell(cluster_file_path, phenotype_file, survival_file, xcell_results_path)
    
    # Input paths for FDR plots.
    #plot_fdrs()
    
    

