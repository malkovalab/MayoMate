# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import os
import pandas as pd
from typing import Literal
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import betterbeeswarm
from mayo.settings import genotype_dict as genotype_dict_master
from mayo.settings import config as cfg
logging.basicConfig(
    handlers=[logging.FileHandler(
        filename="outputs/logs/cluster_coordination.log", 
        encoding='utf-8', 
        mode='w')],
    level=logging.INFO,
    format='%(levelname)s:.%(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

genotype_dict = genotype_dict_master[cfg["genotype_dict"]]

from mayo.settings import haploid_samples as haploid_samples

def check_canonnical_row(row: pd.Series) -> bool:

    #if the reference allele is "C" then the tumor allele must be "T"
    if row['Reference_Allele'] == "C":
        if row['Tumor_Seq_Allele2'] == "T":
            return False
        else:
            return True
    #if the reference allele is "G" then the tumor allele must be "A"
    elif row['Reference_Allele'] == "G":
        if row['Tumor_Seq_Allele2'] == "A":
            return False
        else:
            return True
    
    #else its an error
    else:
        return True

def check_canonnical(df: pd.DataFrame) -> bool:
    #iterate through the rows of the df_cluster
    non_canonical = False
    for index, row in df.iterrows():
        non_canonical = check_canonnical_row(row)
        if non_canonical:
            break
    return non_canonical

def add_genotype_column_PMACD(df: pd.DataFrame, genotype_dict: dict) -> pd.DataFrame:
    df['genotype'] = df['Tumor_Sample_Barcode'].map(genotype_dict)
    return df

def add_genotype_column_JT(df: pd.DataFrame, genotype_dict: dict) -> pd.DataFrame:
    df['genotype'] = df['Sample'].map(genotype_dict)
    return df

def countSwitches(df: pd.DataFrame) -> int:
    #set the current_Reference_Allele to the first row's Reference_Allele
    current_Reference_Allele = df.iloc[0]['Reference_Allele']
    switches = 0
    for index, row in df.iterrows():
        if current_Reference_Allele != row['Reference_Allele']:
            switches += 1
            current_Reference_Allele = row['Reference_Allele']

    return switches

def get_cluster_length(df: pd.DataFrame) -> int:
    #cluster length is defined as the absolute difference in the first and last position of the cluster
    first_position = df.iloc[0]['Start_position']
    last_position = df.iloc[-1]['Start_position']
    return abs(first_position - last_position)

def get_cluster_p_value(df: pd.DataFrame) -> float:
    #returns the value of the p-value column
    return df.iloc[0]['Cluster_Pvalue']

def get_cluster_dict(
        df: pd.DataFrame, 
        genotype, 
        filter_set: list | None = None, 
        kind: Literal["in", "not_in"] = "in",
        remove_rtg: bool = True
        ) -> tuple[dict, pd.DataFrame, dict]:
    
    from mayo.settings import rtg_list  # Assuming this is a valid import
    
    # Count total clusters and max cluster ID
    total_cluster_ids = len(df['Dataset_Cluster_ID'].value_counts())
    max_cluster_id = df['Dataset_Cluster_ID'].max()

    cluster_dict = {}
    clusters_df = pd.DataFrame(columns=['dataset_cluster_id', 'cluster_length', 'cluster_type'])

    for i in range(max_cluster_id):
        df_cluster = df[df['Dataset_Cluster_ID'] == i + 1]
        df_cluster = df_cluster[[
            "Reference_Allele", 
            "Tumor_Seq_Allele2", 
            "genotype", 
            "Tumor_Sample_Barcode", 
            "Start_position", 
            "Cluster_Pvalue", 
            "Dataset_Cluster_ID"]]
        df_cluster = df_cluster[df_cluster['genotype'] == genotype]

        if remove_rtg:
            df_cluster = df_cluster[~df_cluster['Tumor_Sample_Barcode'].isin(rtg_list)]


        if filter_set is not None:
            if kind == "in":
                df_cluster = df_cluster[df_cluster['Tumor_Sample_Barcode'].isin(filter_set)]
            elif kind == "not_in":
                df_cluster = df_cluster[~df_cluster['Tumor_Sample_Barcode'].isin(filter_set)]
            else:
                raise ValueError("kind must be either 'in' or 'not_in'")
            
        cluster_mutation_count = len(df_cluster)
        if cluster_mutation_count == 0:
            continue

        cluster_length = get_cluster_length(df_cluster)
        cluster_p_val = get_cluster_p_value(df_cluster)

        # Determine cluster type
        non_canonical = check_canonnical(df_cluster)
        switches = countSwitches(df_cluster)
        cluster_type = None

        if cluster_mutation_count == 2:
            cluster_type = "just_2_mutations"
        elif non_canonical:
            cluster_type = "non_canonical"
        elif (len(df_cluster['Reference_Allele'].value_counts()) == 1) and \
             (len(df_cluster['Tumor_Seq_Allele2'].value_counts()) == 1):
            allele = df_cluster.iloc[0]['Reference_Allele']
            if allele == "C":
                cluster_type = "non_switching_C"
            elif allele == "G":
                cluster_type = "non_switching_G"
        elif switches == 1:
            value_counts = df_cluster['Reference_Allele'].value_counts()
            starting_allele = df_cluster.iloc[0]['Reference_Allele']
            second_allele = df_cluster.iloc[1]['Reference_Allele']

            if value_counts[1] == 1:
                if starting_allele == "C":
                    cluster_type = "switching_once_end_CG" if second_allele == "C" else "switching_once_beginning_CG"
                elif starting_allele == "G":
                    cluster_type = "switching_once_end_GC" if second_allele == "G" else "switching_once_beginning_GC"
            else:
                cluster_type = "switching_once_middle_CG" if starting_allele == "C" else "switching_once_middle_GC"
        elif switches > 1:
            cluster_type = "multiple_switching"

        if cluster_type:
            #concatenate the cluster_type to the clusters_df
            clusters_df = pd.concat([clusters_df, pd.DataFrame({
                'dataset_cluster_id': [i + 1],
                'cluster_length': [cluster_length],
                'cluster_mutations': [cluster_mutation_count],
                'cluster_type': [cluster_type],
                'genotype': [genotype],
                'cluster_p_val': [cluster_p_val]
            })], ignore_index=True)
            cluster_dict[i + 1] = cluster_type
        else:
            logging.error(f"Cluster ID {i + 1} is not classified")
            print(f"Cluster ID {i + 1} is not classified")

    return clusters_df, cluster_dict

def plot_stacked_bar_chart_clusters(
        master_df_bar: pd.DataFrame,
        title: str = "Proportion of Cluster Types",
        order = ["ung1∆", "ung1∆NAT", "exo1-nd", "pol32∆", "exo1-nd pol32∆"],
        save_name: str|None = None) -> None:
    """
    This function plots a stacked bar chart of cluster types by genotype.

    Parameters:
    master_df_bar (pd.DataFrame): The DataFrame containing the data to be plotted.
    save_name (str|None): The name of the file where the plot will be saved. 
                          If None, the plot will not be saved. Default is None.

    Returns:
    None: The function displays the plot and does not return any value.
    """
    
    #make a colormap for the cluster types
    colors = [
        "tab:blue", 
        "tab:green", 
        "tab:red", 
        "tab:orange", 
        "tab:gray"]
    
    colors = [
        "#4d72b0", 
        "#55a868", 
        "#c44e52", 
        "#dd8452",
        "#8c8c8c"]
    
    # colors = [
    #     "tab:blue", 
    #     "tab:red", 
    #     "tab:purple", 
    #     "tab:brown", 
    #     "tab:pink", 
    #     "tab:gray", 
    #     "tab:green", 
    #     "tab:orange", 
    #     "tab:cyan", 
    #     "tab:olive"]

    #make a stacked bar chart, cluster type is the color, % of clusters is the y axis
    #if exo1-ndpol32∆ rename to exo1-nd pol32∆
    #master_df_bar = master_df_bar.rename(index={"exo1-ndpol32∆": "exo1-nd pol32∆", "exo1-ndsgs1∆C": "exo1-nd sgs1∆C"})
    categories = ['non-switching', 'switching once CG', 'switching once GC', 'multiple switching', 'non canonical']
    #if column is not present, add it with 0s
    for category in categories:
        if category not in master_df_bar.columns:
            master_df_bar[category] = 0
    master_df_bar = master_df_bar[categories]
    master_df_bar = master_df_bar.reindex(order)
    sns.set_style('whitegrid')
    ax = master_df_bar.plot(kind='bar', stacked=True, figsize=(10,8), fontsize=20, rot=45, color=colors, legend=False, grid=True, width=0.75)
    for label in ax.get_xticklabels():
        label.set_ha('right')
    plt.title(title, fontsize=30, fontweight='bold')
    plt.xlabel("Genotype", fontsize=20, fontweight='bold')
    plt.ylabel("% of Clusters", fontsize=20, fontweight='bold')
    plt.xticks(fontstyle='italic', fontsize=20)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    plt.tight_layout()
    if save_name is not None:
        plt.savefig(save_name, dpi=300)
    #plt.show()
    plt.close()
 
def plot_cluster_length(master_df: pd.DataFrame, save_name: str|None = None) -> None:
    """
    This function plots a stacked bar chart of cluster types by genotype.

    Parameters:
    master_df_bar (pd.DataFrame): The DataFrame containing the data to be plotted.
    save_name (str|None): The name of the file where the plot will be saved. 
                          If None, the plot will not be saved. Default is None.

    Returns:
    None: The function displays the plot and does not return any value.
    """
    master_df = master_df.copy()
    master_df["genotype"] = master_df["genotype"].replace({
        "ung1∆": "ung1∆", 
        "exo1-nd": "exo1-nd", 
        "pol32∆": "pol32∆", 
        "exo1-ndpol32∆": "exo1-nd pol32∆",
        "sgs1∆C": "sgs1∆C",
        "exo1-ndsgs1∆C": "exo1-nd sgs1∆C",
        })
    master_df['cluster_length'] = master_df['cluster_length']/1000

    order = [
        "non-switching", 
        "switching once CG", 
        "switching once GC", 
        "multiple switching", 
        "non canonical"]
    
    order_genotype = [
        "ung1∆", 
        "ung1∆NAT", 
        "exo1-nd", 
        "pol32∆", 
        "exo1-nd pol32∆",
        "sgs1∆C", 
        "exo1-nd sgs1∆C"]
    
    #make a colormap for the cluster types
    pallette = {
        "ung1∆": "tab:blue", 
        "ung1∆NAT": "tab:orange", 
        "exo1-nd": "tab:red", 
        "pol32∆": "tab:green", 
        "exo1-nd pol32∆": "tab:purple",
        "sgs1∆C": "tab:pink", 
        "exo1-nd sgs1∆C": "tab:gray"}

    #make a series of violinplots for each cluster type by each genotype
    fig, ax = plt.subplots(figsize=(45, 20))
    sns.violinplot(x='cluster_type', y='cluster_length', hue='genotype', data=master_df, ax=ax, inner="quartile", density_norm="width", order=order, hue_order=order_genotype, palette=pallette)
    sns.swarmplot(x='cluster_type', y='cluster_length', hue='genotype', data=master_df, ax=ax, alpha=0.85, size=4, dodge=True, palette="dark:black", order=order, legend=False)
    plt.xticks(rotation=0, fontsize=30)
    ax.set_yticks(range(0, 71, 5))
    plt.yticks(fontsize=30)
    plt.subplots_adjust(bottom=0.15)
    ax.set_title("Cluster Lengths by Cluster Type and Genotype", fontsize=65)
    ax.set_xlabel("Cluster Type", fontsize=45)
    ax.set_ylabel("Cluster Length (kb)", fontsize=45)
    #place the legend outside of the plot
    ax.legend(fontsize=30, title_fontsize=45, title="Genotype", loc="upper left", bbox_to_anchor=(1, 1))
    for t, l in zip(ax.legend_.texts, ["ung1∆", "ung1∆NAT", "exo1-nd", "pol32∆", "exo1-nd pol32∆", "sgs1∆C", "exo1-nd sgs1∆C"]):
        t.set_text(l)
        t.set_fontsize(30)
        t.set_fontstyle("italic")

    plt.ylim(0, 65)
    plt.tight_layout()
    if save_name is not None:
        plt.savefig(save_name, dpi=300)
    #plt.show()
    plt.close()

def remove_not_assocaited_clusters(
        df: pd.DataFrame, 
        df_association: pd.DataFrame,
        association_feature: Literal["GC", "CO", "both"], 
        inverse=False) -> pd.DataFrame:
    """remove clusters that are not associated with a switch event"""

    #merge the two dataframes
    logging.info(f"Merging the two dataframes... len(df): {len(df)}, len(df_association): {len(df_association)}")
    df = df.merge(df_association, left_on="Dataset_Cluster_ID", right_on="Cluster_ID_", how="left")

    #count the number of cluster IDs by unique Dataset_Cluster_ID
    unique_cluster_ids = df["Dataset_Cluster_ID"].value_counts()
    logging.info(f"Unique Dataset_Cluster_ID in df: {len(unique_cluster_ids)}")

    #remove the rows that are not associated with a switch event
    feature = association_feature
    if association_feature == "both":
        feature = "Clustered"

    if inverse:
        df = df[df[f"{feature}_Proximal_"] == False]
        logging.info("Inverse: removing clusters that are associated with a switch event")
    else:
        df = df[df[f"{feature}_Proximal_"] == True]
        logging.info("Removing clusters that are not associated with a switch event")

    unique_cluster_ids_after = df["Dataset_Cluster_ID"].value_counts()

    logging.info(f"Removed {len(unique_cluster_ids) - len(unique_cluster_ids_after)} clusters that are not associated with a switch event. There are {len(unique_cluster_ids_after)} clusters left.")

    return df

def filter_clusters_by_association_with_rec(
        cluster_caller: str,
        cluster_distance: int,
        df: pd.DataFrame, 
        cluster_association: str, 
        association_feature: Literal["GC", "CO", "both"],
        cluster_data_dir: str,
        snp_support: int = 1,
        ) -> pd.DataFrame:
    
    if cluster_association == "all_clusters":
        pass

    elif cluster_association == "not_associated" or cluster_association == "rec_associated":
        df_association = df_association=pd.read_csv(f'outputs/associations/{cluster_caller}/{snp_support}snp/master_clusters_association_{cluster_data_dir}_{snp_support}snp_{cluster_distance}.csv')

        if cluster_association == "rec_associated":
            df = remove_not_assocaited_clusters(df, df_association, association_feature, inverse=False)
        elif cluster_association == "not_associated":
            df = remove_not_assocaited_clusters(df, df_association, association_feature, inverse=True)

    else:
        logging.error(f"Cluster_association: '{cluster_association}' not found. Choose 'rec_associated', 'not_associated' or 'all_clusters'")
        raise ValueError
    
    return df

def cluster_type_mann_whitney(master_df: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import mannwhitneyu
    import numpy as np

    test_df = pd.DataFrame()
    master_df['cluster_length'] = master_df['cluster_length'].astype(float)

    for cluster_type in master_df['cluster_type'].unique():
        for genotype in master_df['genotype'].unique():
            df1 = master_df[(master_df['genotype'] == genotype) & (master_df['cluster_type'] == cluster_type)]
            df2 = master_df[(master_df['genotype'] == 'ung1∆') & (master_df['cluster_type'] == cluster_type)]
            
            mean1, median1 = df1['cluster_length'].mean(), df1['cluster_length'].median()
            mean2, median2 = df2['cluster_length'].mean(), df2['cluster_length'].median()

            if len(df1) > 0 and len(df2) > 0:

                stat, p = mannwhitneyu(df1['cluster_length'], df2['cluster_length'])
                logging.info(f"cluster_type: {cluster_type:18}, genotype: {genotype:15}, p: {p:.6f}, mean: {mean2:.4f}kb vs. {mean1:.4f}kb, median: {median2:.4f}kb vs. {median1:.4f}kb")

                test_df = pd.concat([
                    test_df, pd.DataFrame({
                        'cluster_type': cluster_type, 
                        'genotype': genotype, 
                        'p_val': p}, 
                    index=[0])],
                    ignore_index=True)

            else:

                logging.info(f"cluster_type: {cluster_type:18}, genotype: {genotype:15}, p: NA, {mean2:.4f}kb vs. {genotype} mean: {mean1:.4f}kb, wt median: {median2:.4f}kb vs. {median1:.4f}kb")
                test_df = pd.concat([
                    test_df, 
                    pd.DataFrame({
                        'cluster_type': cluster_type, 
                        'genotype': genotype, 
                        'p_val': np.nan},
                    index=[0])], 
                    ignore_index=True)
    
    return test_df


if __name__ == "__main__":

    snp_support_vals: list[int] = [
        1, 
        #2
        ]
    
    association_feature: Literal["GC", "CO", "both"] = "both"

    cluster_caller: str = "JT"
    cluster_distance: int = 10000

    sample_sets: list[str] = [
        'allploid',
        # 'haploid',
        # 'aneuploid'
        ]
    
    cluster_associations: list[str] = [
        "all_clusters",
        "rec_associated",
        "not_associated"
        ]

    cluster_data_dirs: list[str] = [
        # '05',
        # '01',
        # '005',
        # '001_def',
        # '0005',
        # '0001'
        0.0001,
        # 0.01,
        # 0.001,
        # 0.00001,
        # 0.000001
        ]
    
    GENOTYPES: list[str] = [
        # "ung1∆ non-selected",
        "ung1∆",
        #"ung1∆NAT",
        "exo1-nd",
        "pol32∆",
        "exo1-ndpol32∆",
        "sgs1∆C",
        "exo1-ndsgs1∆C",
        # WT
        ]
    
    meta_master_df = pd.DataFrame()
    meta_test_df = pd.DataFrame()
    for snp_support in snp_support_vals:
        for cluster_association in cluster_associations:
            for sample_set in sample_sets:
                for cluster_data_dir in cluster_data_dirs:
                    logging.info(f"Working on {snp_support}snp support, {cluster_association}, {sample_set}, pval: {cluster_data_dir}...")
                    
                    if association_feature != "both" and cluster_association == "all_clusters":
                        logging.info(f"association_feature: {association_feature} and cluster_association: {cluster_association} are not compatible. Skipping...")
                        continue

                    if cluster_caller == "PMACD":
                        CLUSTER_DATA_FILE = f"data/cluster_files/clusters_diff_thresholds/master_p{cluster_data_dir}/1_master_para_bed_sorted_anz5.txt"
                        df = pd.read_csv(CLUSTER_DATA_FILE, sep='\t')
                        df = add_genotype_column_PMACD(df, genotype_dict)
                    
                    elif cluster_caller == "JT":
                        CLUSTER_DATA_FILE = f"data/new_clusters/clustered_mutations/clustered_SNVs_JT_10000_{cluster_data_dir}.tsv"
                        df = pd.read_csv(CLUSTER_DATA_FILE, sep='\t', encoding='utf-16')
                        df = add_genotype_column_JT(df, genotype_dict)
                        df["Cluster_Pvalue"] = cluster_data_dir
                        df = df.rename(columns={
                            "Reference": "Reference_Allele",
                            "Allele": "Tumor_Seq_Allele2",
                            "Sample": "Tumor_Sample_Barcode",
                            "Region": "Start_position"
                            })
                        
                    df = filter_clusters_by_association_with_rec(
                        cluster_caller,
                        cluster_distance,
                        df, 
                        cluster_association, 
                        association_feature,
                        cluster_data_dir,
                        snp_support=snp_support)

                    #filter data, include only those where Dataset_Cluster_ID is present
                    df = df[df['Dataset_Cluster_ID'].notna()]

                    #set the df['Dataset_Cluster_ID'] column to be an integer
                    df['Dataset_Cluster_ID'] = df['Dataset_Cluster_ID'].astype(int)

                    master_df = pd.DataFrame()
                    master_df_bar = pd.DataFrame()

                    for GENOTYPE in GENOTYPES:

                        if sample_set == 'allploid':
                            cluster_df, cluster_dict  = get_cluster_dict(df, GENOTYPE)
                        elif sample_set == 'haploid':
                            cluster_df, cluster_dict = get_cluster_dict(df, GENOTYPE, filter_set=haploid_samples, kind="in")
                        elif sample_set == 'aneuploid':
                            cluster_df, cluster_dict = get_cluster_dict(df, GENOTYPE, filter_set=haploid_samples, kind="not_in")
                        else:
                            logging.error(f"Sample Set: {sample_set} not found. Choose between 'haploid' or 'aneuploid' or 'all'")
                            raise ValueError

                        master_df = pd.concat([master_df, cluster_df], ignore_index=True)

                    #rename the columns
                    master_df_bar = master_df.copy()
                    rename_dict = {
                        'non_switching_C':             'C-tract',
                        'non_switching_G':             'G-tract',
                        'switching_once_beginning_CG': 'single C to G-tract',
                        'switching_once_beginning_GC': 'single G to C-tract',
                        'switching_once_end_CG':       'C-tract to single G',
                        'switching_once_end_GC':       'G-tract to single C',
                        'switching_once_middle_CG':    'C-tract to G-tract',
                        'switching_once_middle_GC':    'G-tract to C-tract',
                        'multiple_switching':          'multiple switching',
                        'non_canonical':               'non canonical'}
                    master_df_bar['cluster_type'] = master_df_bar['cluster_type'].replace(rename_dict)
                    master_df_bar = master_df_bar.groupby(['genotype', 'cluster_type'])['cluster_length'].count().reset_index().copy()
                    master_df_bar = master_df_bar.pivot(index='genotype', columns='cluster_type', values='cluster_length')
                    master_df_bar = master_df_bar.fillna(0)

                    master_df_bar["non-switching"] = master_df_bar["C-tract"] + master_df_bar["G-tract"]
                    #if column is not present, add it with 0
                    if "single C to G-tract" not in master_df_bar.columns:
                        master_df_bar["single C to G-tract"] = 0
                    if "single G to C-tract" not in master_df_bar.columns:
                        master_df_bar["single G to C-tract"] = 0
                    if "C-tract to G-tract" not in master_df_bar.columns:
                        master_df_bar["C-tract to G-tract"] = 0
                    if "G-tract to C-tract" not in master_df_bar.columns:
                        master_df_bar["G-tract to C-tract"] = 0
                    if "C-tract to single G" not in master_df_bar.columns:
                        master_df_bar["C-tract to single G"] = 0
                    if "G-tract to single C" not in master_df_bar.columns:
                        master_df_bar["G-tract to single C"] = 0


                    master_df_bar["switching once CG"] = master_df_bar["single C to G-tract"] + master_df_bar["C-tract to G-tract"] + master_df_bar["C-tract to single G"]
                    master_df_bar["switching once GC"] = master_df_bar["single G to C-tract"] + master_df_bar["G-tract to C-tract"] + master_df_bar["G-tract to single C"]
                    master_df_bar = master_df_bar.drop(columns=[
                        'C-tract', 
                        'G-tract', 
                        'single C to G-tract', 
                        'C-tract to G-tract', 
                        'C-tract to single G', 
                        'single G to C-tract', 
                        'G-tract to C-tract', 
                        'G-tract to single C'])
                    
                    master_df_bar = master_df_bar.fillna(0)
                    categories = ['non-switching', 'switching once CG', 'switching once GC', 'multiple switching', 'non canonical']

                    proportion_master_df_bar = master_df_bar.div(master_df_bar.sum(axis=1), axis=0) * 100
                    proportion_master_df_bar = proportion_master_df_bar.fillna(0)
                    for category in categories:
                        if category not in master_df_bar.columns:
                            proportion_master_df_bar[category] = 0
                    proportion_master_df_bar = proportion_master_df_bar[categories]

                    print("propotion of cluster types (number) by genotype")
                    print(proportion_master_df_bar)
                    logging.info("propotion of cluster types (number) by genotype")
                    logging.info(proportion_master_df_bar)

                    plot_stacked_bar_chart_clusters(
                        proportion_master_df_bar, 
                        save_name=f"figures/cluster_types_by_pval/{cluster_caller}/cluster_types_by_genotype_stacked_bar_{sample_set}_{cluster_association}_{association_feature}_{cluster_data_dir}_{snp_support}snp_proportion_number_of_clusters.png",
                        order=GENOTYPES)

                    master_df['cluster_type'] = master_df['cluster_type'].replace({
                        'non_switching_C':             'C-tract',
                        'non_switching_G':             'G-tract',
                        'switching_once_beginning_CG': 'single C G-tract',
                        'switching_once_beginning_GC': 'single G C-tract',
                        'switching_once_end_CG':       'C-tract single G',
                        'switching_once_end_GC':       'G-tract single C',
                        'switching_once_middle_CG':    'C-tract G-tract',
                        'switching_once_middle_GC':    'G-tract C-tract',
                        'multiple_switching':          'multiple switching',
                        'non_canonical':               'non canonical'})
                    
                    master_df['cluster_type'] = master_df['cluster_type'].replace({
                        'C-tract':            'non-switching',
                        'G-tract':            'non-switching',

                        'single C G-tract':   'switching once CG',
                        'C-tract G-tract':    'switching once CG',
                        'C-tract single G':   'switching once CG',

                        'single G C-tract':   'switching once GC',
                        'G-tract C-tract':    'switching once GC',
                        'G-tract single C':   'switching once GC',

                        'multiple switching': 'multiple switching',
                        'non canonical':      'non canonical'})
                    
                    plot_cluster_length(
                        master_df,
                        save_name=f"figures/cluster_length_violin_by_pval/{cluster_caller}/cluster_length_violinplot_{sample_set}_{cluster_association}_{association_feature}_{cluster_data_dir}_{snp_support}snp.png")

                    proportion_ssdna_df = master_df.groupby(['genotype', 'cluster_type'])['cluster_length'].sum().reset_index().copy()
                    proportion_ssdna_df = proportion_ssdna_df.pivot(index='genotype', columns='cluster_type', values='cluster_length')
                    proportion_ssdna_df = proportion_ssdna_df.div(proportion_ssdna_df.sum(axis=1), axis=0) * 100
                    proportion_ssdna_df = proportion_ssdna_df.fillna(0)

                    print("proportion of ssDNA contributed by cluster type by genotype")
                    print(proportion_ssdna_df)
                    logging.info("proportion of ssDNA contributed by cluster type by genotype")
                    logging.info(proportion_ssdna_df)

                    plot_stacked_bar_chart_clusters(
                        proportion_ssdna_df, 
                        save_name=f"figures/cluster_types_by_pval/{cluster_caller}/cluster_types_by_genotype_stacked_bar_{sample_set}_{cluster_association}_{association_feature}_{cluster_data_dir}_{snp_support}snp_proportion_ssdna.png",
                        title="ssDNA Contribution by Cluster Type",
                        order=GENOTYPES)
                    
                    master_df["genotype"] = master_df["genotype"].replace({
                        "ung1∆":         "ung1∆", 
                        "exo1-nd":       "exo1-nd", 
                        "pol32∆":        "pol32∆", 
                        "exo1-ndpol32∆": "exo1-nd pol32∆"})
                    
                    #master_df = master_df[master_df['cluster_length'].notna()]
                    master_df['cluster_caller'] = cluster_caller
                    master_df['cluster_distance'] = cluster_distance
                    master_df["sample_set"] = sample_set
                    master_df["cluster_association"] = cluster_association
                    master_df["association_feature"] = association_feature
                    master_df["snp_support"] = snp_support
                    meta_master_df = pd.concat([meta_master_df, master_df], ignore_index=True)

                    test_df = cluster_type_mann_whitney(master_df)
                    test_df["cluster_caller"] = cluster_caller
                    test_df["cluster_distance"] = cluster_distance
                    test_df["cluster_p_val"] = cluster_data_dir
                    test_df["sample_set"] = sample_set
                    test_df["cluster_association"] = cluster_association
                    test_df["association_feature"] = association_feature
                    test_df["snp_support"] = snp_support
                    meta_test_df = pd.concat([meta_test_df, test_df], ignore_index=True)

    meta_master_df.to_csv(f"outputs/cluster_length_meta_df_{cluster_caller}.csv", encoding='utf-16', sep='\t', index=False)
    meta_test_df.to_csv(f"outputs/cluster_length_meta_test_df_{cluster_caller}.csv", encoding='utf-16', sep='\t', index=False)