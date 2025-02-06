# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
This file contains functions to combine and aggregate cluster data (P-MACD) from multiple files,
which are specified in the FILES variable in the __main__ section at the bottom of the file.
'''
import pandas as pd
import numpy as np
import logging
from typing import List, Dict, Literal

def aggregate_anz5(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates the input DataFrame by grouping rows by 'Dataset_Cluster_ID' and applying various aggregation functions to the other columns.
    
    Args:
        df (pd.DataFrame): The input DataFrame to be aggregated.
        
    Returns:
        pd.DataFrame: The aggregated DataFrame with one row per unique 'Dataset_Cluster_ID'.
    """

    #create an "End_position" column and set it to the same value as "Start_position"; this is needed for the aggregation. Do this in place.
    df.loc[:, "End_position"] = df["Start_position"].copy()


    df = df[[
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_position",
        "End_position",
        #"Strain_Mutation_ID",
        #"Dataset_Mutation_ID",
        "StrainCluster_ID",
        "Dataset_Cluster_ID",
        "Distance_Between_Clusters",
        "Cluster_Size_Mutations",
        "Cluster_Length",
        "Cluster_Pvalue"]]
   
    agg_dic = {
        "Tumor_Sample_Barcode": "max",
        "Chromosome": "max",
        "Start_position": "min",
        "End_position": "max",
        #"Strain_Mutation_ID": "max",
        #"Dataset_Mutation_ID": "max",
        "StrainCluster_ID": ['max'],
        "Distance_Between_Clusters": ['max'],
        "Cluster_Size_Mutations": ['max'],
        "Cluster_Length": ['max'],
        "Cluster_Pvalue": ['max']}
   
    result = df.groupby('Dataset_Cluster_ID').agg(agg_dic)
    result.columns = result.columns.droplevel(1)
    result = result.reset_index()

    return result

def load_anz5_data(file_names: list, genotype_dict_master: dict, files_loc: str = "data/cluster_files", save: bool = False,
                   save_name: str = "outputs/clusters/all_clusters_data_pooled.csv") -> pd.DataFrame:
    """
    Aggregates the input DataFrame by grouping rows by 'Dataset_Cluster_ID' and applying various aggregation functions to the other columns.
    
    Args:
        df (pd.DataFrame): The input DataFrame to be aggregated.
        
    Returns:
        pd.DataFrame: The aggregated DataFrame with one row per unique 'Dataset_Cluster_ID'.
    """
    #create dataframe with all the data
    df_all = pd.DataFrame()

    for file_name in file_names:
        df = pd.read_csv(f'{files_loc}/{file_name}', sep="\t")

        if file_name == '1_master_para_bed_sorted_anz5.txt':
            genotype_dict = genotype_dict_master[file_name]

        elif file_name == '2018_all_snps_for_clusters_PMACD_sorted_anz5.txt':
            genotype_dict = genotype_dict_master[file_name]
        
        elif file_name == 'MeiosisOutcomesNGS062019_SNVs_sorted_anz5.txt':
            genotype_dict = genotype_dict_master[file_name]

        elif file_name == 'NGS102020Clusters_SNP_table_sorted_anz5.txt':
            #those are spo13 spo13spo11...
            continue

        elif file_name == 'NGS122019_all_SNPs_for_clusters_sorted_anz5.txt':
            #those are exo1 and pol32, by Juraj, this experiment will be ommited
            continue

        elif file_name == 'oct_meiosis_SNVs_sorted_anz5.txt':
            genotype_dict = genotype_dict_master[file_name]
            #strip the "_" from the sample names
            df["Tumor_Sample_Barcode"] = df["Tumor_Sample_Barcode"].str.replace("_", "")

        elif file_name == 'spo13spo11_new_sorted_anz5.txt':
            #first set of spo13 spo13spo11 mutants, we are unsure about the selection, omit
            continue

        elif file_name == 'tetrads_sorted_anz5.txt':
            #contains tetrad data, we are not interested in this here
            continue

        elif file_name == 'novogene72_exo1-ndpol32_sorted_anz5.txt':
            genotype_dict = genotype_dict_master[file_name]

        elif file_name == 'novogeneexo1pol32double_sorted_anz5.txt':
            genotype_dict = genotype_dict_master[file_name]

        else:
            continue

        # aggregate the data by aggregate_anz5
        df = aggregate_anz5(df)

        # add genotype column, if the genotype is not in the dictionary, set it to "not found"
        df["genotype"] = df["Tumor_Sample_Barcode"].map(genotype_dict).fillna("not found")

        # get the keys of the genotype_dict
        genotype_dict_keys = list(genotype_dict.keys())
        
        # print the genotype_dict_keys that are not in the df["Tumor_Sample_Barcode"]
        out = set(genotype_dict_keys) - set(df["Tumor_Sample_Barcode"])
        logging.warning(f"Genotype keys that are not in the df (no clusters found): {out}")

        # if there are keys that are not in the df["Tumor_Sample_Barcode"], add them to the df with value np.NaN
        if len(out) > 0:
            for key in out:
                append_dict = {
                    "Tumor_Sample_Barcode": key,
                    "genotype": genotype_dict[key],
                    "Dataset_Cluster_ID": np.NaN,
                    "Chromosome": np.NaN,
                    "Start_position": np.NaN,
                    "End_position": np.NaN,
                    #"Strain_Mutation_ID": np.NaN,
                    #"Dataset_Mutation_ID": np.NaN,
                    "StrainCluster_ID": np.NaN,
                    "Distance_Between_Clusters": np.NaN,
                    "Cluster_Size_Mutations": np.NaN,
                    "Cluster_Length": np.NaN,
                    "Cluster_Pvalue": np.NaN}

                #concat the append_dict to the df
                df = pd.concat([df, pd.DataFrame(append_dict, index=[0])], ignore_index=True)

        # add the file name to the dataframe
        df["dataset_file_name"] = file_name

        # concat the df to the df_all
        df_all = pd.concat([df_all, df], ignore_index=True)     

    if save:
        df_all.to_csv(save_name, index=False)

    return(df_all)

def load_Pcorr_data(file_names: list, files_loc: str = "data/cluster_files", save: bool = False,
                    save_name: str = "outputs/clusters/Pcorr_data.txt") -> pd.DataFrame:
    """
    Loads and concatenates multiple Pcorr cluster files into a single DataFrame.
    
    Args:
        file_names (list): A list of file names (strings) to be loaded and concatenated.
        files_loc (str): The directory location of the cluster files. Defaults to "cluster_files".
        save (bool): Whether or not to save the concatenated DataFrame to a TXT file. Defaults to False.
        save_name (str): The file name (string) to use when saving the concatenated DataFrame to a TXT file. Defaults to "Pcorr_data.txt".
        
    Returns:
        pd.DataFrame: The concatenated DataFrame containing all Pcorr cluster data.
    """
    df_all = pd.DataFrame()

    for file_name in file_names:
        file_name = f'{file_name.replace("anz5.txt", "")}sum_all_fisher_Pcorr.txt'
        df = pd.read_csv(f'{files_loc}\\{file_name}', sep="\t")

        if file_name == 'oct_meiosis_SNVs_sorted_sum_all_fisher_Pcorr.txt':
            df["Sample"] = df["Sample"].str.replace("_", "")

        df["dataset_file_name"] = file_name
        df_all = pd.concat([df_all, df], ignore_index=True)

    if save:
        df_all.to_csv(save_name, sep="\t", index=False)

    return(df_all)

def cross_check(df: pd.DataFrame, df_check: pd.DataFrame) -> bool:
    """
    Compares the 'clusters' column in the input DataFrame 'df' to the 'clusters' column in the input DataFrame 'df_check'.
    If the two columns are not identical, prints a message indicating which samples have different clusters.
    
    Args:
        df (pd.DataFrame): The input DataFrame to be cross-checked.
        df_check (pd.DataFrame): The DataFrame to compare 'df' to.
        
    Returns:
        bool: True if all clusters are the same, False otherwise.
    """

    #rename clusters to clusters_anz5 in df
    df = df.rename(columns={"clusters": "clusters_anz5"})

    df_merge = df.merge(df_check, on="Sample", how="left")
    df_merge["clusters_diff"] = df_merge["clusters"] - df_merge["clusters_anz5"]
    df_merge = df_merge[df_merge["clusters_diff"] != 0]
    
    if len(df_merge) == 0:
        logging.info("Cross-validation OK! All clusters are the same.")
        return(True)
    else:
        logging.warning("Clusters are not the same! Please check your data!")
        logging.info(df_merge[["Sample", "mutations", "clusters", "clusters_anz5", "clusters_diff"]])
        return(False)

def remove_not_assocaited_clusters(
        df: pd.DataFrame, 
        df_association: pd.DataFrame, 
        inverse: bool = False, 
        association_feature: Literal["GC", "CO", "both"] = "both"
        ) -> pd.DataFrame:
    """remove clusters that are not associated with a switch event"""

    #merge the two dataframes
    df = df.merge(df_association, left_on="Dataset_Cluster_ID", right_on="Cluster_ID_", how="left")

    #remove the rows that are not associated with a switch event
    feature = association_feature
    if association_feature == "both":
        feature = "Clustered"

    if inverse:
        df = df[df[f"{feature}_Proximal_"] == False]
    else:
        df = df[df[f"{feature}_Proximal_"] == True]

    return df

def fill_in_missing_samples(df: pd.DataFrame, genotype_dict: dict) -> pd.DataFrame:
    
    #we need to be able include samples that did not have any mutations/clusters
    #get the keys of the genotype_dict and fill in the missing samples with 0s
    genotype_dict_keys = list(genotype_dict.keys())
    df_keys = list(df["Sample"].unique())
    out = set(genotype_dict_keys) - set(df_keys)

    if len(out) > 0:
        logging.warning(f"Genotype keys that are not in the df (no clusters found): {out}")
        logging.warning(f"Adding {len(out)} samples with 0 clusters, 0 mutations")
        for key in out:
            append_dict = {
                "Sample": key,
                "genotype": genotype_dict[key],
                "clusters": 0,
                "Avg_Cluster_Size_Mutations": 0,
                "Avg_Cluster_Length": 0,
                "mutations": 0}

            #concat the append_dict to the df
            df = pd.concat([df, pd.DataFrame(append_dict, index=[0])], ignore_index=True)
    
    return df

def create_combined_cluster_table(
        cluster_caller: str,
        snp_num: int,
        sig: float | str,
        files: List[str],
        genotype_dict_master: Dict[str, str],
        association: bool|None = None,
        association_feature: Literal["GC", "CO", "both"] = "both",
        min_cluster_mutations: int = 2,
        ) -> None:
    """
    Creates a combined cluster table from ANZ5 and Pcorr data.

    Args:
        files: A list of file paths containing ANZ5 and Pcorr data. Needed when using PMACD cluster caller.
        genotype_dict_master: A dictionary containing genotype information.
        output_name: The name of the output file.
        files_loc: The location where intermediate files will be saved (default: "data/cluster_files").
        association: If True, only include clusters that are associated with a mutation (default: None).

    Returns:
        None

    Raises:
        AssertionError: If the number of clusters in ANZ5 and Pcorr data do not match.
    """
    from mayo.settings import config as cfg

    cluster_association_file = f'outputs/associations/{cluster_caller}/{snp_num}snp/master_clusters_association_{sig}_{snp_num}snp_10000.csv' #_7500

    if cluster_caller == "PMACD":

        files_loc = f"data/cluster_files/clusters_diff_thresholds/master_{sig}"
        df = load_anz5_data(files, genotype_dict_master, save=True, files_loc=files_loc, save_name=f"outputs/clusters/{cluster_caller}/all_{cluster_caller}_cluster_data_pooled_master_{sig}.csv")

        if association is not None:
            df = remove_not_assocaited_clusters(df, df_association=pd.read_csv(cluster_association_file), inverse=(not association), association_feature=association_feature)

        df = df[df["Cluster_Size_Mutations"] >= min_cluster_mutations]

        agg_dict = {"Cluster_Size_Mutations": "mean",
                    "Cluster_Length": "mean",
                    "Dataset_Cluster_ID": "count"} #but if value is np.NaN, then count it as 0

        df = df.groupby(["Tumor_Sample_Barcode", "genotype", "dataset_file_name"]).agg(agg_dict).reset_index()

        df = df.rename(columns={"Tumor_Sample_Barcode": "Sample",
                                "Dataset_Cluster_ID": "clusters",
                                "Cluster_Size_Mutations": "Avg_Cluster_Size_Mutations",
                                "Cluster_Length": "Avg_Cluster_Length"})

        df_pcorr = load_Pcorr_data(files, save=True)
        df_pcorr = df_pcorr[["Sample", "mutations", "clusters"]]

        #check if the number of clusters is the same between anz5 and Pcorr
        #if yes, then we can merge the dataframes and add the mutations columns
        if not cross_check(df, df_pcorr):
            logging.warning("The number of clusters in ANZ5 and Pcorr data do not match!")
            #raise AssertionError("The number of clusters in ANZ5 and Pcorr data do not match!")

        #merge the dataframes
        df = df.merge(df_pcorr[["Sample", "mutations"]], on="Sample", how="left").fillna(0)
    
    elif cluster_caller == "JT":

        files_loc = f"data/new_clusters/mutation_clusters/mutation_clusters_pval_{sig}.tsv"
        scattered_mutations_file = f"data/new_clusters/scattered_mutations/scattered_SNVs_JT_10000_{sig}.tsv"

        df = pd.read_csv(files_loc, sep="\t")
        df = df[df["Mutations"] >= min_cluster_mutations]

        #need to add the genotype column
        genotype_dict = genotype_dict_master[cfg["genotype_dict"]]
        df["genotype"] = df["sample"].map(genotype_dict).fillna("not found")

        # CRITICAL! Still need make sure that cluster_ID is called the same... (critical)
        if association is not None:
            df = remove_not_assocaited_clusters(df, df_association=pd.read_csv(cluster_association_file), inverse=(not association), association_feature=association_feature)

        agg_dict = {
            "Mutations": "mean",
            "Length": "mean",
            "Dataset_Cluster_ID": "count",
            "all_mutations": "max"} #but if value is np.NaN, then count it as 0

        df = df.groupby(["sample", "genotype"]).agg(agg_dict).reset_index()

        df = df.rename(columns={
            "sample": "Sample",
            "Dataset_Cluster_ID": "clusters",
            "Mutations": "Avg_Cluster_Size_Mutations",
            "all_mutations": "mutations",
            "Length": "Avg_Cluster_Length"})
        
        df = fill_in_missing_samples(df, genotype_dict)

        #folowing 4 lines, fix issue where samples with 0 clusters would not have any mutations in output
        scattered_mutations_df = pd.read_csv(scattered_mutations_file, sep='\t', encoding='utf-16')
        scattered_mutations_df = scattered_mutations_df.groupby('Sample').size().reset_index(name='scattered_mutations')
        df = pd.merge(df, scattered_mutations_df, on='Sample', how='left').fillna(0)
        df['mutations'] = df.apply(lambda x: x['scattered_mutations'] if x['scattered_mutations'] > x['mutations'] else x['mutations'], axis=1)
    
    else:
        raise ValueError("cluster_caller must be one of 'JT' or 'PMACD'")

    #add total ssDNA column
    df['total_ssDNA_kb'] = (df['Avg_Cluster_Length'] * df['clusters']) / 1000

    #save the dataframe
    output_name = f"outputs/clusters/{cluster_caller}/{snp_num}snp/all_{cluster_caller}_cluster_data_pooled_master_{snp_num}snp_{sig}_{association_feature}.csv"

    if association == None:
        df.to_csv(f"{output_name[:-4]}_summary.csv", index=False)
    elif association == True:
        df.to_csv(f"{output_name[:-4]}_summary_associated.csv", index=False)
    elif association == False:
        df.to_csv(f"{output_name[:-4]}_summary_not_associated.csv", index=False)

##########################################################################################
if __name__ == "__main__":

    OUTPUT_NAME = "all_cluster_data_pooled_summary_test.csv"
    # FILES = [
    #         '2018_all_snps_for_clusters_PMACD_sorted_anz5.txt',
    #         'MeiosisOutcomesNGS062019_SNVs_sorted_anz5.txt',
    #         'NGS102020Clusters_SNP_table_sorted_anz5.txt',
    #         'NGS122019_all_SNPs_for_clusters_sorted_anz5.txt',
    #         'oct_meiosis_SNVs_sorted_anz5.txt',
    #         'spo13spo11_new_sorted_anz5.txt',
    #         'tetrads_sorted_anz5.txt',
    #         'novogene72_exo1-ndpol32_sorted_anz5.txt',
    #         'novogeneexo1pol32double_sorted_anz5.txt']
    
    FILES = ['1_master_para_bed_sorted_anz5.txt']

    from settings import genotype_dict as genotype_dict_master

    create_combined_cluster_table(FILES, genotype_dict_master, OUTPUT_NAME)