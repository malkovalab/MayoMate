# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
from typing import Union, Literal, overload, TypeVar, List
import BioAid as ba

from .ParentalTable import ParentalTable
from .SampleTable import SampleTable
from .OverlapSwitchesTable import OverlapSwitchesTable
from .SwitchesTable import SwitchesTable
from .decorators import fancy_status, time_elapsed
from .TableSNP import TableSNP
from .base_functions import get_samples_with_genotype, create_similarity_matrix

from mayo import *
from mayo.settings import chromosomes as chromosomes
from mayo.settings import config as cfg
from mayo.settings import roman_to_S288C


logger = logging.getLogger('__main__.' + __name__)

def readjust_switches_to_spo11_signal(
        switches_list: list[SwitchesTable],
        spo11_map: pd.DataFrame,
        adjust_window: int = 1000, 
        mode: Literal["strongest_peak", "nearest_peak"] = "strongest_peak"
        ) -> list[SwitchesTable]:
    """
    Readjusts the switch centers of each sample in the switches_list to the SPO11 signal.

    This function iterates over each sample in the switches_list. For each sample, it readjusts the switch centers 
    to the SPO11 signal using the provided spo11_map and adjust_window. The mode of adjustment can be specified 
    (default is "strongest_peak"). After readjustment, it checks for any duplicated rows in the sample DataFrame 
    and raises an assertion error if any are found.

    Parameters:
    switches_list (list): A list of Sample objects. Each Sample object should have a DataFrame (df) with a "Switch_Center" column.
    spo11_map (pandas.DataFrame): A DataFrame representing the SPO11 signal map.
    adjust_window (int, optional): The window size for adjustment. Default is 1000.
    mode (str, optional): The mode of adjustment. Default is "strongest_peak".

    Returns:
    list: The list of Sample objects with readjusted switch centers.
    """    
    
    for sample in switches_list:
        print(f"Processing sample {sample.name}")
        sample.df["Original_Switch"] = sample.df["Switch_Center"]
        sample.readjustSwitchesToSPO11Signal(
            spo11_oligo_map_df=spo11_map,
            adjust_window=adjust_window,
            mode=mode
            )

        #find if there are duplicated rows in sample.df, there shouldn't be any
        #this is to ensure that the readjustment process is correct and no other errors occurr later
        duplicated_rows = sample.df[sample.df.duplicated(subset=["Chromosome", "Switch_Center"], keep=False)]
        assert duplicated_rows.empty, f"Found duplicated rows in sample {sample.name} after readjusting switches to SPO11 signal: {duplicated_rows}"

        sample.df["Switch_Center"] = sample.df["Adjusted_Center"]

    return switches_list

def readjust_GC_CO_to_spo11_signal(ev_dict, spo11_map_filtered, window=2000) -> dict:
    
    GC_dict_out = {}
    
    for chr in ev_dict:
        new_chr_positions = []
        positions_on_chr_list = ev_dict[chr]

        for old_position in positions_on_chr_list:

            df_peaks = spo11_map_filtered[(spo11_map_filtered["Chromosome"] == chr) & (spo11_map_filtered["Position"] < old_position + window) & (spo11_map_filtered["Position"] > old_position - window)].copy()

            if len(df_peaks) > 0:

                peak = df_peaks.loc[df_peaks["Value"].idxmax()]
                new_chr_positions.append(peak['Position'])
            else:
                new_chr_positions.append(old_position)

        #check if there are any duplicates in the new_chr_positions list and remove them
        before = len(new_chr_positions)
        new_chr_positions = list(dict.fromkeys(new_chr_positions))
        after = len(new_chr_positions)
        if before != after:
            print(f"Removed {before - after} duplicates from {chr}")

        GC_dict_out[chr] = new_chr_positions

    return GC_dict_out

def find_possible_sample_duplicacates(samples_parental_SNPs, sample_names, save_path: str = "figures/similarity_matrix_random_spores_nodups.png"):

    samples_subset = [i for i in samples_parental_SNPs if i.name in sample_names]
    print(f"Number of samples in subset: {len(samples_subset)}")

    similarity_matrix = create_similarity_matrix(samples_subset)

    similarity_matrix = np.round(similarity_matrix, 2)

    plt.figure(figsize=(80, 80))    
    mask = np.zeros_like(similarity_matrix)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(similarity_matrix, annot=True, xticklabels=[sample.name for sample in samples_subset], yticklabels=[sample.name for sample in samples_subset], mask=mask, vmax=1, vmin=0, cmap="coolwarm")
    #make y axis labels horizontal
    plt.yticks(rotation=0)

    #save the heatmap
    plt.savefig(save_path, dpi=200)
    plt.show()

def prep_df_master_for_syn_vs_non_syn_analysis(df_master: pd.DataFrame) -> pd.DataFrame:
    """
    Prepares the master DataFrame for synonymous vs non-synonymous analysis by grouping by 
    Chromosome and Region, and counting the number of mutations in each group.

    Parameters:
    df_master (pd.DataFrame): The input DataFrame containing genomic data.

    Returns:
    pd.DataFrame: The processed DataFrame ready for analysis.
    """

    df_analysis = df_master.copy()
    df_analysis = df_analysis[df_analysis["Chromosome"] != "ref|NC_001224|"]
    print(f"There are {len(df_analysis)} genomic mutations in df_analysis.")
    df_analysis = df_analysis.groupby(['Chromosome', 'Region']).size().reset_index(name='Count').sort_values(by=['Chromosome', 'Region'])
    df_analysis["Position"] = df_analysis["Region"]
    return df_analysis

def calculateSynVsNonSynMutations(
        df_analysis, 
        path,
        sample_counter=179,
        just_high_impact=False) -> None:
    """
    Calculate the rate of synonymous vs non-synonymous mutations.

    Parameters:
    df_analysis (pd.DataFrame): DataFrame containing the analysis data.
    path (str): Path to the CSV file containing cleaned data.
    sample_counter (int, optional): Number of samples analyzed. Default is 179.
    just_high_impact (bool, optional): If True, only consider high impact mutations. Default is False.

    Returns:
    None
    """
    
    from mayo.settings import roman_to_S288C

    cleaned_df = pd.read_csv(path, sep='\t')
    
    if just_high_impact:
        cleaned_df = cleaned_df[cleaned_df["IMPACT"] == "HIGH"]

    cleaned_df['Chromosome'] = cleaned_df['Chromosome'].replace(roman_to_S288C)
    cleaned_df = addGenomeReferenceSystemOffesets(cleaned_df)

    #count number of unique Chromosome/Position pairs in cleaned_df
    cleaned_df = cleaned_df[['Chromosome', 'Position']].drop_duplicates()
    unique_sites_count = len(cleaned_df)

    df_master_analysis_cleaned = df_analysis.merge(cleaned_df, on=['Position', 'Chromosome'], how='inner')

    #also show entries in df_analysis that are not contained in cleaned_df by merging on Position and Chromosome
    not_df_master_analysis_cleaned = df_analysis.merge(cleaned_df, on=['Position', 'Chromosome'], how='outer', indicator=True)
    not_df_master_analysis_cleaned = not_df_master_analysis_cleaned[not_df_master_analysis_cleaned['_merge'] == 'left_only']

    # if path == "cleaned_data2.tsv":
    #     not_df_master_analysis_cleaned.to_csv(f"not_in_cleaned_data.tsv", sep='\t', index=False)

    print(not_df_master_analysis_cleaned["Count"].sum())
    
    mutations_in_sites = df_master_analysis_cleaned["Count"].sum()

    #calculate rate
    rate = mutations_in_sites / unique_sites_count / sample_counter
    print(f"Analyzed {sample_counter} samples with {mutations_in_sites} mutations in {unique_sites_count} unique sites in {path}.")
    print(f"  The rate of mutations in selected sites is {rate:.4e} mutations per unique site per sample.")

@time_elapsed
def pull_SNP_vs_mutation_vs_recombination_data(df_mutations, df_SNPs, df_all_rec_events, window_size=200, step_size=None):

    if step_size is None:
        step_size = window_size

    df_SNPs.reset_index(inplace=True)
    df_SNPs = df_SNPs[df_SNPs["Parent"]=="parent_1"]
    df_SNPs = df_SNPs[df_SNPs["Chromosome"] != "ref|NC_001224|"].reset_index(drop=True)
    df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])

    df_all_rec_events = df_all_rec_events[df_all_rec_events["Chromosome"] != "ref|NC_001224|"].reset_index(drop=True)

    df_mutations = df_mutations.copy()
    df_mutations = removeTelomereSNPs(df_mutations)
    df_mutations = removeRepetitiveDNASNPs(df_mutations, repetitive_DNA_fasta_path=cfg["repetitive_regions_fasta_path"])

    df_SNPs_mutations_window = pd.DataFrame(columns=["Chromosome", "Position_start", "Position_end", "SNP_count", "Mutation_count"])

    for chr in df_SNPs["Chromosome"].unique():

        chr_len = df_chr[df_chr["chromosome"] == chr]["end_position"].values[0]

        df_SNPs_subset = df_SNPs[df_SNPs["Chromosome"] == chr]
        df_mutations_subset = df_mutations[df_mutations["Chromosome"] == chr]
        df_recev_subset = df_all_rec_events[df_all_rec_events["Chromosome"] == chr]

        for i in range(0, chr_len, step_size):
    
            if i + window_size > chr_len: #make sure the window does not go over the end of the chromosome
                continue
            
            mutations_in_window = df_mutations_subset[(df_mutations_subset["Region"] >= i) & (df_mutations_subset["Region"] < i + window_size)]
            SNPs_in_window = df_SNPs_subset[(df_SNPs_subset["Region"] >= i) & (df_SNPs_subset["Region"] < i + window_size)]
            recombination_in_window = df_recev_subset[(df_recev_subset["Region"] >= i) & (df_recev_subset["Region"] < i + window_size)]

            window_dict = {
                "Chromosome": chr,
                "Position_start": i,
                "Position_end": i + window_size,
                "SNP_count": len(SNPs_in_window),
                "Mutation_count": len(mutations_in_window),
                "Rec_event_count": len(recombination_in_window)
            }
            #concatenate the window_dict to the df_SNPs_mutations_window dataframe
            df_SNPs_mutations_window = pd.concat([df_SNPs_mutations_window, pd.DataFrame([window_dict])])

    return df_SNPs_mutations_window

def pull_SNP_vs_recombination_data(df_rec_events, df_SNPs, window_size=200, step_size=None):

    if step_size is None:
        step_size = window_size

    df_SNPs.reset_index(inplace=True)
    df_SNPs = df_SNPs[df_SNPs["Parent"]=="parent_1"]
    df_SNPs = df_SNPs[df_SNPs["Chromosome"] != "ref|NC_001224|"].reset_index(drop=True)

    df_rec_events = df_rec_events[df_rec_events["Chromosome"] != "ref|NC_001224|"].reset_index(drop=True)

    df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])

    df_SNPs_recev_window = pd.DataFrame(columns=["Chromosome", "Position_start", "Position_end", "SNP_count"])

    for chr in df_SNPs["Chromosome"].unique():

        chr_len = df_chr[df_chr["chromosome"] == chr]["end_position"].values[0]

        df_SNPs_subset = df_SNPs[df_SNPs["Chromosome"] == chr]
        df_recev_subset = df_rec_events[df_rec_events["Chromosome"] == chr]

        for i in range(0, chr_len, step_size):
    
            if i + window_size > chr_len: #make sure the window does not go over the end of the chromosome
                continue
            
            recombination_in_window = df_recev_subset[(df_recev_subset["Region"] >= i) & (df_recev_subset["Region"] < i + window_size)]
            SNPs_in_window = df_SNPs_subset[(df_SNPs_subset["Region"] >= i) & (df_SNPs_subset["Region"] < i + window_size)]

            window_dict = {
                "Chromosome": chr,
                "Region": i,
                "Position_start": i,
                "Position_end": i + window_size,
                "SNP_count": len(SNPs_in_window),
                "RecEv_count": len(recombination_in_window)
            }
            #concatenate the window_dict to the df_SNPs_mutations_window dataframe
            df_SNPs_recev_window = pd.concat([df_SNPs_recev_window, pd.DataFrame([window_dict])])
    
    
    df_SNPs_recev_window["Region"] = df_SNPs_recev_window["Region"].astype(int)
    df_SNPs_recev_window["SNP_count"] = df_SNPs_recev_window["SNP_count"].astype(int)
    df_SNPs_recev_window["RecEv_count"] = df_SNPs_recev_window["RecEv_count"].astype(int)
    df_SNPs_recev_window["Position_start"] = df_SNPs_recev_window["Position_start"].astype(int)
    df_SNPs_recev_window["Position_end"] = df_SNPs_recev_window["Position_end"].astype(int)
    
    df_SNPs_recev_window = removeTelomereSNPs(df_SNPs_recev_window)
    df_SNPs_recev_window = removeRepetitiveDNASNPs(df_SNPs_recev_window, repetitive_DNA_fasta_path=cfg["repetitive_regions_fasta_path"])
    
    df_SNPs_recev_window.reset_index(drop=True, inplace=True)

    return df_SNPs_recev_window

def findTranscriptAndDirection(
        df: pd.DataFrame, 
        features_path: str,
        include_UTRs: bool = False,
        UTR_data_path: str = None
		) -> pd.DataFrame:

    #sort the df by Chromosome and Region columns
    df = df.sort_values(by=["Chromosome", "Region"])

    df["Transcript"] = ""
    df["Direction"] = ""
    df["Feature_length"] = 0

    df_features_summary = pd.DataFrame()

    #open the features fasta file with BioAid
    DNA_list = ba.extractSeqFromFastaToList(features_path)
    if include_UTRs:
        df_UTR, df_UTR_3, df_UTR_5 = prep_yeast_mine_data(path=UTR_data_path)

    analyzed_features = 0

    for seq_name, seq in DNA_list:

        #extract name from seq_name
        feature_name = seq_name.split(" ")[0].strip(">").strip()
        chromosome_name = seq_name.split(",")[1].split(" ")[2].strip()            

        # if feature_name not in non_essential_genes:
        #     continue

        # if "ARS" not in seq_name: #tRNA_gene
        #     continue

        if chromosome_name not in roman_to_S288C.keys():
            logging.warning(f"Chromosome {chromosome_name} of {feature_name} is not in roman_to_S288C.keys(). This behavior is normal for Mitochondrial Sequences. Skipping...")
            continue
        
        analyzed_features += 1
        print(f"Evaluating {seq_name}")

        tct_count, tca_count, tcg_count, aga_count, tga_count, cga_count = getA3A_motifs((seq_name, seq))

        chromosome_name = roman_to_S288C[chromosome_name]
        seq_coords = seq_name.split(",")[1].split(" ")[-1].strip()
        seq_start = int(seq_coords.split("-")[0])
        seq_end = int(seq_coords.split("-")[1])

        if seq_start > seq_end:
            feature_direction = "-"
        else:
            feature_direction = "+"
            
        if include_UTRs:

            #extract the 5' UTR start and end positions
            if feature_name in df_UTR_5['gene_name'].values:

                utr_seq_start = df_UTR_5[df_UTR_5['gene_name'] == feature_name]['start'].values[0]
                utr_seq_end = df_UTR_5[df_UTR_5['gene_name'] == feature_name]['end'].values[0]
                if feature_direction == "+":
                    seq_start = utr_seq_start
                if feature_direction == "-":
                    seq_start = utr_seq_end

            # #extract the 3' UTR start and end positions
            # if feature_name in df_UTR_3['gene_name'].values:
    
            #     utr_seq_start = df_UTR_3[df_UTR_3['gene_name'] == feature_name]['start'].values[0]
            #     utr_seq_end = df_UTR_3[df_UTR_3['gene_name'] == feature_name]['end'].values[0]
            #     if feature_direction == "+":
            #         seq_end = utr_seq_end
            #     if feature_direction == "-":
            #         seq_end = utr_seq_start
        
            if feature_direction == "+":
                seq_end = seq_start + 100
                seq_start = seq_start - 100
            if feature_direction == "-":
                seq_end = seq_start - 100
                seq_start = seq_start + 100

        reporter_start_pos_chr3 = 211813
        offset_reporter = 2723
        if (chromosome_name == "ref|NC_001135|") & (seq_start >= reporter_start_pos_chr3):
            seq_start += offset_reporter
            seq_end += offset_reporter
            logging.info(f"Added offset of ura3-29 reporter ({offset_reporter}) to {feature_name} on chromosome {chromosome_name} (region: {seq_start}-{seq_end})")

        #ura3_start_pos_chr5 = 116166
        ura3_end_pos_chr5 = 116970
        ura3_del_len = 805
        if (chromosome_name == "ref|NC_001137|") & (seq_start > ura3_end_pos_chr5):
            #remove the length of the ura3 deletion from the start and end positions
            seq_start -= ura3_del_len
            seq_end -= ura3_del_len
            logging.info(f"Removed length of ura3 ({ura3_del_len}) from {feature_name} on chromosome {chromosome_name} (region: {seq_start}-{seq_end})")

        #sort the genomic coordinates, this is needed because the repetitive DNA file has some coordinates in reverse order (transcription direction)
        seq_genomic_coord = [seq_start, seq_end]
        seq_genomic_coord.sort()

        #find the SNVs that are in the feature and add the feature name to the Transcript column and the direction to the Direction column
        df_feature = df[(df["Chromosome"] == chromosome_name) & (df["Region"] >= seq_genomic_coord[0]) & (df["Region"] <= seq_genomic_coord[1])].copy()
        df_feature["Transcript"] = feature_name
        df_feature["Direction"] = feature_direction
        df.update(df_feature)

        feature_count = len(df_feature)
        feature_length = seq_genomic_coord[1] - seq_genomic_coord[0] + 1
        feature_G_count = df_feature[df_feature["Reference"] == "G"]["Reference"].count()
        feature_C_count = df_feature[df_feature["Reference"] == "C"]["Reference"].count()
        feature_A_count = df_feature[df_feature["Reference"] == "A"]["Reference"].count()
        feature_T_count = df_feature[df_feature["Reference"] == "T"]["Reference"].count()

        feature_mut_rate = feature_count / feature_length
        feature_A3A_rate = (feature_G_count + feature_C_count) / feature_length

        df_feature["Feature_Length"] = feature_length

        concat_dict = {
            "feature_name": feature_name,
            "feature_length": feature_length,
            "feature_count": feature_count,
            "feature_mut_rate": feature_mut_rate,
            "feature_A3A_rate": feature_A3A_rate,
            "feature_G_count": feature_G_count,
            "feature_C_count": feature_C_count,
            "feature_A_count": feature_A_count,
            "feature_T_count": feature_T_count,
            "feature_direction": feature_direction,
            "feature_tct_count": tct_count,
            "feature_tca_count": tca_count,
            "feature_tcg_count": tcg_count,
            "feature_aga_count": aga_count,
            "feature_tga_count": tga_count,
            "feature_cga_count": cga_count,
            "feature_A3A_motifs_top": tct_count + tca_count + tcg_count,
            "feature_A3A_motifs_bottom": aga_count + tga_count + cga_count,
            "feature_chromosome": chromosome_name,
            "feature_start": seq_start,
            "feature_end": seq_end,

        }
        df_concat = pd.DataFrame(concat_dict, index=[0])

        df_features_summary = pd.concat([df_features_summary, df_concat], ignore_index=True)

    print(f"Evaluated {analyzed_features} features in {features_path}")

    return df, df_features_summary

#Null hypothesis simulations functions:

def create_summary_null_df(simulation_pval_sample, min_clust_mutations = 3) -> pd.DataFrame:
    from numpy import NaN


    df_cluster_null = pd.DataFrame(columns=["cluster_sim_id", "total_clusters", "total_ssDNA", "average_cluster_length", "average_cluster_mutations"])
    df_cluster_null = df_cluster_null.astype({"cluster_sim_id": int, "total_clusters": int, "total_ssDNA": int, "average_cluster_length": float, "average_cluster_mutations": float})
    for sim_index in simulation_pval_sample:
        #print(f"sim_index: {sim_index}")

        df = simulation_pval_sample[sim_index]

        #make sure that the df has at least n mutations
        df = df[df["Mutations"] >= min_clust_mutations]

        #get the total number of clusters (length of the df)
        total_clusters = len(df)

        #get the total amount of ssDNA by summing up th Length column
        total_ssDNA = df["Length"].sum()

        average_cluster_length = df["Length"].mean()
        if average_cluster_length is NaN:
            average_cluster_length = 0

        average_cluster_mutations = df["Mutations"].mean()
        if average_cluster_mutations is NaN:
            average_cluster_mutations = 0

        #get total clustered mutations by summing up the Mutations column
        total_clustered_mutations = df["Mutations"].sum()

        #add the data to the df_cluster_null
        concat_dict = {
            "cluster_sim_id": sim_index,
            "total_clusters": total_clusters, 
            "total_ssDNA": total_ssDNA, 
            "average_cluster_length": average_cluster_length, 
            "average_cluster_mutations": average_cluster_mutations,
            "total_clustered_mutations": total_clustered_mutations}
        
        df_concat = pd.DataFrame(concat_dict, index=[0])

        df_cluster_null = pd.concat([df_cluster_null, df_concat], ignore_index=True)

    return df_cluster_null
    
def draw_distribution_of_null_total_ssDNA(null_total_ssDNA_list):
    import numpy as np
    plt.figure(figsize=(10, 5))
    sns.distplot(null_total_ssDNA_list, kde=False, bins=np.arange(0, 150000, 5000))
    plt.title("Distribution of average null total ssDNA")
    plt.show()

def draw_distribution_of_null_cluster_mutations(null_cluster_mutations_list):
    import numpy as np
    plt.figure(figsize=(10, 5))
    sns.distplot(null_cluster_mutations_list, kde=False, bins=np.arange(0, 15, 1))
    plt.title("Distribution of average null cluster mutations")
    plt.show()

def draw_distribution_of_null_cluster_length(null_cluster_length_list):
    import numpy as np
    #draw distribution of cluster_length 95th percentile
    plt.figure(figsize=(10, 5))
    sns.distplot(null_cluster_length_list, kde=False, bins=np.arange(0, 50000, 500))
    plt.title("Distribution of average null cluster length")
    plt.show()

def draw_sample_null_summary_plots(df_all_clust_sample):
    
    import numpy as np
    #plot the distributions of the total_ssDNA, Length, and Mutations columns
    plt.figure(figsize=(10, 5))
    sns.distplot(df_all_clust_sample["Length"], kde=False, bins=np.arange(0, 50000, 500))
    plt.title("Distribution of cluster length")
    plt.show()

    plt.figure(figsize=(10, 5))
    sns.distplot(df_all_clust_sample["Mutations"], kde=False, bins=np.arange(0, 15, 1))
    plt.title("Distribution of cluster mutations")
    plt.show()

    #draw a scatter plot of the Length vs Mutations with a regression line
    plt.figure(figsize=(10, 10))
    sns.regplot(x="Mutations", y="Length", data=df_all_clust_sample, scatter=True, order=2, ci=99, scatter_kws={"s": 60, "alpha": 0.5})
    plt.title("Mutations vs Length")
    plt.show()

    #draw a jointplot of the Length vs Mutations
    plt.figure(figsize=(20, 20))
    sns.jointplot(x="Length", y="Mutations", data=df_all_clust_sample, kind="kde")
    plt.show()

    #draw a stripplot of the Length vs Mutations
    plt.figure(figsize=(10, 10))
    sns.stripplot(x="Mutations", y="Length", data=df_all_clust_sample, alpha=0.5)
    plt.show()

    plt.figure(figsize=(10, 10))
    sns.violinplot(x="Mutations", y="Length", data=df_all_clust_sample)
    plt.show()

def getSampleNames(combined_simulations) -> list:
    #extract all sample names (keys) from the combined_simulations dictionary
    sample_names = []
    for i in combined_simulations:
        sample_names.append(i)
    return sample_names

def describe_sample_null_df(df_cluster_null):

    print(df_cluster_null.sort_values(by="total_ssDNA", ascending=True))
    print(df_cluster_null["total_ssDNA"].unique()) 
    print(df_cluster_null.describe())
    print(df_cluster_null.quantile([0.9, 0.95, 0.99]))