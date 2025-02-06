# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

sig = 'p001_def' # p05, p01, p005, p001_def, p0005, p0001 

config: dict = {
    # "parents_path": r"C:\Users\twaro\development\MayoSwitcher\TheMayoSwitcher\data\sequencing_calls\CombinedSet\Parents\haploid_final",
    # "parental_SNPs_path": r"C:\Users\twaro\My Drive\Meiosis Project\CombinedSet\Switches",
    # "SNVs_path": r"C:\Users\twaro\My Drive\Meiosis Project\CombinedSet\SNVs",
    # "reference_fasta_path": r"C:\Users\twaro\development\MayoSwitcher\TheMayoSwitcher\S288C_reference_sequencewoUra3wUra329wNAT.fa",
    "parents_path": "data/sequencing_calls/CombinedSet/Parents/haploid_final",
    "parental_SNPs_path": "data/sequencing_calls/CombinedSet/Switches",
    "SNVs_path": "data/sequencing_calls/CombinedSet/SNVs",
    "reference_fasta_path": "S288C_reference_sequencewoUra3wUra329wNAT.fa",
    "repetitive_regions_fasta_path": "data/other_features_genomic.fasta",
    "cluster_table_path": rf"C:\Users\twaro\development\MayoSwitcher\TheMayoSwitcher\data\cluster_files\clusters_diff_thresholds\master_{sig}\1_master_para_bed_sorted_anz5.txt",
    "diploid": True,
    "mitochondrial_chr_name": ["ref|NC_001224|"],
    "SNP_NUM": 1,
    "genotype_dict": "final_set", #or "1_master_para_bed_sorted_anz5.txt" (full)
    
    "CLUSTER_TYPE_ANALYSIS": "JT", # JT or PMACD,
    "SIG": 0.0001, # 0.01, 0.001, 0.0001, 0.00001

    #CLUSTER_TYPE_ANALYSIS = "PMACD"
    #SIG = 'p001_def' # p05, p01, p005, p001_def, p0005, p0001

    "CLUST_INTER_MUT_MAX": 10000,
    "mitochondrial_chr_name": ["ref|NC_001224|"],
    }

config_cluster: dict = {
    "PATH": "data/sequencing_calls/CombinedSet/SNVs",
    "PATH_CLUSTER_RESULTS": f"data/cluster_files/clusters_diff_thresholds/master_{sig}",
    "FILTER_INPUT": True,
    "MAX_DUPLICATES":  15,
    "SAVE_DUPLICATES_FILE": "outputs/clusters/1_all_duplicates.csv",
    "SAVE_MASTER_FILE": "outputs/clusters/1_master_para_bed.txt",
    "CLUSTER_DATA_FILE": f"outputs/clusters/all_cluster_data_pooled_master_{sig}.csv"
}