# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
This code is used to analyze the PMACD data. It is divided into three parts:
1. Create clean pmacd input (mostly remove duplicates). Can be skipped if the input is already clean.
2. Create combined cluster table. Aggregates the data from multiple files into one table for analysis.
3. Conduct cluster data analysis. Conducts the statistical analysis (specified genotypes) and creates the plots.
'''
import logging
from mayo.settings import config_cluster as cfg_c
from mayo.settings import genotype_dict as genotype_dict
from mayo.clusters import create_clean_pmacd_input
from mayo.clusters import create_combined_cluster_table
from mayo.clusters import conduct_cluster_data_analysis

logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s:.%(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

SNP_NUM = 1
ASSOCIATION = None # None, True or False
ASSOCIATION_FEATURE = "both" # "GC", "CO", "both"
CLUSTER_CALLER = "JT" # JT or PMACD (deprecated)
SIG = 0.0001 #JT; 0.01, 0.001, 0.0001, 0.00001
PLOIDY = 'haploid' #all, haploid, aneuploid
min_cluster_mutations = 3


if __name__ == '__main__':

    ## Part 1: Create clean pmacd input
    FILTER_INPUT: bool = cfg_c["FILTER_INPUT"]
    #remove poor quality data in create_clean_pmacd_input? What is poor quality data?
    # create_clean_pmacd_input(
    #     path = cfg_c["PATH"], 
    #     max_duplicates = cfg_c["MAX_DUPLICATES"], 
    #     save_master_file = cfg_c["SAVE_MASTER_FILE"], 
    #     save_duplicates_file = cfg_c["SAVE_DUPLICATES_FILE"])

    ## Part 2: Create combined cluster table
    #make sure to change below function so that the right cluster association is read in
    create_combined_cluster_table(
        cluster_caller = CLUSTER_CALLER, #JT or PMACD
        snp_num=SNP_NUM,
        sig=SIG,
        files = ['1_master_para_bed_sorted_anz5.txt'],
        genotype_dict_master = genotype_dict, 
        association = ASSOCIATION, # None, True or False
        association_feature = ASSOCIATION_FEATURE,
        min_cluster_mutations = min_cluster_mutations
        )

    ## Part 3: Conduct cluster data analysis
    # genotype_list = [
        # 'ung1∆',
        # "ung1∆NAT",
        # "exo1-nd",
        # "pol32∆",
        # "exo1-ndpol32∆",
        # 'UNG1',
        # "ung1∆ EV",
        # 'ung1∆ non-selected',
        # 'ung1∆ premeiotic non-selected',
        # 'inc. Tetrad ung1∆ non-selected',
        # "spo13∆",
        # "spo13∆spo11∆"
        # ]

    #don't need genotypes_include, just specify the order

    #order, figsize = ['ung1∆', "exo1-nd", "pol32∆", "exo1-ndpol32∆", "sgs1∆C", "exo1-ndsgs1∆C", "ung1∆NAT"], (8, 5) "ung1∆NAT",
    #order, figsize = ["inc. Tetrad ung1∆ non-selected",'ung1∆ premeiotic non-selected', 'UNG1', 'ung1∆ non-selected', 'ung1∆', "ung1∆NAT"], (10, 5) #For fig 1
    #order, figsize = ['spo13∆', 'spo13∆spo11∆'], (4,5) #'ung1∆',
    #order, figsize = ['ung1∆ premeiotic non-selected','UNG1', 'ung1∆ EV', 'ung1∆ non-selected', 'ung1∆', "ung1∆NAT"], (10, 5)# For fig 1
    order, figsize = ['ung1∆ non-selected', 'ung1∆', "ung1∆NAT"], (5, 5) #For fig 2
    #order, figsize = ['UNG1', 'ung1∆ EV', 'ung1∆ non-selected', 'ung1∆', "ung1∆NAT", "exo1-nd", "pol32∆", "exo1-ndpol32∆"], (10, 5) #For fig 2 supplement
    #order, figsize = ['ung1∆ premeiotic non-selected', 'UNG1', 'ung1∆ EV', 'ung1∆ non-selected', 'ung1∆', "ung1∆NAT", "exo1-nd", "pol32∆", "exo1-ndpol32∆", "sgs1∆C", "exo1-ndsgs1∆C", 'spo13∆', 'spo13∆spo11∆'], (20, 5)


    palette = {
        'ung1∆ premeiotic non-selected': 'tab:orange',
        'UNG1': 'tab:olive',
        'ung1∆ EV': 'tab:brown',
        'ung1∆ non-selected': 'green',
        'ung1∆': 'tab:blue', # or 'tab:red' when with spo13∆ and spo13∆spo11∆
        'ung1∆NAT': 'tab:purple',

        'exo1-nd': 'tab:red',
        'pol32∆': 'tab:green',
        'exo1-ndpol32∆': 'tab:purple',
        'sgs1∆C': 'tab:pink',
        'exo1-ndsgs1∆C': 'tab:gray',

        'spo13∆': 'tab:blue', 
        'spo13∆spo11∆': 'tab:orange'
        }

    #if you specify association != None, then we can't use the: ["ung1∆ premeiotic non-selected", "spo13∆", "spo13∆spo11∆"]
    # if ASSOCIATION != None:
    #     #remove the genotypes that don't have the association feature since they are diploid/premeiotic
    #     genotypes.remove("ung1∆ premeiotic non-selected")
    #     genotypes.remove("spo13∆")
    #     genotypes.remove("spo13∆spo11∆")

    conduct_cluster_data_analysis(
        
        cluster_caller = CLUSTER_CALLER, #JT or PMACD
        snp_num = SNP_NUM,
        sig = SIG,
        association = ASSOCIATION, # None, True or False
        association_feature = ASSOCIATION_FEATURE,
        ploidy = "all", #PLOIDY #all, haploid, anueploid
        genotype1 = 'ung1∆', #'spo13∆', #'ung1∆',          #for stats: UNG1 ung1∆
        genotype2 = 'ung1∆NAT', #'ung1∆NAT', #'exo1-ndsgs1∆C',  #for stats: ung1∆NAT', UNG1, exo1-nd, pol32∆, exo1-ndpol32∆, ung1∆ premeiotic non-selected
        save_plt = True, # save the plot as a png file 
        save_plotly = True, # save the plot as an interactive html file, 
        genotypes_include = order,
        figsize = figsize,
        order = order,
        palette = palette, #remove palette param to make all samples tab:blue
        )