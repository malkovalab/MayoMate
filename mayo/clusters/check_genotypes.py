# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
This script finds the number of samples in each genotype by ploidy class and prints the results.
'''
import logging
try:
    from ..settings import genotype_dict as genotype_dict
except ImportError:
    pass

logging.basicConfig(level=logging.INFO, format='%(levelname)s:. %(message)s')
logger = logging.getLogger(__name__)

haploid_list = ['JK2_10B_', 'JK2_11B_', 'JK2_12B_', 'JK2_13_G', 'JK2_14_A', 'JK2_15_C', 'JK2_19B_', 'JK2_1_CT', 'JK2_21_T', 'JK2_22_A', 'JK2_23_A', 'JK2_24_A', 'JK2_25_T', 'JK2_26_G', 'JK2_27_A', 'JK2_28_T', 'JK2_29_G', 'JK2_2_TA', 'JK2_30_A', 'JK2_31_G', 'JK2_32_G', 'JK2_33_A', 'JK2_34_G', 'JK2_35_C', 'JK2_36_A', 'JK2_37_C', 'JK2_38_G', 'JK2_39_G', 'JK2_40_G', 'JK2_41_T', 'JK2_42_C', 'JK2_43_A', 'JK2_44_G', 'JK2_45_T', 'JK2_47_A', 'JK2_49_C', 'JK2_4_TG', 'JK2_51_T', 'JK2_52_G', 'JK2_53_T', 'JK2_56_C', 'JK2_57_A', 'JK2_58_A', 'JK2_59_C', 'JK2_5_GA', 'JK2_60_A', 'JK2_61_C', 'JK2_62_C', 'JK2_63_A', 'JK2_64_G', 'JK2_65_C', 'JK2_68_A', 'JK2_70_G', 'JK2_71_C', 'JK2_72_G', 'JK2_73_G', 'JK2_74_G', 'JK2_75_A', 'JK2_77_C', 'JK2_78_T', 'JK2_7_GT', 'JK2_91_A', 'JK2_98_A', 'JK2_9B_C', 'JMT49_CK', 'JMT50_CK', 'JMT51_CK', 'JMT52_CK', 'JMT53_CK', 'JMT54_CK', 'JMT55_CK', 'JMT56_CK', 'JMT73_CK', 'JMT74_CK', 'JMT76_CK', 'JMT77_CK', 'JMT78_CK', 'JMT78_CK', 'JMT79_CK', 'JMT79_CK', 'JMT80_CK', 'JMT81_CK', 'JMT81_CK', 'JMT82_CK', 'JMT83_CK', 'JMT84_CK', 'JMT85_CK', 'JMT85_CK', 'JMT86_CK', 'JMT87_CK', 'JMT89_CK', 'JMT90_CK', 'JMT90_CK', 'JMT91_CK', 'JMT92_CK', 'JMT94_CK', 'JMTD_10_', 'JMTD_12_', 'JMTD_13_', 'JMTD_15_', 'JMTD_16_', 'JMTD_17_', 'JMTD_18_', 'JMTD_19_', 'JMTD_1_C', 'JMTD_21_', 'JMTD_22_', 'JMTD_23_', 'JMTD_2_C', 'JMTD_3_C', 'JMTD_5_C', 'JMTD_6_C', 'JMTD_8_C', 'JMTD_9_C', 'JT13_CKD', 'JT14_CKD', 'JT15_CKD', 'JT16_CKD', 'JT17_CKD', 'JT21_CKD', 'JT24_CKD', 'JT25_CKD', 'JT27_CKD', 'JT28_CKD', 'JT31_CKD', 'JT36_CKD', 'JT37_CKD', 'JT38_CKD', 'JT39_Mer', 'JT44_Mer', 'JT49_CKD', 'JT50_CKD', 'JT53_CKD', 'JT55_CKD', 'JT57_CKD', 'JT61_CKD', 'JT62_CKD', 'JT63_CKD', 'Mal_JK11', 'Mal_JK14', 'Mal_JK16', 'Mal_JK17', 'Mal_JK19', 'Mal_JK21', 'Mal_JK23', 'Mal_JK24', 'Mal_JK26', 'Mal_JK27', 'Mal_JK29', 'Mal_JK30', 'Mal_JK31', 'Mal_JK35', 'Mal_JK38', 'Mal_JK40', 'Mal_JK41', 'Mal_JK42', 'Mal_JK44', 'Mal_JK45', 'Mal_JK5-', 'Mal_JK58', 'Mal_JK59', 'Mal_JK6-', 'Mal_JK60', 'Mal_JK62', 'Mal_JK63', 'Mal_JK64', 'Mal_JK66', 'Mal_JK68', 'Mal_JK69', 'Mal_JK7-', 'Mal_JK70', 'Mal_JK71', 'Mal_JK72', 'Mal_JK73', 'Mal_JK74', 'Mal_JK75', 'Mal_JK79', 'Mal_JK8-', 'Mal_JK80', 'Mal_JK85', 'Mal_JK87', 'Mal_JK88', 'Mal_JK9-', 'Mal_JK90', 'Mal_JK92', 'Mal_JK93', 'Mal_JK96']
only_heterozygous_list = ['A50_CSFP', 'A51_CSFP', 'A52_CSFP', 'A56_CSFP', 'A58_CSFP', 'A5_CSFP2', 'A63_CSFP', 'A64_CSFP', 'B100_CSF', 'B101_CSF', 'B110_CSF', 'B112_CSF', 'B113_CSF', 'B22_CSFP', 'JK2_48_G', 'JK2_55_C', 'JK2_66_T', 'JK2_79_A', 'JK2_80_T', 'JK2_82_A', 'JK2_83_T', 'JK2_84_T', 'JK2_86_A', 'JK2_87_G', 'JK2_88_T', 'JK2_90_G', 'JK2_97_C', 'JK2_99_G', 'JK5_11-D', 'JK5_12-E', 'JK5_13-F', 'JK5_14-G', 'JK5_15-H', 'JK5_16-A', 'JK5_18-C', 'JK5_19-D', 'JK5_20-E', 'JK5_21-F', 'JK5_22-G', 'JK5_23-H', 'JK5_25-B', 'JK5_26-C', 'JK5_34-B', 'JK5_36-D', 'JK5_38-F', 'JK5_41-E', 'JK5_42-F', 'JK5_44-H', 'JK5_45-A', 'JK5_46-B', 'JK5_47-C', 'JK5_48-D', 'JK5_50-F', 'JK5_51-A', 'JK5_52-C', 'JK5_58-B', 'JK5_62-C', 'JK5_66-D', 'JK5_69-H', 'JK5_7-H4', 'JK5_70-A', 'JK5_74-E', 'JK5_76-G', 'JK5_77-H', 'JK5_78-A', 'JK5_8-A5', 'JK5_85-H', 'JK5_9-B5', 'JMT10_CK', 'JMT12_CK', 'JMT17_CK', 'JMT18_CK', 'JMT19_CK', 'JMT1_CKD', 'JMT20_CK', 'JMT21_CK', 'JMT22_CK', 'JMT23_CK', 'JMT24_CK', 'JMT29_CK', 'JMT2_CKD', 'JMT30_CK', 'JMT31_CK', 'JMT32_CK', 'JMT33_CK', 'JMT34_CK', 'JMT35_CK', 'JMT36_CK', 'JMT37_CK', 'JMT38_CK', 'JMT39_CK', 'JMT3_CKD', 'JMT40_CK', 'JMT41_CK', 'JMT42_CK', 'JMT43_CK', 'JMT44_CK', 'JMT45_CK', 'JMT46_CK', 'JMT47_CK', 'JMT48_CK', 'JMT4_CKD', 'JMT57_CK', 'JMT58_CK', 'JMT59_CK', 'JMT5_CKD', 'JMT62_CK', 'JMT63_CK', 'JMT64_CK', 'JMT65_CK', 'JMT66_CK', 'JMT67_CK', 'JMT68_CK', 'JMT69_CK', 'JMT6_CKD', 'JMT72_CK', 'JMT7_CKD', 'JMT8_CKD', 'JMTD_4_C', 'JT40_CKD', 'JT48_CKD', 'JT52_CKD', 'JT56_CKD', 'JT5_CKDN', 'Mal_JK15', 'Mal_JK18', 'Mal_JK20', 'Mal_JK25', 'Mal_JK34', 'Mal_JK65', 'Mal_JK83', 'Mal_JK89', 'Mal_JK91', 'Mal_JK95']
only_anuploid_list = ['JK2_50_C', 'JK2_67_G', 'JK2_69_A', 'JK2_6_TG', 'JMT75_CK', 'JMT88_CK', 'JMT93_CK', 'JMT95_CK', 'JMT96_CK', 'JMTD_11_', 'JT12_CKD', 'JT20_CKD', 'JT23_CKD', 'JT2_CKDN', 'JT30_CKD', 'JT32_CKD', 'JT34_CKD', 'JT35_CKD', 'JT42_Mer', 'JT43_CKD', 'JT51_CKD', 'JT60_CKD', 'JT6_CKDN', 'JT7_CKDN', 'Mal_JK10', 'Mal_JK13', 'Mal_JK22', 'Mal_JK32', 'Mal_JK33', 'Mal_JK36', 'Mal_JK37', 'Mal_JK39', 'Mal_JK67', 'Mal_JK76', 'Mal_JK77', 'Mal_JK78', 'Mal_JK81', 'Mal_JK84', 'Mal_JK86', 'Mal_JK94']
anu_het_list = ['JK2_18B_', 'JK2_3_CC', 'JMTD_14_', 'JMTD_20_', 'JMTD_24_', 'JMTD_7_C', 'JT18_CKD', 'JT19_CKD', 'JT26_CKD', 'JT41_CKD', 'JT46_CKD', 'JT47_CKD', 'JT54_Mer', 'JT9_CKDN', 'Mal_JK12', 'Mal_JK28']
chimeras = ['A53_CSFP', 'A57_CSFP', 'B108_CSF', 'B24_CSFP', 'B36_CSFP', 'JK2_76_C', 'JK2_81_C', 'JK2_85_T', 'JK2_89_A', 'JK5_1-B4', 'JK5_10-C', 'JK5_17-B', 'JK5_24-A', 'JK5_27-D', 'JK5_28-E', 'JK5_29-F', 'JK5_31-G', 'JK5_32-H', 'JK5_33-A', 'JK5_35-C', 'JK5_37-E', 'JK5_39-G', 'JK5_40-H', 'JK5_43-G', 'JK5_49-E', 'JK5_5-F4', 'JK5_53-D', 'JK5_57-E', 'JK5_6-G4', 'JK5_68-G', 'JK5_71-B', 'JK5_72-C', 'JK5_73-D', 'JK5_75-F', 'JK5_79-B', 'JK5_80-C', 'JK5_81-D', 'JK5_82-E', 'JK5_83-F', 'JK5_84-G', 'JK_92_CA', 'JMT11_CK', 'JMT13_CK', 'JMT14_CK', 'JMT15_CK', 'JMT16_CK', 'JMT25_CK', 'JMT26_CK', 'JMT27_CK', 'JMT28_CK', 'JMT60_CK', 'JMT61_CK', 'JMT70_CK', 'JMT71_CK', 'JMT9_CKD', 'JT10_CKD', 'JT11_CKD', 'JT1_CKDN', 'JT22_CKD', 'JT29_CKD', 'JT3_CKDN', 'JT45_CKD', 'JT58_CKD', 'JT59_CKD', 'JT8_CKDN', 'Mal_JK61', 'Mal_JK82']

genotypes = ["WT", 'ung1∆', "exo1-nd", "pol32∆", "exo1-ndpol32∆", "ung1∆ EV", "WT EV"]

def find_list_duplicates(list_1) -> list:
    '''
    This function takes a list and returns a list of unique elements that appear more than once in that list.
    '''
    from collections import Counter
    counts = Counter(list_1)
    duplicates = [item for item, count in counts.items() if count > 1]
    return duplicates

def check_sample_in_genotype_dict(sample: str, genotype_dict: dict) -> str | None:
    """
    Checks if a sample is present in the genotype dictionary.

    This function checks if a given sample is present in the provided genotype dictionary.
    If the sample is not present, it logs an informational message and returns None.
    If the sample is present, it returns the corresponding genotype.

    Args:
        sample (str): The sample to check.
        genotype_dict (dict): The dictionary of genotypes to check against.

    Returns:
        str or None: The genotype corresponding to the sample if present, None otherwise.
    """
    if sample not in genotype_dict:
        print(f"Sample {sample} is not in genotype_dict. Skipping...")
        return None
    else:
        return genotype_dict[sample]
    
def count_genotypes_in_list(genotypes: list, in_genotype: list) -> None:
    """
    Counts and prints the occurrences of each genotype in a list.

    This function iterates over a list of genotypes and for each genotype, it counts the number of occurrences 
    in a second list (in_genotype). It then prints the genotype and its count.

    Args:
        genotypes (list): The list of genotypes to count.
        in_genotype (list): The list in which to count the occurrences of each genotype.

    Returns:
        None
    """
    for genotype in genotypes:
        print(genotype, in_genotype.count(genotype))

def show_genotype_breakdown(sample_name_list: list, genotype_dict: dict, genotypes: list) -> None:
    """
    Shows the breakdown of genotypes for a given list of samples.

    This function takes a list of sample names, a dictionary mapping sample names to genotypes, 
    and a list of genotypes. It checks each sample name in the list against the dictionary, 
    and if the sample name is found, it appends the corresponding genotype to a list. 
    It then counts the occurrences of each genotype in the list and prints the counts.

    Args:
        sample_name_list (list): The list of sample names to check.
        genotype_dict (dict): The dictionary mapping sample names to genotypes.
        genotypes (list): The list of genotypes to count.

    Returns:
        None
    """
    in_genotype = []
    for sample in sample_name_list: 
        sample_genotype = check_sample_in_genotype_dict(sample, genotype_dict)
        if sample_genotype is not None:
            in_genotype.append(sample_genotype)
    count_genotypes_in_list(genotypes, in_genotype)

def combine_and_deduplicate_lists(lists: list) -> tuple:
    """
    Combines multiple lists into one and removes duplicates.

    This function takes a list of lists, combines them into a single list, and then removes any duplicate elements.
    It returns two lists: the combined list (with duplicates) and the deduplicated list.

    Args:
        lists (list of list): The list of lists to combine and deduplicate.

    Returns:
        tuple: A tuple containing the combined list (with duplicates) and the deduplicated list.
    """
    combined_list = sum(lists, [])
    deduplicated_list = list(set(combined_list))
    return combined_list, deduplicated_list

def check_list_lengths(original_list: list, deduplicated_list: list) -> None:
    """
    Checks if the lengths of the original and deduplicated lists are equal.

    This function takes two lists: the original list and a deduplicated version of it. 
    It checks if the lengths of the two lists are equal. If they are not, it logs a warning 
    and prints the duplicates in the original list.

    Args:
        original_list (list): The original list before deduplication.
        deduplicated_list (list): The list after deduplication.

    Returns:
        None
    """
    if len(original_list) != len(deduplicated_list):
        logging.warning(f"Length of list before deduplication: {len(original_list)}")
        logging.warning(f"Length of list after deduplication: {len(deduplicated_list)}")
        logging.warning("Lengths are not equal")
        logging.warning(f"Duplicates: {find_list_duplicates(original_list)}")

def check_genotypes_in_lists(lists_to_check: list, genotype_dict: dict, genotypes: list) -> None:
    """
    Checks the breakdown of genotypes in multiple lists.

    This function takes a list of tuples, where each tuple contains a list to check and a name for that list. 
    For each tuple, it prints the name of the list and then shows the breakdown of genotypes in the list.

    Args:
        lists_to_check (list of tuple): The list of tuples to check. Each tuple should contain a list and a name for that list.
        genotype_dict (dict): The dictionary mapping sample names to genotypes.
        genotypes (list): The list of genotypes to count.

    Returns:
        None
    """
    for list_to_check, list_name in lists_to_check:
        print(f"Breakdown of genotypes in {list_name}")
        show_genotype_breakdown(list_to_check, genotype_dict, genotypes)
        print("\n")

if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
    from settings import config as cfg
    from settings import genotype_dict as genotype_dict_master

    genotype_dict = genotype_dict_master[cfg["genotype_dict"]]

    non_haploids_list, non_haploids_list_dedup = combine_and_deduplicate_lists(
        [only_heterozygous_list, only_anuploid_list, anu_het_list, chimeras]
    )

    check_list_lengths(non_haploids_list, non_haploids_list_dedup)

    lists_to_check = [
        [haploid_list, "haploid_list"], 
        [non_haploids_list_dedup, "non_haploids_list_dedup"]
    ]

    check_genotypes_in_lists(lists_to_check, genotype_dict, genotypes)
