# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
from typing import Union, Literal, overload, TypeVar, List

from .ParentalTable import ParentalTable
from .SampleTable import SampleTable
from .OverlapSwitchesTable import OverlapSwitchesTable
from .SwitchesTable import SwitchesTable
from .decorators import fancy_status, time_elapsed
from .TableSNP import TableSNP

logger = logging.getLogger('__main__.' + __name__)

#starts are needed in case some chromosomes have 0 calls
#and still need to be drawn on a map
starts = pd.DataFrame( [['ref|NC_001133|', 0, "NaN", "NaN"],
                        ['ref|NC_001134|', 0, "NaN", "NaN"],
                        ['ref|NC_001135|', 0, "NaN", "NaN"],
                        ['ref|NC_001136|', 0, "NaN", "NaN"],
                        ['ref|NC_001137|', 0, "NaN", "NaN"],
                        ['ref|NC_001138|', 0, "NaN", "NaN"],
                        ['ref|NC_001139|', 0, "NaN", "NaN"],
                        ['ref|NC_001140|', 0, "NaN", "NaN"],
                        ['ref|NC_001141|', 0, "NaN", "NaN"],
                        ['ref|NC_001142|', 0, "NaN", "NaN"],
                        ['ref|NC_001143|', 0, "NaN", "NaN"],
                        ['ref|NC_001144|', 0, "NaN", "NaN"],
                        ['ref|NC_001145|', 0, "NaN", "NaN"],
                        ['ref|NC_001146|', 0, "NaN", "NaN"],
                        ['ref|NC_001147|', 0, "NaN", "NaN"],
                        ['ref|NC_001148|', 0, "NaN", "NaN"]],
                        columns=["Chromosome", "Region", "Reference","Allele"])
starts_switches = starts.copy()[["Chromosome", "Reference", "Allele"]].rename(columns={"Reference": "Position_1", "Allele": "Position_2"})

def printLog(
        string: str,
        log: bool
        ) -> None:
    '''Prints a string to the console if logging is enabled.

    Args:
        string (str): A string to be printed to the console.
        log (bool): A boolean indicating whether logging is enabled.

    Returns:
        None
    '''

    if log:
        print(string)

@overload
def create_table_object(
    file_data: pd.DataFrame,
    kind: Literal["parental"],
    sample_name: str
    ) -> ParentalTable: ...
@overload
def create_table_object(
    file_data: pd.DataFrame,
    kind: Literal["sample"],
    sample_name: str
    ) -> SampleTable: ...
@overload
def create_table_object(
    file_data: pd.DataFrame,
    kind: Literal["switches_overlap"],
    sample_name: str
    ) -> OverlapSwitchesTable: ...

def create_table_object(
        file_data: pd.DataFrame,
        kind: Literal["parental", "sample", "switches_overlap"],
        sample_name: str
        ) -> Union[ParentalTable, SampleTable, OverlapSwitchesTable]:
    '''
    Creates a table object of the appropriate type based on the specified kind.

    Args:
        file_data (pd.DataFrame): A pandas DataFrame containing the data to be used 
            to create the table object.
        kind (Literal["parental", "sample", "switches_overlap"]): A string literal 
            indicating the kind of table object to be created. Possible values are 
            "parental", "sample", and "switches_overlap".
        sample_name (str): A string representing the name of the sample.

    Returns:
        Union[ParentalTable, SampleTable, OverlapSwitchesTable]: A table object of 
            the appropriate type based on the specified kind.
    '''
    if kind == "parental":
        return ParentalTable(file_data, sample_name)
    elif kind == "sample":
        return SampleTable(file_data, sample_name)
    elif kind == "switches_overlap":
        return OverlapSwitchesTable(file_data, sample_name)
    
def read_csv_file(
        path: str,
        format: str
        ) -> pd.DataFrame:
    '''
    Reads a CSV or text file from the specified path and returns a pandas DataFrame.

    Args:
        path (str): A string representing the path to the file to be read.
        format (str): A string representing the format of the file to be read. 
            Possible values are "csv" and "txt".

    Returns:
        pd.DataFrame: A pandas DataFrame containing the data from the file.
    '''
    if format == "txt":
        return pd.read_table(path)
    else:
        return pd.read_csv(path)

@overload
def csvDataFrameImport(
    directory: str,
    kind: Literal["parental"],
    log: bool = False,
    format: str = "csv"
    ) -> list[ParentalTable]: ...
@overload
def csvDataFrameImport(
    directory: str, 
    kind: Literal["sample"], 
    log: bool = False, 
    format: str = "csv",
    genotype_dict: dict | None = None
    ) -> list[SampleTable]: ...
@overload
def csvDataFrameImport(
    directory: str, 
    kind: Literal["switches_overlap"], 
    log: bool = False, 
    format: str = "csv",
    genotype_dict: dict | None = None
    ) -> list[OverlapSwitchesTable]: ...

def csvDataFrameImport(
        directory: str,
        kind: Literal["parental", "sample", "switches_overlap"],
        log: bool = False,
        format: str = "csv",
        genotype_dict: dict | None = None
        ) -> Union[list[ParentalTable], list[SampleTable], list[OverlapSwitchesTable], list[ParentalTable|SampleTable|OverlapSwitchesTable]]:
    '''Imports CSV or TXT files containing variant calls and returns a list of 
    TableSNP objects.

    Args:
        directory (str): A string representing the path to the directory 
            containing the CSV or TXT files.
        kind (str): A string indicating the kind of TableSNP subclass to be 
            used. Possible values are "parental", "sample", and 
            "switches_overlap".
        log (bool, optional): A boolean indicating whether logging is enabled. 
            Defaults to False.
        format (str, optional): A string indicating the format of the input 
            files. Possible values are "csv" and "txt". Defaults to "csv".

    Returns:
        list: A list of TableSNP objects, where each object corresponds to a 
            CSV or TXT file in the input directory.
    '''
    frames_list_samples: List[Union[ParentalTable, SampleTable, OverlapSwitchesTable]] = []

    for filename in os.listdir(directory):
        if filename.endswith(f".{format}"):
            sample_name = filename[:8] #.split("-")[0] skip split for now...

            #skip the sample if it is not in the genotype_dict (if provided)
            if genotype_dict is not None and sample_name not in genotype_dict:
                logging.warning(f"Sample {sample_name} not found in genotype_dict. Skipping...")
                continue

            path = os.path.join(directory, filename)

            file_data = read_csv_file(path, format)
            table = create_table_object(file_data, kind, sample_name)
            if log:
                table.printRows()
            table.extractSNVs(log=log)

            frames_list_samples.append(table)
    

    logging.info(f"found {len(frames_list_samples)} samples in {directory}")
    logging.info("Extracted SNVs from all rows")

    #check that frames_list_samples is of the csvDataFrameImportType type), if not, raise an error
    return(frames_list_samples)

def removeCommonDataFrameRows(
        parent_df1: ParentalTable,
        parent_df2: ParentalTable
        ) -> list[TableSNP]:
    '''Removes rows present in both parents as a form of quality control and 
    returns a list of ParentalTable objects.

    Args:
        parent_df1 (ParentalTable): A ParentalTable object representing the 
            first parental genome.
        parent_df2 (ParentalTable): A ParentalTable object representing the 
            second parental genome.

    Returns:
        list: A list of ParentalTable objects, where each object corresponds to 
            a parental genome with the common rows removed.
    '''

    subtracted_df_list = []

    #first, remove all Heterozygous rows from both parents
    
    len_before_p1 = len(parent_df1.df)
    len_before_p2 = len(parent_df2.df)
    parent_df1.df = parent_df1.df[parent_df1.df['Zygosity'] != "Heterozygous"]
    parent_df2.df = parent_df2.df[parent_df2.df['Zygosity'] != "Heterozygous"]
    logging.info(f'Removed {len_before_p1 - len(parent_df1.df)} Heterozygous records from parent 1')
    logging.info(f'Removed {len_before_p2 - len(parent_df2.df)} Heterozygous records from parent 2')


    #show the common rows
    common_df = pd.merge(parent_df1.df[["Chromosome", "Region"]],parent_df2.df[["Chromosome", "Region"]], indicator=True, how='inner').drop('_merge', axis=1)
    logging.info(f'Found {len(common_df)} common rows between parents')
    logging.info(f'Common rows between parents are: {common_df}')

    df1 = (pd.merge(parent_df1.df[["Chromosome", "Region"]],parent_df2.df[["Chromosome", "Region"]], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1))
    df1 = ParentalTable(df1, "parent1")
    df2 = (pd.merge(parent_df2.df[["Chromosome", "Region"]],parent_df1.df[["Chromosome", "Region"]], indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1))
    df2 = ParentalTable(df2, "parent2")

    subtracted_df_list.append(df1)
    subtracted_df_list.append(df2)

    logging.info(f'Removed {len(parent_df1.df) - len(df1.df)} rows from parent 1')
    logging.info(f'Removed {len(parent_df2.df) - len(df2.df)} rows from parent 2')

    return subtracted_df_list

def createDataFrameWithCompleteParentalSNPs(
        parent_df1: ParentalTable, 
        parent_df2: ParentalTable
        ) -> pd.DataFrame:
    '''Creates a pandas DataFrame containing a mapping of SNP positions and 
    alleles to parents.

    Args:
        parent_df1 (ParentalTable): A ParentalTable object representing the 
            first parental genome.
        parent_df2 (ParentalTable): A ParentalTable object representing the 
            second parental genome.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing the following columns: 
            "Chromosome", "Region", "Reference", "Allele", and "Parent". Each 
            row represents a SNP position and allele in one of the parents.
    '''

    parent1_final=parent_df1.df
    parent1_final['Parent']='parent_1'
    parent2_final=parent_df2.df
    parent2_final['Parent']='parent_2'

    #inverse frames to show all SNP positions in both parents
    parent1_inverse=parent1_final.copy()
    parent2_inverse=parent2_final.copy()
    parent1_inverse=parent1_inverse.rename(columns={"Reference": "Allele", "Allele": "Reference"})[['Chromosome', 'Region', 'Reference', 'Allele', 'Parent']]
    parent1_inverse['Parent']='parent_2'
    parent2_inverse=parent2_inverse.rename(columns={"Reference": "Allele", "Allele": "Reference"})[['Chromosome', 'Region', 'Reference', 'Allele', 'Parent']]
    parent2_inverse['Parent']='parent_1'

    #create final DF

    finalFrames = [parent1_final,parent2_final,parent1_inverse,parent2_inverse]
    df_switches_pos = pd.concat(finalFrames)
    df_switches_pos = df_switches_pos.astype({'Region': int})
    df_switches_pos = df_switches_pos.sort_values(by=["Chromosome","Region"],ascending=(True,True),ignore_index=True)
    df_switches_pos = df_switches_pos.astype({'Region': int})

    logging.info(f'Length of parental SNPs db_df is: {len(df_switches_pos)} with {len(df_switches_pos)/2} SNP positions')

    return df_switches_pos

def clusterTableImport(
        path: str,
        end_position_fix: bool = False
        ) -> pd.DataFrame:
    '''Imports a cluster table from a TSV file and returns a pandas DataFrame 
    with the relevant columns.

    Args:
        path (str): A string representing the path to the TSV file.
        end_position_fix (bool, optional): A boolean indicating whether to fix 
            the end position of each cluster to be the same as the start 
            position. Defaults to False.

    Returns:
        pandas.DataFrame: A pandas DataFrame with the following columns: 
            "Tumor_Sample_Barcode", "Variant_Type", "Chromosome", 
            "Start_position", "End_position", "StrainCluster_ID", 
            "Dataset_Cluster_ID", and "Cluster_Length". Each row represents a 
            cluster of variants in a sample.
    '''

    df = pd.read_csv(path, sep='\t')

    if end_position_fix:
        df["End_position"] = df["Start_position"]
    df = df[["Tumor_Sample_Barcode","Variant_Type","Chromosome","Start_position","End_position","StrainCluster_ID","Dataset_Cluster_ID","Cluster_Length"]]
    df = df[df['Dataset_Cluster_ID'].notna()]
    df = df.groupby(["Tumor_Sample_Barcode", "Dataset_Cluster_ID", "StrainCluster_ID", "Chromosome","Cluster_Length"], as_index=False) \
             .agg(cluster_start=("Start_position", "min"), cluster_end=("End_position", "max"))
    
    logging.info(f"Imported {len(df)} PMACD clusters among {df['Tumor_Sample_Barcode'].unique()} samples from {path}")

    return df

def renameChrToRoman(
        df: pd.DataFrame,
        chromosome_column_name: str, 
        naming_scheme: str = "S288C"
        ) -> pd.DataFrame:
    '''Renames chromosome names in a pandas DataFrame from numeric or S288C 
    naming schemes to Roman numerals.

    Args:
        df (pandas.DataFrame): A pandas DataFrame containing the chromosome 
            names to be renamed.
        chromosome_column_name (str): A string representing the name of the 
            column containing the chromosome names.
        naming_scheme (str, optional): A string indicating the naming scheme 
            used in the input DataFrame. Possible values are "numeric" and
            "S288C". Defaults to "S288C".

    Returns:
        pandas.DataFrame: A pandas DataFrame with the chromosome names renamed 
            to Roman numerals.
    '''
    from mayo.settings import numeric_to_roman
    from mayo.settings import S288C_to_roman

    numeric_to_roman: dict[str,str] = {
        "chr1": "I",
        "chr2": "II",
        "chr3": "III",
        "chr4": "IV",
        "chr5": "V",
        "chr6": "VI",
        "chr7": "VII",
        "chr8": "VIII",
        "chr9": "IX",
        "chr10": "X",
        "chr11": "XI",
        "chr12": "XII",
        "chr13": "XIII",
        "chr14": "XIV",
        "chr15": "XV",
        "chr16": "XVI"}

    S288C_to_roman : dict[str,str] = {
        'ref|NC_001133|':"I",
        'ref|NC_001134|':"II",
        'ref|NC_001135|':"III",
        'ref|NC_001136|':"IV",
        'ref|NC_001137|':"V",
        'ref|NC_001138|':"VI",
        'ref|NC_001139|':"VII",
        'ref|NC_001140|':"VIII",
        'ref|NC_001141|':"IX",
        'ref|NC_001142|':"X",
        'ref|NC_001143|':"XI",
        'ref|NC_001144|':"XII",
        'ref|NC_001145|':"XIII",
        'ref|NC_001146|':"XIV",
        'ref|NC_001147|':"XV",
        'ref|NC_001148|':"XVI"}
    
    if naming_scheme == "numeric":
        df.replace({chromosome_column_name: numeric_to_roman}, inplace=True)

    elif naming_scheme == "S288C":
        df.replace({chromosome_column_name: S288C_to_roman}, inplace=True)

    else:
        logging.error(f"Naming scheme '{naming_scheme}' is not recognized")

    return(df)

def get_dens_regions(
        SNP_df: pd.DataFrame,
        max_distance: int = 500
        ) -> dict:
    """
    Given a DataFrame of SNPs and a maximum distance, returns a dictionary of regions where SNPs are close together.

    Args:
        SNP_df (pd.DataFrame): A DataFrame containing SNP data, with columns 'Chromosome' and 'Region'.
        max_distance (int, optional): The maximum distance between SNPs to consider them part of the same region. Defaults to 500.

    Returns:
        dict: A dictionary where the keys are chromosome names and the values are lists of contiguous regions where SNPs are close together.
    """
    SNP_df = SNP_df.sort_values(by=['Chromosome', 'Region'])
    chromosomes = SNP_df['Chromosome'].unique()

    # create a dictionary to store the regions where SNPs are close together
    close_regions = {}

    # loop through each chromosome
    for chromosome in chromosomes:
        df_chrom = SNP_df[SNP_df['Chromosome'] == chromosome]
        positions = df_chrom['Region'].tolist()
        close_regions[chromosome] = []
        # loop through each position
        current_pos = positions[0]
        in_contig = False
        contig_start = None

        for i in range(len(positions) - 1):
            # check if the next position is close enough to the current position
            distance = positions[i+1] - current_pos

            # if the next position is close enough, it is part of the same contig
            if distance <= max_distance and contig_start is None:
                in_contig = True
                contig_start = current_pos

            # if the next contig position is not close enough, but we are in a contig, save the contig
            elif distance > max_distance and in_contig:
                # save the contig
                close_regions[chromosome].append([contig_start, current_pos])
                contig_start = None
                in_contig = False
                
            current_pos = positions[i+1]

        # if we are still in a contig at the end of the chromosome, save the contig
        if in_contig:
            close_regions[chromosome].append([contig_start, current_pos])

    return close_regions

def get_complement_regions(regions: dict) -> dict:
    from mayo.settings.chromosome_lengths import chromosomes
    chromosomes = {chromosome[0]: chromosome[1] for chromosome in chromosomes} #convert chromosomes to a dictionary
    complement_regions = {}

    for chromosome in chromosomes:
        chromosome_length = chromosomes[chromosome]

        #if there are no regions on this chromosome, add a single region that covers the whole chromosome
        if chromosome not in regions or len(regions[chromosome]) == 0:
            complement_regions[chromosome] = [[0, chromosome_length]]
            continue
        
        #if there are regions on this chromosome, create a list to store the complement regions
        list_complement_regions = []
        chromosomal_regions = regions[chromosome]

        #make sure the regions are sorted
        chromosomal_regions.sort(key=lambda x: x[0])

        for index, region in enumerate(chromosomal_regions):
            region_start = region[0]
            region_end = region[1]

            #if the first region starts after the start of the chromosome, add a region from the start of the chromosome to the start of the first region
            if index == 0 and region_start > 0:
                list_complement_regions.append([0, region_start])

            #if there is a next region, add a region from the end of the current region to the start of the next region
            if index < len(regions[chromosome]) - 1:
                list_complement_regions.append([region_end, regions[chromosome][index+1][0]])
            
            #if there is no next region, add a region from the end of the current region to the end of the chromosome
            else:
                list_complement_regions.append([region_end, chromosome_length])
            
        complement_regions[chromosome] = list_complement_regions

    return complement_regions

def loadDensityMap(
        lowDensity: bool = False, 
        max_distance: int = 500, 
        data_path: str = 'outputs/all_parental_SNPs.tsv'
        ) -> dict:
    """
    Load a density map from a TSV file.

    Args:
        lowDensity (bool): If True, return a low-density map. If False, return a high-density map. Defaults to False.
        max_distance (int): The maximum distance between SNPs to consider when generating the density map. Defaults to 500.
        data_path (str): The path to the TSV file containing the density map data. Defaults to 'outputs/all_parental_SNPs.tsv'.

    Returns:
        dict: A dictionary containing the density map data.
    """
    df = pd.read_csv(data_path, sep='\t', header=0)
    df = df[df['Parent'] == "parent_1"]

    logging.info(df['Chromosome'].value_counts())
    logging.info(df['Zygosity'].value_counts())
    logging.info(df['Parent'].value_counts())
    logging.info(df[df['Zygosity'] != "Homozygous"])
    
    highDensityMap = get_dens_regions(df, max_distance=max_distance)

    if lowDensity:
        return get_complement_regions(highDensityMap)
    else:
        return highDensityMap

def plot_chr_SNP_coverage_by_SNP_distance(
        df: pd.DataFrame,
        distance_max: int = 2000,
        distace_step: int = 10,
        inverse: bool = False,
        by_chromosome: bool = False,
        save_path: str = 'outputs/SNP_density.tsv',
        fig_path: str = 'figures/SNP_density.png',
        show: bool = False
        ) -> None:
    """
    Given a DataFrame of SNP data, plots the fraction of each chromosome covered by SNP-dense regions as a function of the maximum distance between SNPs.

    Args:
        df (pd.DataFrame): A DataFrame containing SNP data, with columns 'Chromosome' and 'Region'.
        distance_max (int, optional): The maximum distance between SNPs to consider them part of the same region. Defaults to 2000.
        distace_step (int, optional): The step size to use when iterating over maximum distances. Defaults to 10.
        save_path (str, optional): The file path to save the output data to. Defaults to 'outputs/SNP_density.tsv'.

    Returns:
        None
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    from mayo.settings.chromosome_lengths import chromosomes
    from mayo.settings import S288C_to_roman

    chromosomes = {chromosome[0]: chromosome[1] for chromosome in chromosomes} #convert chromosomes to a dictionary
    df_density = pd.DataFrame(columns=['chr', 'max_distance', 'density'])

    for max_distance in range(0, distance_max, distace_step):
        logging.info(f"Max distance: {max_distance}")

        regions = get_dens_regions(df, max_distance=max_distance)

        if inverse:
            regions = get_complement_regions(regions)
        
        for chromosome in regions:
            if chromosome not in chromosomes:
                continue
            total_length = 0
            for contig in regions[chromosome]:
                total_length += contig[1] - contig[0]

            df_density = pd.concat([df_density, pd.DataFrame([[chromosome, max_distance, total_length/chromosomes[chromosome]]], columns=['chr', 'max_distance', 'density'])], ignore_index=True)
    
    df_density['chr'] = df_density['chr'].map(S288C_to_roman)

    #plot the curve max_distance (x) vs density (y) for each chromosome
    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    sns.set_context("talk")
    if by_chromosome: 
        sns.lineplot(data=df_density, x='max_distance', y='density', hue='chr')
    else:
        sns.lineplot(data=df_density, x='max_distance', y='density')

    plt.title("SNP coverage by distance", fontsize=20, fontweight='bold')
    plt.xlabel("Max distance between SNPs", fontsize=16, fontweight='bold')

    plt.ylabel("Fraction of chromosome covered by SNP-dense regions", fontsize=14, fontweight='bold')
    if inverse:
        plt.ylabel("Fraction of chromosome covered by SNP-poor regions", fontsize=14, fontweight='bold')

    #make legent small
    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., title='Chromosome')

    plt.subplots_adjust(right=0.80) #make room for the legend
    plt.savefig(fig_path, dpi=300)
    logging.info(f"Saved SNP density plot to {fig_path}")

    if show:
        plt.show()
    
    else:
        plt.close()

    #save the data
    df_density.to_csv(save_path, sep='\t', index=False)
    logging.info(f"Saved SNP density data to {save_path}")

def findSpectraStrandwise(colA: str, colB: str) -> str:
    """
    Given two nucleotide bases, returns the type of substitution that occurs when
    transitioning from the first base to the second base.

    Args:
        colA (str): The first nucleotide base.
        colB (str): The second nucleotide base.

    Returns:
        str: The type of substitution that occurs when transitioning from colA to colB.
    """
    if ((colA == 'C') & (colB == 'T')):
        return('C_to_T')
    elif ((colA == 'C') & (colB == 'A')):
        return('C_to_A')
    elif ((colA == 'C') & (colB == 'G')):
        return('C_to_G')
    elif ((colA == 'T') & (colB == 'C')):
        return('T_to_C')
    elif ((colA == 'T') & (colB == 'G')):
        return('T_to_G')
    elif ((colA == 'T') & (colB == 'A')):
        return('T_to_A')
    
    elif ((colA == 'G') & (colB == 'A')):
        return('G_to_A')
    elif ((colA == 'G') & (colB == 'T')):
        return('G_to_T')
    elif ((colA == 'G') & (colB == 'C')):
        return('G_to_C')
    elif ((colA == 'A') & (colB == 'G')):
        return('A_to_G')
    elif ((colA == 'A') & (colB == 'C')):
        return('A_to_C')
    elif ((colA == 'A') & (colB == 'T')):
        return('A_to_T')
    
def add_legend_to_ax(ax, pval_list, cmap_dict):
    import matplotlib.patches as mpatches
    legend_handles = []
    for pval in pval_list:
        color = cmap_dict[pval]
        patch = mpatches.Patch(color=color, label=f'p<{pval:.1e}')
        legend_handles.append(patch)

    new_legend = ax.legend(handles=legend_handles, title='Cluster p-value', loc='center right', fontsize=10, title_fontsize=12, bbox_to_anchor=(0.98, 0.5))
    ax.add_artist(new_legend)

def add_clusters_patch_to_ax(ax, pval_list, SampleTable_Object, clusters_kind="JT"):

    from mayo.settings import S288C_to_numeric
    import numpy as np

    #assign a descrete colormap for pvalues from the seaborn's colormap
    #generate a list of colors from the "rocket" colormap
    cmap = sns.color_palette("rocket", as_cmap=True)
    cmap_colors = cmap(np.linspace(0, 1, len(pval_list)))
    #convert the list of colors to a dictionary
    cmap_dict = {pval_list[i]: cmap_colors[i] for i in range(len(pval_list))}

    pval_counter = 0
    for pval in pval_list:

        if clusters_kind == "JT":
            cluster_df = SampleTable_Object.cluster_dict_JT[pval].copy()
            
            #add a "Sample_Cluster_ID" column to the cluster_df from index
            cluster_df["Sample_Cluster_ID"] = cluster_df.index

            #convert the chromosome names in the cluster_df["Chromosome"] column to numeric
            cluster_df["Chromosome"] = cluster_df["Chromosome"].replace(S288C_to_numeric)

            cluster_dict = {}
            for index, row in cluster_df.iterrows():
                #get the Dataset_Cluster_ID, Chromosome, Start_position, Cluster_Length
                Cluster_ID = row["Sample_Cluster_ID"]
                Chromosome = row["Chromosome"]
                Start_position = row["Start"]
                Cluster_Length = row["Length"]
                
                #add the cluster to the cluster_dict
                cluster_dict[Cluster_ID] = [Chromosome, Start_position, Cluster_Length]
        
        elif clusters_kind == "PMACD":
            cluster_df = SampleTable_Object.cluster_dict_PMACD[pval].copy()

            cluster_dict = {}
            for index, row in cluster_df.iterrows():
                #get the Dataset_Cluster_ID, Chromosome, Start_position, Cluster_Length
                Dataset_Cluster_ID = row["StrainCluster_ID"]
                Chromosome = row["Chromosome"]
                Start_position = row["Start_position"]
                Cluster_Length = row["Cluster_Length"]
                
                #add the cluster to the cluster_dict
                cluster_dict[Dataset_Cluster_ID] = [Chromosome, Start_position, Cluster_Length]

        for key in cluster_dict:
            cluster_chromosome = cluster_dict[key][0]

            #check if the cluster_chromosome returned is a string
            if not isinstance(cluster_chromosome, str):
                logging.warning(f"Cluster chromosome is not a string. Skipping cluster {key}. for sample {SampleTable_Object.name}.")
                continue
            
            cluster_chromosome = cluster_chromosome.replace("chr", "")
            cluster_chromosome = int(cluster_chromosome)
            cluster_chromosome_coord = cluster_chromosome - (1.45) + (pval_counter*0.04) #(1.5 - offset to move rectangle lower/higher)
            cluster_start = cluster_dict[key][1]
            cluster_length = cluster_dict[key][2]
            
            #print(f"Drawing cluster {key} with pvalue {pval} and color {cmap_dict[pval]} at {cluster_chromosome}:{cluster_start}-{cluster_start+cluster_length}.")
            ax.add_patch(plt.Rectangle((cluster_start+500, cluster_chromosome_coord), cluster_length+500, 0.5, edgecolor='black', linewidth=0.1, fill=True, facecolor=cmap_dict[pval], alpha=1, zorder=-2))
        pval_counter += 1
        
    add_legend_to_ax(ax, pval_list, cmap_dict)

def drawSNPMap(
        SampleTable_Object,
        df_chr_lengths: pd.DataFrame,
        chr_starts_df: pd.DataFrame,
        saveMap: bool = True,
        map_type: str = "switches",
        showMap: bool = False,
        save_name_suffix: str = "",
        cluster_kind: str | None = None
        ) -> None:
    '''Draws a SNP map for a given SampleTable object and saves it as a PNG file.

    Args:
        SampleTable_Object (SampleTable): A SampleTable object representing the 
            sample to be plotted.
        df_chr_lengths (pandas.DataFrame): A pandas DataFrame containing the 
            lengths of each chromosome.
        chr_starts_df (pandas.DataFrame): A pandas DataFrame containing the 
            start positions of each chromosome.
        saveMap (bool, optional): A boolean indicating whether to save the map 
            as a PNG file. Defaults to True.
        map_type (str, optional): A string indicating the type of SNP map to 
            draw. Possible values are "switches", "SNPs", "SNVs_A3A_cluster", and "SNVs_A3A_all"
            Defaults to "switches".
        showMap (bool, optional): A boolean indicating whether to show the map 
            in a window. Defaults to False.

    Returns:
        None
    '''

    # style with ticks
    sns.set_style("ticks")
    sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
    
    fig, ax = plt.subplots()
    title = f"{SampleTable_Object.name}" #{round(SampleTable_Object.percent,3)}
    
    df_chr_lengths = renameChrToRoman(df_chr_lengths, "chromosome")
    #use #E6E6E6 for gray color
    df_chr_lengths.plot(kind='barh',legend=False, ax=ax, x="chromosome", y="end_position", fontsize=15, figsize=(20,6), edgecolor="black", linewidth=1, color="#E6E6E6" , zorder=2, title=title, label="Position")

    #To rename legend elements, change plot markers/colors, modify here
    label_dict = {"REF": "Reference Allele"}
    
    if map_type=="switches":
        palette = {"parent_1": "tab:green", "parent_2": "tab:red"}
        markers = {"parent_1": "|", "parent_2": "|"}

        pd_df = SampleTable_Object.output_df.copy()
        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True)
        pd_df = renameChrToRoman(pd_df, "Chromosome")
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", hue="Parent", palette=palette, style="Parent", markers=markers, s=120, alpha=0.5,zorder=3, linewidth=1.5)

    #renamed from SNPs to SNVs
    elif map_type=="SNVs":
        pd_df = SampleTable_Object.df.copy()
        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True)
        pd_df = renameChrToRoman(pd_df, "Chromosome")
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", marker="|", s=120, alpha=0.5,zorder=3, linewidth=1.5)
    
    elif map_type=="SNVs_A3A_no_switch":
        pd_df = SampleTable_Object.df.copy()
        pd_df = pd_df[(((pd_df["Reference"] == "C") & (pd_df["Allele"] == "T")) | ((pd_df["Reference"] == "G") & (pd_df["Allele"] == "A")))]
        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True).sort_values(by=["Chromosome", "Region"])
        pd_df = renameChrToRoman(pd_df, "Chromosome")
        palette = {"C": "tab:blue", "G": "tab:red", 'NaN': "tab:gray", "T": "tab:gray", "A": "tab:gray"}
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", marker="|", s=120, alpha=1, zorder=3, linewidth=1.5, hue="Reference", palette=palette)

    elif map_type=="SNVs_A3A_cluster":
        pd_df = SampleTable_Object.df_SNPs.copy()
        pd_df = pd_df[(((pd_df["Reference"] == "C") & (pd_df["Allele"] == "T")) | ((pd_df["Reference"] == "G") & (pd_df["Allele"] == "A")))]
        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True).sort_values(by=["Chromosome", "Region"])
        pd_df = renameChrToRoman(pd_df, "Chromosome")
        palette = {"C": "tab:blue", "G": "tab:red", 'NaN': "tab:gray", "T": "tab:gray", "A": "tab:gray"}
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", marker="|", s=120, alpha=1, zorder=3, linewidth=1.5, hue="Reference", palette=palette)

        switches_df = SampleTable_Object.switchesTable.copy()
        switches_df = switches_df[switches_df["Chromosome"] != "ref|NC_001224|"]
        switches_df = pd.concat([switches_df, starts_switches], ignore_index=True)
        switches_df = renameChrToRoman(switches_df, "Chromosome")
        sns.scatterplot(ax=ax,data=switches_df, x="Switch_Center", y="Chromosome", color="black", marker=2, s=120, zorder=0, linewidth=2.5, alpha=0.8)

    elif map_type=="SNVs_A3A":
        pd_df = SampleTable_Object.df.copy()
        pd_df = pd_df[(((pd_df["Reference"] == "C") & (pd_df["Allele"] == "T")) | ((pd_df["Reference"] == "G") & (pd_df["Allele"] == "A")))]
        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True).sort_values(by=["Chromosome", "Region"])
        pd_df = renameChrToRoman(pd_df, "Chromosome")
        palette = {"C": "tab:blue", "G": "tab:red", 'NaN': "tab:gray", "T": "tab:gray", "A": "tab:gray"}
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", marker="|", s=120, alpha=1, zorder=3, linewidth=1.5, hue="Reference", palette=palette)

        switches_df = SampleTable_Object.switchesTable.copy()
        switches_df = switches_df[switches_df["Chromosome"] != "ref|NC_001224|"]
        switches_df = pd.concat([switches_df, starts_switches], ignore_index=True)
        switches_df = renameChrToRoman(switches_df, "Chromosome")
        sns.scatterplot(ax=ax,data=switches_df, x="Switch_Center", y="Chromosome", color="black", marker=2, s=120, zorder=0, linewidth=2.5, alpha=0.8)

    elif map_type=="combined":
        import numpy as np
        from natsort import index_natsorted
        from mayo.settings import centromeres as chr_centromeres

        #first deal with the SNPs (swaps of parent 1 and parent 2)
        pd_df = SampleTable_Object.parentalSNPs.copy()
        pd_df_heterozygous = SampleTable_Object.df_SNPs_heterozygous.copy()
        pd_df_heterozygous = pd.concat([pd_df_heterozygous, chr_starts_df], ignore_index=True)
        pd_df_heterozygous = pd_df_heterozygous.sort_values(by="Chromosome", key=lambda x: np.argsort(index_natsorted(pd_df_heterozygous["Chromosome"])))
        pd_df_heterozygous = renameChrToRoman(pd_df_heterozygous, "Chromosome")

        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True)
        pd_df = pd_df.sort_values(by="Chromosome", key=lambda x: np.argsort(index_natsorted(pd_df["Chromosome"])))  
        pd_df = renameChrToRoman(pd_df, "Chromosome")

        palette = {"parent_1": sns.color_palette("deep")[2], "parent_2": sns.color_palette("deep")[3]}
        #palette = {"parent_1": "tab:green", "parent_2": "tab:red"}
        markers = {"parent_1": 3, "parent_2": 3}

        sns.scatterplot(ax=ax, data=pd_df_heterozygous, x="Region", y="Chromosome", color="tab:gray", marker=3, s=60, alpha=0.6,zorder=-1, linewidth=0.75)
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", hue="Parent", palette=palette, style="Parent", markers=markers, s=70, alpha=0.65,zorder=-1, linewidth=0.75)

        #second deal with the SNVs (A3A)
        label_dict = {
            "cen": 'Centromere',
            "SPECTRA_STRANDWISE": "Reference Allele",
            "Position": "Position",
            "G": "G→N",
            "C": "C→N",
            "A": "A→N",
            "T": "T→N",
            "C_to_T": "C→T",
            "C_to_A": "C/G→V/B",
            "C_to_G": "C/G→V/B",

            "G_to_A": "G→A",
            "G_to_C": "C/G→V/B",
            "G_to_T": "C/G→V/B",

            'T_to_A': "A/T→N",
            'T_to_C': "A/T→N",
            'T_to_G': "A/T→N",
            'A_to_C': "A/T→N",
            'A_to_G': "A/T→N",
            'A_to_T': "A/T→N",
            }
        markers = {
            'cen': "o",
            "C_to_T": "$|$",
            "C_to_A": "$|$",
            "C_to_G": "$|$",

            "G_to_A": "$|$",
            "G_to_C": "$|$",
            "G_to_T": "$|$",

            'T_to_A': "$|$",
            'T_to_C': "$|$",
            'T_to_G': "$|$",
            'A_to_C': "$|$",
            'A_to_G': "$|$",
            'A_to_T': "$|$",
            }
        palette = {
            'cen': "white",
            "C_to_T": "tab:red",
            "C_to_A": "tab:olive",
            "C_to_G": "tab:olive",

            "G_to_A": "tab:blue",
            "G_to_C": "tab:olive",
            "G_to_T": "tab:olive",

            'A_to_C': "tab:gray",
            'A_to_G': "tab:gray",
            'A_to_T': "tab:gray",

            'T_to_A': "tab:gray",
            'T_to_C': "tab:gray",
            'T_to_G': "tab:gray"
            }

        df_SNVs = SampleTable_Object.df.copy()

        try:
            df_SNVs["SPECTRA_STRANDWISE"] = df_SNVs.apply(lambda x: findSpectraStrandwise(x["Reference"], x["Allele"]), axis=1)
        except:
            logging.error(f"Error in findSpectraStrandwise function. Skipping sample {SampleTable_Object.name}. The sample dataframe might be empty (size: {len(df_SNVs)}).")
            plt.close()
            return


        #df_SNVs = pd.concat([df_SNVs, chr_starts_df], ignore_index=True).sort_values(by=["Chromosome", "Region"])
        df_SNVs = renameChrToRoman(df_SNVs, "Chromosome")
        chr_centromeres = renameChrToRoman(chr_centromeres, "Chromosome", "numeric")
        df_SNVs = pd.concat([df_SNVs, chr_centromeres], ignore_index=True)
        df_SNVs = df_SNVs.sort_values(by="Chromosome", key=lambda x: np.argsort(index_natsorted(df_SNVs["Chromosome"])))    

        sns.scatterplot(ax=ax, data=df_SNVs[df_SNVs["Reference"].isin(["cen"])], x="Region", y="Chromosome", hue="SPECTRA_STRANDWISE", palette=palette, style="SPECTRA_STRANDWISE", markers=markers, s=100, alpha=1,zorder=2, linewidth=0.75, edgecolor="black")
        sns.scatterplot(ax=ax, data=df_SNVs[~df_SNVs["Reference"].isin(["cen"])], x="Region", y="Chromosome", hue="SPECTRA_STRANDWISE", palette=palette, style="SPECTRA_STRANDWISE", markers=markers, s=100, alpha=0.9,zorder=3, linewidth=0.05, edgecolor="black")
        #sns.scatterplot(ax=ax, data=df_SNVs, x="Region", y="Chromosome", marker="|", s=120, alpha=1, zorder=3, linewidth=1.5, hue="Reference", palette=palette)

        switches_df = SampleTable_Object.switchesTable.copy()
        switches_df = switches_df[switches_df["Chromosome"] != "ref|NC_001224|"]
        switches_df = pd.concat([switches_df, starts_switches], ignore_index=True)
        switches_df = renameChrToRoman(switches_df, "Chromosome")
        sns.scatterplot(ax=ax,data=switches_df, x="Switch_Center", y="Chromosome", color="black", marker=3, s=90, zorder=-2, linewidth=0.75, alpha=0.70)
        
        df_recEvents = RecEventsTable_from_GC_CO_dict(SampleTable_Object.GC_dict, SampleTable_Object.CO_dict)
        #remove heterozygous chromosomes from the recombination events; they are not informative/reliable there
        heterozygous_chr = SampleTable_Object.heterozygous_chromosomes

        if len(df_recEvents) != 0:
            df_recEvents = df_recEvents[~df_recEvents["Chromosome"].isin(heterozygous_chr)]
            df_recEvents = pd.concat([df_recEvents, chr_starts_df], ignore_index=True)
            df_recEvents = renameChrToRoman(df_recEvents, "Chromosome")
            palette_recEvents = {"GC": "tab:blue", "CO": "tab:red"}
            sns.scatterplot(ax=ax, data=df_recEvents, x="Region", y="Chromosome", hue="Type", palette=palette_recEvents, marker=11, s=16, zorder=2, linewidth=0.5, alpha=0.7)
            add_GC_CO_legend_to_ax(ax, palette_recEvents)

        if cluster_kind is None:
            pass

        else:
            #pval_list must be ordered from largest to smallest pvalue to make sense in the legend

            if cluster_kind == "PMACD":
                pval_list = ['p05', 'p01', 'p005', 'p001_def'] #'p0005', #'p0001']
                pval_list.reverse()
            
            elif cluster_kind == "JT":
                pval_list = [0.01, 0.001, 0.0001, 0.00001]
                pval_list.reverse()

            else:
                raise ValueError(f"Cluster kind {cluster_kind} is not supported. Please choose from 'JT' or 'PMACD'.")

            save_name_suffix += f"_{cluster_kind}_clusters"
            add_clusters_patch_to_ax(ax, pval_list, SampleTable_Object, clusters_kind=cluster_kind)

    else: #just clusters
        pd_df = SampleTable_Object.df_SNPs.copy()

        pd_df = pd.concat([pd_df, chr_starts_df], ignore_index=True).sort_values(by=["Chromosome", "Region"])
        pd_df = renameChrToRoman(pd_df, "Chromosome")
        palette = {True: "blue", False: "red"}
        sns.scatterplot(ax=ax, data=pd_df, x="Region", y="Chromosome", marker="|", s=120, alpha=1, zorder=3, linewidth=1.5, hue='Clustered_Proximal', palette=palette) #hue='Switch_Proximal', Is_Clustered, Clustered_Proximal

        switches_df = SampleTable_Object.switchesTable.copy()
        switches_df = switches_df[switches_df["Chromosome"] != "ref|NC_001224|"]
        switches_df = pd.concat([switches_df, starts_switches], ignore_index=True)
        switches_df = renameChrToRoman(switches_df, "Chromosome")
        sns.scatterplot(ax=ax,data=switches_df, x="Switch_Center", y="Chromosome", color="black", marker=2, s=120, zorder=0, linewidth=1.5, alpha=0.8)

    ax.set_xlabel("POSITION")
    ax.set_ylabel("CHROMOSOME")

    if map_type=="combined":

        #dont repeat the same item in the legend
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        #order the legend
        legend_order = ["cen", "Position", "C_to_T", "G_to_A", "C_to_A", "C_to_G", "G_to_C", "G_to_T", "T_to_A", "T_to_C", "T_to_G", "A_to_C", "A_to_G", "A_to_T"]

        #relable the legend texts according to the legend_order

        while True:
            try:
                by_label = {label_dict.get(key): by_label[key] for key in legend_order}
                break

            except Exception as e:
                try:
                    logging.warning(f"Could not relabel legend. {e}. Removing {e} from legend_order and trying again.")
                    #extract the string from the error message
                    e = str(e).split("'")[1]
                    legend_order = [x for x in legend_order if x != str(e)]
                except Exception as e:
                    logging.warning(f"Cound not relabel legend. {e}.")
                    break

        ax.legend(by_label.values(), by_label.keys(), fontsize=10, loc='lower right', bbox_to_anchor=(0.98, 0.05), title=None)
        
        plt.tight_layout()

    else:

        legend_texts = plt.legend().get_texts()
        for i in range(len(legend_texts)):
            label_str = legend_texts[i].get_text()
            #remove NaNs from legend
            if label_str == "NaN":
                legend_texts[i].set_text("")
            if label_str in label_dict:
                new_label = label_dict.get(label_str)

                legend_texts[i].set_text(new_label)
    
    if saveMap:
        if not os.path.exists(f'figures/tests/{map_type}'):
            os.makedirs(f'figures/tests/{map_type}')
        plt.savefig(f'figures/{map_type}/{title}_{map_type}{save_name_suffix}.png', transparent=False, dpi=600)

    if showMap:
        plt.show()

    else:
        # close the plot
        plt.close()

def countClusteredSNPs(
        df_SNPs: pd.DataFrame, 
        df_Switches, 
        window_size: int
        ) -> pd.DataFrame:
    '''Counts the number of SNPs within a given window size of each switch in a 
    pandas DataFrame. Not currently used.

    Args:
        df_SNPs (pandas.DataFrame): A pandas DataFrame containing SNP 
            information, with columns "Chromosome" and "Region".
        df_Switches (pandas.DataFrame): A pandas DataFrame containing switch 
            information, with columns "Chromosome" and "Switch_Center".
        window_size (int): An integer representing the size of the window around 
            each switch to count SNPs in.

    Returns:
        pandas.DataFrame: A pandas DataFrame with the same columns as 
            df_Switches, plus two additional columns: "Closeby_SNPs" (the number 
            of SNPs within the window) and "Window_Size" (the size of the window).
    '''

    df_Switches["Closeby_SNPs"] = 0
    df_Switches["Window_Size"] = window_size

    for index, row in df_Switches.iterrows():
        chromosome = row["Chromosome"]
        bound_left = row["Switch_Center"] - window_size
        bound_right = row["Switch_Center"] + window_size

        for snp_index, snp_row in df_SNPs.iterrows():
            position=snp_row["Region"]
            snp_chromosome=snp_row["Chromosome"]

            if (snp_chromosome == chromosome):
                if (position > bound_left) & (position < bound_right):
                    df_Switches.loc[index, "Closeby_SNPs"] += 1

    return(df_Switches)

def countA3ASNPs(df: pd.DataFrame) -> tuple[int, int]:
    '''Counts the number of A3A SNPs in a pandas DataFrame and returns the 
    counts as a tuple. Note: here, A3A SNPs are defined as SNPs where the
    reference allele is C and the alternate allele is T, or where the reference
    allele is G and the alternate allele is A. Genomic context is not taken into
    account.

    Args:
        df (pandas.DataFrame): A pandas DataFrame containing SNP information, 
            with columns "Reference" and "Allele".

    Returns:
        tuple[int, int]: A tuple containing two integers: the number of A3A SNPs 
            and the total number of SNPs in the input DataFrame.
    '''

    SNPs = len(df)
    if SNPs == 0:
        return 0, 0
    A3A_SNPs_CT = len(df[df.Reference.isin(["C"]) & df.Allele.isin(["T"])])
    A3A_SNPs_GA = len(df[df.Reference.isin(["G"]) & df.Allele.isin(["A"])])

    A3A_SNPs = A3A_SNPs_CT + A3A_SNPs_GA

    logging.debug(f"SNPs: {SNPs}, A3A_SNPs: {A3A_SNPs}, A3A_SNPs_CT: {A3A_SNPs_CT}, A3A_SNPs_GA: {A3A_SNPs_GA}")

    return A3A_SNPs, SNPs

def countA3ASNPs_list(
        df_list: Union[
            list[ParentalTable],
            list[SampleTable],
            list[OverlapSwitchesTable]]
        ) -> None:
    '''Counts the number of A3A SNPs in a list of SampleTable objects and logs 
    the results. A3A context is not taken into account, just the C/T and G/A

    Args:
        df_list (list[SampleTable]): A list of SampleTable objects representing 
            the samples to be analyzed.

    Returns:
        None
    '''
    All_SNPs = 0
    All_A3A_SNPs = 0
    for sample in df_list:

        A3A_SNPs, SNPs = countA3ASNPs(sample.df)
        if SNPs == 0:
            logging.warning(f"{sample.name} has no SNPs... Moving on...")
            continue

        All_SNPs += SNPs
        All_A3A_SNPs += A3A_SNPs

        logging.info(f"{sample.name} has {(A3A_SNPs/SNPs)*100}% of A3A SNPs (out of {SNPs} total)")

    logging.info(f"out of {All_SNPs}SNPs, there are {(All_A3A_SNPs/All_SNPs)*100}% A3A SNPs")

def createRecombinationMaps(
        frames_list_samples: list[SampleTable],
        parental_origin_db_df: pd.DataFrame,
        df_chr_lengths: pd.DataFrame, 
        starts: pd.DataFrame,
        draw_map: bool = True,
        saveMap: bool = True,
        min_SNPs_switch_support: int = 1,
        snps_and_rec_ev_file: str = "outputs/parent_snps_and_rec_ev.txt"
        ) -> list[SwitchesTable]:

    '''Creates recombination maps for a list of SampleTable objects and returns 
    a list of SwitchesTable objects representing the switches in each sample.

    Args:
        frames_list_samples (list): A list of SampleTable objects representing 
            the samples to be analyzed.
        parental_origin_db_df (pandas.DataFrame): A pandas DataFrame containing 
            the parental origin information for each SNP.
        df_chr_lengths (pandas.DataFrame): A pandas DataFrame containing the 
            lengths of each chromosome.
        starts (pandas.DataFrame): A pandas DataFrame containing the start 
            positions of each chromosome. Needed in case no SNPs are found in a
            chromosome.
        draw_map (bool, optional): A boolean indicating whether to draw a SNP 
            map for each sample. Defaults to True.
        saveMap (bool, optional): A boolean indicating whether to save the map 
            as a PNG file. Defaults to True.
        min_SNPs_switch_support (int, optional): An integer indicating the 
            minimum number of SNPs required to support a switch. Defaults to 1.
        snps_and_rec_ev_file (str, optional): A string representing the path to 
            the output file for SNP and recombination event information. 
            Defaults to "outputs/parent_snps_and_rec_ev.txt".

    Returns:
        list[SwitchesTable]: A list of SwitchesTable objects representing the 
            switches in each sample.
    '''

    # make figures directory if it doesn't exist
    if not os.path.exists("figures"):
        os.makedirs("figures")
        logging.info("Created figures directory")

    if os.path.exists(snps_and_rec_ev_file):
        os.remove(snps_and_rec_ev_file)
        logging.info(f"Removed old {snps_and_rec_ev_file} file")

    with open(snps_and_rec_ev_file, "a") as f:

        switches_list = []

        for sample in frames_list_samples:
            sample.regionAsInt()
            sample.innerMergeWithParentalOrigin(
                parental_origin_db_df, 
                min_SNPs_switch_support=min_SNPs_switch_support
                )

            switchesTable_object = SwitchesTable(
                df=sample.switches_df,
                name=sample.name,
                percent=sample.percent,
                aneuploid = sample.aneuploid,
                heterozygous = sample.heterozygous,
                aneuploid_chromosomes = sample.aneuploid_chromosomes,
                heterozygous_chromosomes = sample.heterozygous_chromosomes,
                usable = sample.usable,
                )
            switches_list.append(switchesTable_object)

            f.write(sample.output_str)

            if draw_map:
                drawSNPMap(sample, df_chr_lengths, starts, saveMap=saveMap, map_type="switches") #SNPs


    return switches_list

@time_elapsed
def createDatasetNullHypothesis( switches_list: list[SwitchesTable], chromosome_lengths: list[list]) -> None:
    '''Simulates random switches for each SwitchesTable object in a list, 
    creating a null hypothesis dataset.

    Args:
        switches_list (list[SwitchesTable]): A list of SwitchesTable objects 
            representing the switches in each sample.
        df_chr_lengths (pandas.DataFrame): A pandas DataFrame containing the 
            lengths of each chromosome.

    Returns:
        None
    '''
    for i in switches_list:
        i.simulateRandomSwitches(chromosome_lengths)

@fancy_status
@time_elapsed
def addSwitchesAndClustersToNormalSamples(
        switches_overlap_samples: list[OverlapSwitchesTable], 
        switches_list: list,
        df_clusters_type: Literal["PMACD", "JT"] = "JT",
        cluster_significance = 0.001
        ) -> None:
    '''Adds switches and clustered SNPs to a list of OverlapSwitchesTable objects 
    representing normal samples.

    Args:
        switches_overlap_samples (list): A list of OverlapSwitchesTable objects 
            representing the normal samples to be analyzed.
        switches_list (list): A list of SwitchesTable objects representing the 
            switches in each sample.
        df_clusters (pandas.DataFrame): A pandas DataFrame containing the 
            clustered SNP information.

    Returns:
        None
    '''
    for index, sample in enumerate(switches_overlap_samples):

        sample.addSwitchesTable(switches_list[index])
        sample.regionAsInt()

        sample.findClusteredSNPs(df_clusters_type=df_clusters_type, cluster_significance=cluster_significance)

        #remove the unclustered SNPs
        sample.df_SNPs = sample.df_SNPs[sample.df_SNPs["Is_Clustered"] == True]

        logging.info(f"Saved {len(sample.df_SNPs)} clustered SNVs for {sample.name}... into sample.df_SNPs")
        logging.debug(f"sample.df_SNPs: {sample.df_SNPs}")

@time_elapsed
def findAssociatedClusters(
        frames_list_samples_normal_ref: list[OverlapSwitchesTable],
        window_size: int,
        generation: int = 0,
        exclude_low_SNP_density_regions: bool = True
        ) -> pd.DataFrame:
    '''Finds clusters of SNPs within a given window size of each switch in a 
    list of SampleTable objects representing normal samples.

    Args:
        frames_list_samples_normal_ref (list[OverlapSwitchesTable]): A list of SampleTable 
            objects representing the normal samples to be analyzed.
        window_size (int): An integer representing the size of the window around 
            each switch to search for clusters in.
        generation (int, optional): An integer representing the generation of 
            simulated data to use for comparison. Defaults to 0. 0 is different
            because we want to find the association in the original data and
            then compare between that and the simulated data.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing information about the 
            clusters found in each sample, including the sample name, chromosome, 
            cluster ID, and the number of SNPs in the cluster.
    '''
    logging.debug(f"Looking for clusers within {window_size}bp of the switch...")

    if exclude_low_SNP_density_regions:
        #still needs to be established if max_distance should be fixed or dynamic with window_size...
        lowDensityMap = loadDensityMap(lowDensity=True, max_distance=1000)
    else:
        #needed to bind the variable
        lowDensityMap = dict()

    counted_clusters_df_list = []
    for i, sample in enumerate(frames_list_samples_normal_ref):
        logging.debug(f"The i in findAssociatedClusters is: {i}")

        if generation == 0:
            logging.debug(f"Finding associated clusters in {sample.name}...")
            sample.findClusters(window_size)

        logging.debug(f"Finding random associated clusters in {sample.name}... (generation {generation})")
        sample.findRandomClusters(window_size)

        if exclude_low_SNP_density_regions:
            sample.markLowSNPDensityClusters(lowDensityMap, window_size=500)
            counted_clusters_df = sample.df_SNPs[sample.df_SNPs["Low_Density"] == False].copy()
        else:
            counted_clusters_df = sample.df_SNPs.copy()

        counted_clusters_df = counted_clusters_df[[
            "Cluster_ID", 
            "Clustered_Proximal", 
            "Random_Clustered_Proximal",
            "Chromosome",
            "Region"]].groupby([
                "Cluster_ID", 
                "Clustered_Proximal", 
                "Random_Clustered_Proximal",
                "Chromosome"], as_index=False).agg({'Region': ['min', 'max', 'count']})
        counted_clusters_df.columns = ["_".join(a) for a in counted_clusters_df.columns.to_flat_index()]

        counted_clusters_df["Sample_Name"] = sample.name
        counted_clusters_df["Percent_Parental_SNPs"] = sample.percent

        counted_clusters_df_list.append(counted_clusters_df)

    #concat the dataframes
    master_counted_clusters_df = pd.concat(counted_clusters_df_list, ignore_index=True)

    return(master_counted_clusters_df)

def overlapClustersWithSwitches(
        samples_SNVs: list[OverlapSwitchesTable],
        window_size: int,
        switches_list: list,
        df_chr: pd.DataFrame,
        drawMap: bool = True,
        df_clusters_type: Literal["JT", "PMACD"] = "JT",
        cluster_significance: Union[float, str] = 0.001
        ) -> pd.DataFrame:
    '''Finds clusters of SNPs within a given window size of each switch in a 
    list of SampleTable objects representing normal samples, and returns a 
    pandas DataFrame containing information about the clusters found in each 
    sample.

    Args:
        samples_SNVs (list[OverlapSwitchesTable]): A list of OverlapSwitchesTable 
            objects representing the normal samples to be analyzed.
        window_size (int): An integer representing the size of the window around 
            each switch to search for clusters in.
        switches_list (list[SwitchesTable]): A list of SwitchesTable objects 
            representing the switches in each sample.
        df_chr (pandas.DataFrame): A pandas DataFrame containing the lengths and 
            start positions of each chromosome.
        drawMap (bool, optional): A boolean indicating whether to draw a SNP map 
            for each sample. Defaults to True.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing information about the 
            clusters found in each sample, including the sample name, chromosome, 
            cluster ID, and the number of SNPs in the cluster.
    '''

    logging.info(f"Looking for clusers within {window_size}bp of the switch...")

    master_counted_clusters_df = pd.DataFrame()
    cluster_id_index = 0

    for i, sample in enumerate(samples_SNVs):

        sample.addSwitchesTable(switches_list[i])
        sample.addRandomSwitchesTable(switches_list[i])

        sample.findClusteredSNPs(df_clusters_type=df_clusters_type, cluster_significance=cluster_significance)

        sample.findClusters(window_size=window_size)
        sample.findRandomClusters(window_size=window_size)
        sample.df_SNPs = sample.df_SNPs[sample.df_SNPs["Is_Clustered"] == True]

        print(f"{sample.name}, {len(sample.df_SNPs)} clustered SNPs, heterozygous: {sample.heterozygous}, aneuploid: {sample.aneuploid}")

        if df_clusters_type == "JT":

            sample.find_GC_Clusters(window_size=window_size)
            sample.find_CO_Clusters(window_size=window_size)

            #fill NaNs with Fallse for GC_Proximal and CO_Proximal columns
            sample.df_SNPs["GC_Proximal"].fillna(False, inplace=True)
            sample.df_SNPs["CO_Proximal"].fillna(False, inplace=True)

            counted_clusters_df = sample.df_SNPs[[
                "Cluster_ID", 
                "Clustered_Proximal", 
                "Random_Clustered_Proximal",
                "Chromosome",
                "GC_Proximal",
                "CO_Proximal",
                "Region"]].groupby([
                    "Cluster_ID", 
                    "Clustered_Proximal",
                    "GC_Proximal",
                    "CO_Proximal", 
                    "Random_Clustered_Proximal", 
                    "Chromosome"], as_index=False).agg(
                        {'Region': ['min', 'max', 'count']})
            
            counted_clusters_df.columns = ["_".join(a) for a in counted_clusters_df.columns.to_flat_index()]
            all_clusters = len(counted_clusters_df["Cluster_ID_"].unique())

            if sample.aneuploid == True or sample.heterozygous == True:
                logging.warning(f"Cannot use complete GC & CO analysis for sample {sample.name}. Removing aneuploid chromosomes: {sample.aneuploid_chromosomes}, and heterozygous chromosomes {sample.heterozygous_chromosomes} from analysis.")

            for chromosome in sample.heterozygous_chromosomes:
                #remove the chromosome from the counted_clusters_df
                counted_clusters_df = counted_clusters_df[counted_clusters_df["Chromosome_"] != chromosome]

            for chromosome in sample.aneuploid_chromosomes:
                #remove the chromosome from the counted_clusters_df
                counted_clusters_df = counted_clusters_df[counted_clusters_df["Chromosome_"] != chromosome]
        

        elif df_clusters_type == "PMACD":
            logging.warning(f"df_clusters_type {df_clusters_type} is not supported for GC & CO analysis. Skipping GC and CO clusters.")

            counted_clusters_df = sample.df_SNPs[[
                "Cluster_ID", 
                "Clustered_Proximal", 
                "Random_Clustered_Proximal",
                "Chromosome",
                "Region"]].groupby([
                    "Cluster_ID", 
                    "Clustered_Proximal", 
                    "Random_Clustered_Proximal", 
                    "Chromosome"], as_index=False).agg(
                        {'Region': ['min', 'max', 'count']})
            
            counted_clusters_df.columns = ["_".join(a) for a in counted_clusters_df.columns.to_flat_index()]
            all_clusters = len(counted_clusters_df["Cluster_ID_"].unique())
                               

        counted_clusters_df["Sample_Name"] = sample.name
        counted_clusters_df["Percent_Parental_SNPs"] = sample.percent
        #use df.index to make unique Dataset_Cluster_ID
        counted_clusters_df["Dataset_Cluster_ID"] = counted_clusters_df.index + 1 + cluster_id_index
        cluster_id_index = cluster_id_index + all_clusters

        #sample.df_SNPs = pd.concat([sample.df_SNPs, starts], ignore_index=True).sort_values(by=["Chromosome", "Region"], ascending=(True, True), ignore_index=True)
        #sample.df_SNPs = sample.df_SNPs.sort_values(by=["Chromosome", "Region"], ascending=(True, True), ignore_index=True)

        if drawMap:
            drawSNPMap(sample, df_chr, starts, saveMap=True, showMap=False, map_type=f"clustersSwitchesOverlaps_{window_size}bp")

        master_counted_clusters_df = pd.concat([master_counted_clusters_df, counted_clusters_df], ignore_index=True)

    #make sure that Dataset_Cluster_ID is unique, if not print a warning
    if len(master_counted_clusters_df["Dataset_Cluster_ID"].unique()) != len(master_counted_clusters_df):
        logging.warning(f"Dataset_Cluster_ID is not unique in the master_counted_clusters_df. There are {len(master_counted_clusters_df)} rows in the dataframe, but only {len(master_counted_clusters_df['Dataset_Cluster_ID'].unique())} unique values in the Dataset_Cluster_ID column.")


    #drop the _Cluster_ID column and rename Dataset_Cluster_ID to _Cluster_ID
    master_counted_clusters_df = master_counted_clusters_df.drop(columns=["Cluster_ID_"])
    master_counted_clusters_df = master_counted_clusters_df.rename(columns={"Dataset_Cluster_ID": "Cluster_ID_"})

    return(master_counted_clusters_df)

def saveSwitchesToFile(
        switches_list: list,
        path: str = "outputs/all_switches.tsv"
        ) -> None:
    '''Saves information about switches in a list of SwitchesTable objects to 
    a TSV file.

    Args:
        switches_list (list[SwitchesTable]): A list of SwitchesTable objects 
            representing the switches in each sample.
        path (str, optional): A string representing the path to the output file. 
            Defaults to "outputs/all_switches.tsv".

    Returns:
        None
    '''

    all_switches_df = pd.DataFrame()
    for i in switches_list:
        dfi = i.df.copy()
        dfi["Sample_Name"] = i.name
        dfi["Percent_Parental"] = i.percent
        all_switches_df = pd.concat([all_switches_df, dfi], ignore_index=True)

    all_switches_df.to_csv(path, sep='\t', index = False)
    logging.info(f"Saved all switches to {path}")

@time_elapsed
def conductNullHypothesisSimulation(
        switches_list: list, 
        chromosomes, 
        frames_list_samples_normal_ref: list, 
        simulation_generations: int = 10,
        window_params: tuple[int, int, int] = (0, 50001, 500),
        save_file_name: str = "null_hypothesis_simulation_na.csv"
        ) -> pd.DataFrame:
    '''Conducts a null hypothesis simulation by randomly changing position of
    switches in a list of SampleTable objects representing normal samples, 
    and counting the number of clusters of SNPs within a given window size of
    each switch.

    Args:
        switches_list (list): A list of SwitchesTable objects representing the 
            switches in each sample.
        chromosomes: The chromosomes to be analyzed.
        frames_list_samples_normal_ref (list): A list of SampleTable objects 
            representing the normal samples to be analyzed.
        simulation_generations (int, optional): The number of simulation generations 
            to run. Defaults to 10.
        window_params (tuple[int, int, int], optional): A tuple containing three
            integers representing the start, stop, and step values for the window
            size. Defaults to (0, 50001, 500).
        save_file_name (str, optional): The name of the file to save the results to. 
            Defaults to "null_hypothesis_simulation_na.csv".

    Returns:
        pandas.DataFrame: A pandas DataFrame containing information about the 
            clusters found in each sample, including the sample name, chromosome, 
            cluster ID, and the number of SNPs in the cluster.
    '''
        
    percent_parental_min = 70

    logging.info(f"Considering clusters with at least {percent_parental_min}% parental SNPs")
    #remove samples with less than 70% parental SNPs from
    switches_list[:] = [i for i in switches_list if i.percent >= percent_parental_min]
    frames_list_samples_normal_ref[:] = [i for i in frames_list_samples_normal_ref if i.percent >= percent_parental_min]
    logging.debug(f"The samples have been filtered by % parental SNPs. The following samples are left: {[i.name for i in switches_list]}")

    proximal_clusters = pd.DataFrame()

    for generation in range(simulation_generations):
        logging.info(f"###$$$### GENERATION:{generation:<5} ###$$$###")

        createDatasetNullHypothesis(switches_list, chromosomes)
        for i, frame in enumerate(frames_list_samples_normal_ref):
            frames_list_samples_normal_ref[i].addRandomSwitchesTable(switches_list[i])

        window_proximal_clusters = pd.DataFrame()
        window_start, window_stop, window_step = window_params
        
        for window_size in range(window_start, window_stop, window_step): #(0, 50001, 500) is default
            logging.debug(f"GEN: {generation:<5}, WINDOW SIZE: {window_size:<5}")

            master_counted_clusters_df = findAssociatedClusters(frames_list_samples_normal_ref, window_size, generation=generation, exclude_low_SNP_density_regions=False)

            master_len = len(master_counted_clusters_df)
            random_filtered_master_len = master_counted_clusters_df["Random_Clustered_Proximal_"].sum()
            filtered_master_len = master_counted_clusters_df["Clustered_Proximal_"].sum()
            
            concat_dict = {
                'Generation' : generation,
                'Window_Size' : window_size,
                'Associated_Clusters' : filtered_master_len,
                'Other_Clusters' : (master_len - filtered_master_len),
                'Random_Associated_Clusters' : random_filtered_master_len,
                'Random_Other_Clusters' : (master_len - random_filtered_master_len)}

            #make the dictionary into a dataframe
            concat_pd = pd.DataFrame(concat_dict, index=[0])
            window_proximal_clusters = pd.concat([window_proximal_clusters, concat_pd], ignore_index=True)
        
        #if its the first generation, save the file (overwrite), else append
        if generation == 0:
            window_proximal_clusters.to_csv(save_file_name, index=False)
        else:
            window_proximal_clusters.to_csv(save_file_name, mode='a', header=False, index=False)

        logging.info(f"Saved Generation:{generation} to file: {save_file_name}")

        proximal_clusters = pd.concat([proximal_clusters, window_proximal_clusters], ignore_index=True)
            
    return proximal_clusters

def drawContextA3A(
        sample_list: list[OverlapSwitchesTable], 
        show: bool = False, 
        save: bool = False,
        focus_set: list | None = None,
        save_name_suffix: str = ""
        ) -> None:
    '''Draws an A3A context plot for combined samples SNVs. Requires the
    BioAid package.

    Args:
        sample_list (list): A list of SampleTable objects 
            representing the normal samples to be analyzed.
        show (bool, optional): A boolean indicating whether to show the plot. 
            Defaults to False.
        save (bool, optional): A boolean indicating whether to save the plot. 
            Defaults to False.

    Returns:
        None
    '''

    from BioAid import pullGenomicContext, drawGenomicContext, rev_compl

    context_list_all_C = []
    context_list_all_G = []

    for sample in sample_list:

        if focus_set is not None:
            if sample.name not in focus_set:
                continue

        df = sample.df.copy()
        sample_name = sample.name

        logging.info(f"Analyzing {sample_name}...")

        df['Region'] = df['Region'].astype(int)

        df_C = df[df["Reference"] == 'C']
        positions_C = df_C[['Chromosome','Region']].values.tolist()
        context_list_C = pullGenomicContext(positions_C, 'S288C_reference_sequencewoUra3wUra329wNAT.fa', context_flank=4)

        df_G = df[df["Reference"] == 'G']        
        positions_G = df_G[['Chromosome','Region']].values.tolist()
        context_list_G = pullGenomicContext(positions_G, 'S288C_reference_sequencewoUra3wUra329wNAT.fa', context_flank=4)
        
        try:
            if save:
                drawGenomicContext(context_list_C, save_path=f'figures/contextFigs/context_C_{sample_name}.png')
            else:
                drawGenomicContext(context_list_C, show=show)

            context_list_all_C.extend(context_list_C)
        
        except Exception as e:
            logging.warning(f"Could not draw context plot for {sample_name}. Skipping...")
            logging.warning(e)
        
        try:
            if save:
                drawGenomicContext(context_list_G, save_path=f'figures/contextFigs/context_G_{sample_name}.png')
            else:
                drawGenomicContext(context_list_G, show=show)

            context_list_all_G.extend(context_list_G)
        
        except Exception as e:
            logging.warning(f"Could not draw context plot for {sample_name}. Skipping...")
            logging.warning(e)
        

    #get a reverse complement of the items in the list, also reverse the order of the list items
    context_list_all_G_rev_compl = []
    for left, center, right, chr in context_list_all_G.copy():
        length_context = len(left)
        seq_str = left + center + right
        assert (len(seq_str) == ((length_context*2) + 1)), "The length of the sequence string is not correct."
        seq_str_rev_compl = rev_compl(seq_str)
        new_left = seq_str_rev_compl[:length_context]
        new_center = seq_str_rev_compl[length_context]
        new_right = seq_str_rev_compl[length_context+1:]
        context_list_all_G_rev_compl.append([new_left, new_center, new_right, chr])

    context_list_all_both = context_list_all_C + context_list_all_G_rev_compl

    if save:
        drawGenomicContext(context_list_all_C, save_path=f'figures/contextFigs/context_C_all{save_name_suffix}.png', show=show)
        drawGenomicContext(context_list_all_G, save_path=f'figures/contextFigs/context_G_all{save_name_suffix}.png', show=show)
        drawGenomicContext(context_list_all_both, save_path=f'figures/contextFigs/context_both_all{save_name_suffix}.png', show=show)
    else:
        drawGenomicContext(context_list_all_C, show=show)
        drawGenomicContext(context_list_all_G, show=show)
        drawGenomicContext(context_list_all_both, show=show)

def removeA3AmotifSNPs(
        df: pd.DataFrame, 
        reference_genome_path: str = 'S288C_reference_sequencewoUra3wUra329wNAT.fa'
        ) -> pd.DataFrame:
    '''Removes SNPs with A3A motif from a pandas DataFrame containing SNP information.
    Requires BioAid.

    Args:
        df (pandas.DataFrame): A pandas DataFrame containing SNP information, 
            including the chromosome, region, reference allele, and alternate allele.
        reference_genome_path (str, optional): A string representing the path to 
            the reference genome file. Defaults to 'S288C_reference_sequencewoUra3wUra329wNAT.fa'.

    Returns:
        pandas.DataFrame: A pandas DataFrame containing SNP information with A3A 
            motif SNPs removed.
    '''
    from BioAid import pullGenomicContext

    df_C = df[df["Reference"] == 'C']
    df_C = df_C[df_C["Allele"] == 'T']
    positions_C = df_C[['Chromosome','Region']].values.tolist()
    context_list_C = pullGenomicContext(positions_C, reference_genome_path, context_flank=1)

    df_G = df[df["Reference"] == 'G']
    df_G = df_G[df_G["Allele"] == 'A']
    positions_G = df_G[['Chromosome','Region']].values.tolist()
    context_list_G = pullGenomicContext(positions_G, reference_genome_path, context_flank=1)
    bad_positions = []

    for i in context_list_C:
        if (i[0] == 'T') & (i[2] != 'C'):
            bad_positions.append([i[3], i[4]])

    for i in context_list_G:
        if (i[0] != 'G') & (i[2] == 'A'):
            bad_positions.append([i[3], i[4]])

    logging.info(f"Found {len(bad_positions)} A3A motif SNPs... Removing them from parental origin database")
    df = df[~df.set_index(['Chromosome','Region']).index.isin(bad_positions)]
    logging.info(f"Removed {len(bad_positions)} A3A motif SNPs from parental origin database")
    logging.info(f"Parental origin database now has {len(df)/2} SNPs positions")
    
    return(df)

def correlateSNVsWithSwitches(frames_list_samples_normal_ref: list) -> None:
    '''Correlates SNVs and switches from a list of SampleTable objects representing 
    normal samples, and plots their distribution on each chromosome as a histogram.

    Args:
        frames_list_samples_normal_ref (list): A list of SampleTable objects 
            representing the normal samples to be analyzed.

    Returns:
        None
    '''

    sample_names = [
        'JK10','JK11','JK13', 'JK14','JK16',
        'JK17','JK19','JK21','JK22','JK23',
        'JK24','JK26','JK27','JK29','JK30',
        'JK31','JK32','JK33','JK35','JK5',
        'JK58','JK59','JK6','JK62','JK64',
        'JK65','JK66','JK67', 'JK68',
        'JK7','JK8','JK81','JK84','JK9',
        ]

    #create a dataframe with all SNVs
    df_SNPs_all = pd.DataFrame()
    for i in frames_list_samples_normal_ref:
        
        #only take samples from the list
        if i.name.strip("Mal_") not in sample_names:
            print(f"{i.name} not in sample_names list")
            continue

        df_temp = i.df.copy()
        df_temp["sample_name"] = i.name
        df_temp["Event"] = "SNV"
        print(i)
        df_temp.rename(columns={"Region": "Event_Position"}, inplace=True)
        df_SNPs_all = pd.concat([df_SNPs_all, df_temp], ignore_index=True)
    #print out total number of SNVs
    print("Total number of SNVs: ", len(df_SNPs_all))

    #create a dataframe with all switches
    df_switches_all = pd.DataFrame()
    for i in frames_list_samples_normal_ref:

        #only take samples from the list
        if i.name.strip("Mal_") not in sample_names:
            continue

        df_temp = i.switchesTable.copy()
        df_temp["sample_name"] = i.name
        df_temp["Event"] = "Switch"
        print(i)
        df_temp.rename(columns={"Switch_Center": "Event_Position"}, inplace=True)
        df_switches_all = pd.concat([df_switches_all, df_temp], ignore_index=True)

    #print out total number of switches
    print("Total number of switches: ", len(df_switches_all))

    #merge the two dataframes
    df_combined = pd.concat([df_SNPs_all, df_switches_all], ignore_index=True)
    
    #rename the Chromosomes to chr1, chr2, etc.
    df_combined = renameChrToRoman(df_combined, "Chromosome")
    print(df_combined.head())


    #plot the distribution of SNVs and switches on each chromosome as a histogram (x-axis: chromosome position, y-axis: number of events)
    for i in df_combined["Chromosome"].unique():
        
        #make the histogram of SNVs and switches, on the same plot, with different colors and separate y-axes
        plt.rcParams['figure.figsize'] = [30, 5]
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.hist(df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]=="SNV")]["Event_Position"], bins=300, color="blue", alpha=0.5)
        ax2.hist(df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]=="Switch")]["Event_Position"], bins=300, color="red", alpha=0.5)

        #set the y-axis limits
        ax1.set_ylim(0, 50)
        ax2.set_ylim(0, 50)

        ax1.set_ylabel("SNVs", color="blue")
        ax2.set_ylabel("Switches", color="red")
        ax1.set_xlabel("Chromosome position")
        ax1.set_title(i)
        plt.show()

def drawNullHypothesisSimulationResults(
        proximal_clusters_path: str, 
        show: bool = True, 
        x_ax: str = "window_size_kb",
        as_percent:bool = False, 
        other_clusters: bool = False, 
        bootstrap: bool = False,
        save_path = None
        ) -> None:
    """
    Draw a line plot of the number of clusters associated with a recombination hotspot.

    Args:
        proximal_clusters_path (str): Path to the CSV file containing the proximal clusters data.
        show (bool, optional): Whether to show the plot. Defaults to True.
        x_ax (str, optional): The x-axis variable. Defaults to "Window_Size_kb".
        as_percent (bool, optional): Whether to plot the y-axis as a percentage. Defaults to False.
        other_clusters (bool, optional): Whether to plot the number of other clusters. Defaults to False.
        bootstrap (bool, optional): Whether to use bootstrap resampling. Defaults to False.
        save_path (str, optional): Path to save the plot. Defaults to None.
    """
    clusters_df = pd.read_csv(proximal_clusters_path)

    clusters_df["total_clusters"] = clusters_df["Associated_Clusters"] + clusters_df["Other_Clusters"]
    clusters_df["associated_clusters_percent"] = clusters_df["Associated_Clusters"] / clusters_df["total_clusters"] * 100
    clusters_df["other_clusters_percent"] = clusters_df["Other_Clusters"] / clusters_df["total_clusters"] * 100
    clusters_df["random_associated_clusters_percent"] = clusters_df["Random_Associated_Clusters"] / clusters_df["total_clusters"] * 100
    clusters_df["random_other_clusters_percent"] = clusters_df["Random_Other_Clusters"] / clusters_df["total_clusters"] * 100
    clusters_df["window_size_kb"] = clusters_df["Window_Size"] / 1000

    actual_df = clusters_df[clusters_df['Generation'] == 0]

    y_ax_actual = "Associated_Clusters"
    y_ax_other = "Other_Clusters"
    y_ax_random = "Random_Associated_Clusters"
    y_ax_random_other = "Random_Other_Clusters"
    ylabel = "Number of clusters"
    xlabel = "Window size (kb)"
    legend_labels = ["Actual", "Random"]

    if as_percent:
        ylabel = "Percent of clusters"
        y_ax_actual = "associated_clusters_percent"
        y_ax_other = "other_clusters_percent"
        y_ax_random = "random_associated_clusters_percent"
        y_ax_random_other = "random_other_clusters_percent"


    sns.set_context(context="poster", font_scale=2, rc={"lines.linewidth": 8})
    sns.set_style("whitegrid")
    plt.figure(figsize=(16, 16))

    if bootstrap:
        #legend_labels = ["Actual", "Random (bootstrap)"]
        sns.lineplot(data=actual_df, x=x_ax, y=y_ax_actual, label="Actual", color="tab:blue")
        #sns.lineplot(data=clusters_df, x=x_ax, y=y_ax_random, palette="flare", label="Random (bootstrap)", linestyle="--")
        sns.lineplot(data=clusters_df, x=x_ax, y=y_ax_random, errorbar=("pi", 95), label="Random", err_style="band", linestyle="--", color="tab:green")
    else:
        sns.lineplot(data=actual_df, x=x_ax, y=y_ax_actual, label="Actual", color="tab:blue")
        sns.lineplot(data=clusters_df, x=x_ax, y=y_ax_random, palette="flare", hue="Generation", legend="", alpha=0.05, linewidth=3)

    if other_clusters:
        sns.lineplot(data=actual_df, x=x_ax, y=y_ax_other, label="Actual")
        if bootstrap:
            sns.lineplot(data=clusters_df, x=x_ax, y=y_ax_random_other, palette="flare", label="Random (bootstrap)", legend="")
        else:
            sns.lineplot(data=clusters_df, x=x_ax, y=y_ax_random_other, palette="flare", label="Random", hue="Generation", legend="")

    plt.ylim(bottom=0)
    plt.title("Proportion of recombination associated clusters", fontweight='bold', fontsize=46)
    plt.xlabel(xlabel, fontweight='bold', fontsize=34)
    plt.ylabel(ylabel, fontweight='bold', fontsize=34)
    #increase the size of the x and y axis ticks
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.legend(fontsize=28)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def plot_genomic_SNP_coverage(
        df: pd.DataFrame,
        title: str,
        max_coverage_show: int|float = 5, 
        save: bool = False, 
        show: bool = True
        ) -> None:
    '''
    Plots a scatter plot of the genomic SNP coverage.

    Args:
        df (pd.DataFrame): A pandas DataFrame containing the SNP data.
        title (str): A string representing the title of the plot.
        max_coverage_show (int|float): An integer or float representing the maximum coverage to show on the plot.
            Default value is 5.
        save (bool): A boolean representing whether to save the plot to a file.
            Default value is False.
        show (bool): A boolean representing whether to show the plot.
            Default value is True.

    Returns:
        None
    '''
    from mayo.settings import chromosomes as chromosomes
    df_chr = pd.DataFrame(chromosomes, columns=["Chromosome", "end_position"])

    position_column = "Region"

    #we need to linearize all chromosomes onto one axis (column), named position_x_axis
    #we will use the end position of each chromosome to linearize it
    #we will use the start position of each chromosome to draw the chromosome line
    df_graph = df.copy()
    df_graph[position_column] = df_graph[position_column].astype(int)

    #for each chromosomal position we need to add the length of chromosomes that are before it
    df_graph["position_x_axis"] = df_graph[position_column]
    for i in range(len(df_chr)):
        df_graph.loc[df_graph["Chromosome"] == df_chr["Chromosome"][i], "position_x_axis"] += df_chr["end_position"].shift(1).fillna(0).cumsum()[i]

    logging.debug(df_graph)

    #if coverage_norm is > than max_coverage_show, set it to max_coverage_show
    df_graph["Coverage_norm"] = df_graph["Coverage_norm"].apply(lambda x: max_coverage_show if x > max_coverage_show else x)

    #create a scatter plot of the SNPs
    plt.figure(figsize=(40, 5))
    #marker pallete for zygosity_from_freq column Homozygous = 0, Heterozygous = 1, Aneuploid = 2
    markers = {
        "Homozygous": "o", 
        "Heterozygous": "X", 
        "Unknown": "D"}

    sns.scatterplot(data=renameChrToRoman(df_graph, "Chromosome"), x="position_x_axis", y="Coverage_norm", hue="Chromosome", palette="tab20", s=12, alpha=0.75, edgecolor="none", style="zygosity_from_freq", markers=markers)

    #draw vertical lines for each chromosome end position, must be cumulative sum of chromosome lengths
    for i in range(len(df_chr)):
        plt.axvline(x=df_chr["end_position"].shift(1).fillna(0).cumsum()[i], color="grey", linestyle=(0, (1, 1)))
    
    #draw horizontal lines for the coverage at 1x, 2x and 3x
    plt.axhline(y=1, color="black", linestyle="--")
    plt.axhline(y=2, color="black", linestyle="--")
    plt.axhline(y=3, color="black", linestyle="--")

    plt.ylim(0, max_coverage_show)
    
    #set legend outside of plot to the right
    plt.legend(bbox_to_anchor=(1.001, 1), loc=2)
    
    #set title and axis labels
    plt.title(title, fontsize=20)
    plt.xlabel("Genomic position", fontsize=15)
    plt.ylabel("Coverage (normalized)", fontsize=15)

    #remove x axis ticks and labels
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    #save plot
    if save:

        #if the folder does not exist, create it
        if not os.path.exists("figures/coverage"):
            os.makedirs("figures/coverage")

        plt.savefig(f"figures/coverage/{title}_SNP_coverage_map.png", dpi=400, bbox_inches='tight')
    if show:
        plt.show()

    else:
        plt.close()

def removeTelomereSNPs(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Removes SNPs from the telomeres from a pandas DataFrame containing SNP information.
    '''
    from mayo.settings import telomeres as telomeres
    from mayo.settings import chromosomes as chromosomes

    df["Region"] = df["Region"].astype(int)

    len_df = len(df)
    for record in telomeres:
        chromosome_name, l_tel_len, r_tel_len = record
        #find chromosome length of the chromosome in the chromosome list
        chr_len = [x[1] for x in chromosomes if x[0] == chromosome_name][0]
        logging.debug(chromosome_name, chr_len, l_tel_len, r_tel_len)
        right_tel_start = chr_len - r_tel_len

        #remove SNPs from the left telomere
        df = df[~((df["Chromosome"] == chromosome_name) & (df["Region"] <= l_tel_len))]
        #remove SNPs from the right telomere
        df = df[~((df["Chromosome"] == chromosome_name) & (df["Region"] >= right_tel_start))]

    logging.info(f"Removed {len_df - len(df)} telomeric SNPs from dataframe")
    return(df)

def removeRepetitiveDNASNPs(parental_origin_db_df, repetitive_DNA_fasta_path):
    '''
    Removes SNPs from repetitive DNA from a pandas DataFrame containing SNP information.
    '''
    from mayo.settings import roman_to_S288C
    import BioAid as ba

    #open the repetitive DNA file with BioAid
    repetitive_DNA_list = ba.extractSeqFromFastaToList(repetitive_DNA_fasta_path)

    len_parental_origin_db_df = len(parental_origin_db_df)

    for record in repetitive_DNA_list:
        seq_name, seq = record
        #extract name from seq_name
        feature_name = seq_name.split(" ")[0].strip(">").strip()
        chromosome_name = seq_name.split(",")[1].split(" ")[2].strip()
        if chromosome_name not in roman_to_S288C.keys():
            logging.debug(f"Chromosome {chromosome_name} of {feature_name} is not in roman_to_S288C.keys(). This behavior is normal for Mitochondrial Sequences. Skipping...")
            continue
        chromosome_name = roman_to_S288C[chromosome_name]
        seq_coords = seq_name.split(",")[1].split(" ")[-1].strip()
        seq_start = int(seq_coords.split("-")[0])
        seq_end = int(seq_coords.split("-")[1])

        reporter_start_pos_chr3 = 211813
        offset_reporter = 2723
        if (chromosome_name == "ref|NC_001135|") & (seq_start >= reporter_start_pos_chr3):
            seq_start += offset_reporter
            seq_end += offset_reporter
            logging.debug(f"Added offset of ura3-29 reporter ({offset_reporter}) to {feature_name} on chromosome {chromosome_name} (region: {seq_start}-{seq_end})")

        #ura3_start_pos_chr5 = 116166
        ura3_end_pos_chr5 = 116970
        ura3_del_len = 805
        if (chromosome_name == "ref|NC_001137|") & (seq_start > ura3_end_pos_chr5):
            #remove the length of the ura3 deletion from the start and end positions
            seq_start -= ura3_del_len
            seq_end -= ura3_del_len
            logging.debug(f"Removed length of ura3 ({ura3_del_len}) from {feature_name} on chromosome {chromosome_name} (region: {seq_start}-{seq_end})")

        #sort the genomic coordinates, this is needed because the repetitive DNA file has some coordinates in reverse order (transcription direction)
        seq_genomic_coord = [seq_start, seq_end]
        seq_genomic_coord.sort()

        #remove SNPs that are between the seq_start and seq_end from the parental_origin_db_df for the corresponding chromosome
        parental_origin_db_df = parental_origin_db_df[~((parental_origin_db_df["Chromosome"] == chromosome_name) & (parental_origin_db_df["Region"] >= seq_genomic_coord[0]) & (parental_origin_db_df["Region"] <= seq_genomic_coord[1]))]
        logging.debug(f"Removed SNPs from {feature_name} on chromosome {chromosome_name} from dataframe (region: {seq_start}-{seq_end}), genomic coordinates: {seq_genomic_coord}")


    logging.info(f"Removed {len_parental_origin_db_df - len(parental_origin_db_df)} repetitive DNA SNPs from dataframe")
    return(parental_origin_db_df)

def createParentalOriginDF(
        parents_path,
        removeA3A_motif: bool = False, 
        remove_problematic_positions: bool =True,
        remove_telomeres: bool = True,
        remove_all_repetitive_DNA: bool = True,
        parental_origin_db_df_path_save: str = "outputs/parental_origin_db.tsv"
        ) -> pd.DataFrame:
    '''
    Creates a DataFrame with complete parental SNPs.

    Args:
        parents_path (str): A string representing the path to the parental SNP data.
        removeA3A_motif (bool): A boolean representing whether to remove SNPs with A3A motif.
            Default value is False.
        parental_origin_db_df_path_save (str): A string representing the path to save the parental origin database.
            Default value is "outputs/parental_origin_db.tsv".

    Returns:
        pd.DataFrame: A pandas DataFrame containing the complete parental SNPs.
    '''
    from mayo.settings import problematic_positions as problematic_positions
    from mayo.settings import config as cfg


    frames_list_parents = csvDataFrameImport(parents_path, kind="parental")
    subtracted_parent_list = removeCommonDataFrameRows(frames_list_parents[0], frames_list_parents[1])

    logging.info("Removing duplicated SNPs from parental SNPs (2 parents)")
    for i in range(len(subtracted_parent_list)):
        subtracted_parent_list[i].df = subtracted_parent_list[i].df[subtracted_parent_list[i].df.duplicated()==False]
        subtracted_parent_list[i].df = pd.merge(subtracted_parent_list[i].df, frames_list_parents[i].df, on=["Chromosome", "Region"], validate="1:1")

    logging.info("Performing quality filter on parental SNPs (2 parents)...")
    for i in subtracted_parent_list:
        i.qualitySNVFilter(stringent=True)

    parental_origin_db_df = createDataFrameWithCompleteParentalSNPs(subtracted_parent_list[0],subtracted_parent_list[1])
    logging.info(f"Unfiltered Parental origin database has {len(parental_origin_db_df)/2} SNPs positions")

    if remove_problematic_positions:
        problematic_positions_df = pd.DataFrame(data=problematic_positions, columns = ["Chromosome", "Region"])
        logging.debug(len(parental_origin_db_df))
        logging.debug(len(parental_origin_db_df[~parental_origin_db_df['Region'].isin(problematic_positions_df['Region'])]))
        logging.debug(len(parental_origin_db_df[~parental_origin_db_df.set_index(['Chromosome','Region']).index.isin(problematic_positions_df.set_index(['Chromosome','Region']).index)]))
        
        parental_origin_db_df = parental_origin_db_df[~parental_origin_db_df.set_index(['Chromosome','Region']).index.isin(problematic_positions_df.set_index(['Chromosome','Region']).index)]
        logging.info("Removed problematic positions from parental origin database")

    # Remove SNPs with A3A motif (if desired, however, this is not recommended as: (1) it will remove a lot of SNPs, (2) will remove SNPs that are not necessarily A3A-edited, (3) in tests, does not improve results)
    if removeA3A_motif:
        parental_origin_db_df = removeA3AmotifSNPs(parental_origin_db_df)

    if remove_telomeres:
        parental_origin_db_df = removeTelomereSNPs(parental_origin_db_df)

    if remove_all_repetitive_DNA:
        parental_origin_db_df = removeRepetitiveDNASNPs(parental_origin_db_df, repetitive_DNA_fasta_path=cfg["repetitive_regions_fasta_path"])

    logging.info(f"Parental origin database has {len(parental_origin_db_df)/2} SNPs positions")

    parental_origin_db_df.to_csv(parental_origin_db_df_path_save, sep='\t', index=False)
    logging.info(f"Saved parental origin database to {parental_origin_db_df_path_save}")

    return parental_origin_db_df

def performHeterozygosityAndPloidyAnalysis(
        samples_parental_SNPs: list[SampleTable], 
        show_details: bool = False, 
        plot_coverage: bool = True, 
        save_plot: bool = True, 
        save_high_coverage_SNP_positions: bool = True, 
        remove_aneuploid_from_high_coverage_SNP_positions: bool = True,
        parental_origin_db_df: pd.DataFrame | None = None
        ) -> tuple[
                list[str],
                list[str], 
                list[str], 
                list[str], 
                list[str]]:
    '''
    Performs heterozygosity and ploidy analysis on a list of SampleTable objects.

    Args:
        samples_parental_SNPs (List[SampleTable]): A list of SampleTable objects representing the parental SNP data.
        show_details (bool): A boolean representing whether to show detailed information about the analysis.
            Default value is False.
        plot_coverage (bool): A boolean representing whether to plot the coverage of the SNPs.
            Default value is True.
        save_plot (bool): A boolean representing whether to save the plot of the SNP coverage.
            Default value is True.
        save_high_coverage_SNP_positions (bool): A boolean representing whether to save the positions of high coverage SNPs.
            Default value is True.
        remove_aneuploid_from_high_coverage_SNP_positions (bool): A boolean representing whether to remove aneuploid SNPs from the high coverage SNP positions.
            Default value is True.
        parental_origin_db_df (pd.DataFrame): A pandas DataFrame containing the complete parental SNPs.
            Default value is None. It is used to remove non-SNP positions from the samples.

    Returns:
        Tuple[List[str], List[str], List[str], List[str], List[str]]: A tuple containing five lists of sample names:
            - only_haploid_sample_list: samples with only haploid chromosomes
            - only_heterozygous_sample_list: samples with only heterozygous chromosomes
            - only_aneuploid_sample_list: samples with only aneuploid chromosomes
            - both_heterozygous_and_aneuploid_sample_list: samples with both heterozygous and aneuploid chromosomes
            - chimeric_sample_list: samples with chimeric chromosomes
    '''
    only_haploid_sample_list = []
    only_heterozygous_sample_list = []
    only_aneuploid_sample_list = []
    both_heterozygous_and_aneuploid_sample_list = []
    chimeric_sample_list = [] #chimeric for aneuploid and heterozygous


    #before doing anything, we need to remove positions that are not SNP positions by referencing the parental origin databaseDF
    if parental_origin_db_df is not None:
        logging.info("Removing non-SNP positions from samples...")
        parental_origin_db_df = parental_origin_db_df.copy()
        parental_origin_db_df = parental_origin_db_df[["Chromosome", "Region", "Allele"]]
        for i in samples_parental_SNPs:
            #make sure that Region is an integer
            i.df["Region"] = i.df["Region"].astype(int)
            i.df = pd.merge(i.df, parental_origin_db_df, on=["Chromosome", "Region", "Allele"], validate="1:1")
            logging.info(f"Sample {i.name} has {len(i.df)} SNP positions")

    for i in samples_parental_SNPs:

        heterozygous_chr_list = i.findHeterozygousChromosomes()
        aneuploid_chr_list = i.findAnuploidChromosomesBySNPCoverage(
            show_details=show_details,
            plot_coverage=plot_coverage, 
            save_plot=save_plot, 
            save_high_coverage_SNP_positions=save_high_coverage_SNP_positions, 
            remove_aneuploid_from_high_coverage_SNP_positions=remove_aneuploid_from_high_coverage_SNP_positions)
        
        #if both lists are empty, then the sample is presumed to be haploid
        if not heterozygous_chr_list and not aneuploid_chr_list:
            logging.info(f"Sample {i.name} is haploid")
            only_haploid_sample_list.append(i.name)
        else:
            #check what chromosomes are heterozygous, and what chromosomes are anuploid, and what chromosomes are both
            only_heterozygous_chr_list = [i for i in heterozygous_chr_list if i not in aneuploid_chr_list]
            only_aneuploid_chr_list = [i for i in aneuploid_chr_list if i not in heterozygous_chr_list]
            both_heterozygous_and_aneuploid_chr_list = [i for i in heterozygous_chr_list if i in aneuploid_chr_list]
            if only_heterozygous_chr_list:
                logging.warning(f"Sample {i.name} is heterozygous for chromosomes {only_heterozygous_chr_list}")
            if only_aneuploid_chr_list:
                logging.warning(f"Sample {i.name} is aneuploid for chromosomes {only_aneuploid_chr_list}")
            if both_heterozygous_and_aneuploid_chr_list:
                logging.warning(f"Sample {i.name} is heterozygous and aneuploid for chromosomes {both_heterozygous_and_aneuploid_chr_list}")

            if heterozygous_chr_list and not aneuploid_chr_list and not both_heterozygous_and_aneuploid_chr_list:
                only_heterozygous_sample_list.append(i.name)
            elif aneuploid_chr_list and not heterozygous_chr_list and not both_heterozygous_and_aneuploid_chr_list:
                only_aneuploid_sample_list.append(i.name)
            elif both_heterozygous_and_aneuploid_chr_list and (set(heterozygous_chr_list) == set(aneuploid_chr_list)):
                both_heterozygous_and_aneuploid_sample_list.append(i.name)
            else:
                chimeric_sample_list.append(i.name)

    logging.info(f"Samples with only haploid chromosomes ({len(only_haploid_sample_list)}): {only_haploid_sample_list}")
    logging.info(f"Samples with only heterozygous chromosomes ({len(only_heterozygous_sample_list)}): {only_heterozygous_sample_list}")
    logging.info(f"Samples with only aneuploid chromosomes ({len(only_aneuploid_sample_list)}): {only_aneuploid_sample_list}")
    logging.info(f"Samples with both heterozygous and aneuploid chromosomes ({len(both_heterozygous_and_aneuploid_sample_list)}): {both_heterozygous_and_aneuploid_sample_list}")
    logging.info(f"Samples with chimeric chromosomes ({len(chimeric_sample_list)}): {chimeric_sample_list}")
    logging.info(f"Total number of samples: {len(samples_parental_SNPs)}, normal:{len(only_haploid_sample_list)}, other:{len(only_heterozygous_sample_list) + len(only_aneuploid_sample_list) + len(both_heterozygous_and_aneuploid_sample_list) + len(chimeric_sample_list)}")

    return only_haploid_sample_list, only_heterozygous_sample_list, only_aneuploid_sample_list, both_heterozygous_and_aneuploid_sample_list, chimeric_sample_list

def getHaplotypeBlockLengths(
        switches_list: list[SwitchesTable],
        filter_set: list[str],
        random: bool = False,
        inverse=False,
        excluded_chromosomes: list[str] | None = ["ref|NC_001224|"],
        ) -> pd.DataFrame:
    """
    Returns a DataFrame containing the lengths of haplotype blocks between switch centers.

    Args:
        switches_list: A list of Switches objects.
        filter_set: A set of sample names to include or exclude, depending on the value of `inverse`.
        random: If True, use the `random_switch_df` attribute of each Switches object instead of the `df` attribute (default: False).
        inverse: If True, exclude samples in `filter_set` instead of including them (default: False).

    Returns:
        A Distance sorted (ascending) pandas DataFrame containing the following columns:
        - "Chromosome": The chromosome where the haplotype block is located.
        - "Switch_Center": The position of the switch center.
        - "Distance": The distance between the current and previous switch centers.
        - "samples_name": The name of the sample where the haplotype block is located.
    """
    #import chromosome lengths
    from mayo.settings import chromosomes as chromosomes
    chr_df = pd.DataFrame(chromosomes, columns=["Chromosome", "Length"])

    master_df = pd.DataFrame()
    for i in switches_list:
        
        #default
        if not inverse:
            if i.name not in filter_set:
                continue
        else:
            if i.name in filter_set:
                continue
        
        if not random:
            df = i.df.copy()
        else:
            df = i.random_switch_df.copy()

        #iterate through each chr in i.df and find all distances bwtween Switch_Centers
        df["samples_name"] = i.name

        #remove chromosomes that are not considered, e.g. mitochondria (default: ref|NC_001224|)
        if excluded_chromosomes is not None:
            df = df[~df["Chromosome"].isin(excluded_chromosomes)]

        for chr in df["Chromosome"].unique():
            df_temp = df[df["Chromosome"] == chr].copy()

            # create a "switch center" at Distance = 0 and Distance = [Chromosome length]
            df_temp = pd.concat([df_temp, pd.DataFrame({"Chromosome": chr, "Switch_Center": 0, "samples_name": i.name}, index=[0])], ignore_index=True)
            df_temp = pd.concat([df_temp, pd.DataFrame({"Chromosome": chr, "Switch_Center": chr_df[chr_df["Chromosome"] == chr]["Length"].values[0], "samples_name": i.name}, index=[0])], ignore_index=True)
            df_temp.sort_values(by=["Switch_Center"], inplace=True)
            df_temp["Distance"] = df_temp["Switch_Center"].diff()
            df_temp.dropna(inplace=True)
            master_df = pd.concat([master_df, df_temp], ignore_index=True)

    master_df.sort_values(by=["Distance"], inplace=True)
    
    return master_df

@time_elapsed
def getFractionOfGapsAtDistance(df: pd.DataFrame, distance: int = 50000, increment: int = 100) -> pd.DataFrame:
    """
    Returns a DataFrame containing the cumulative fraction of gaps (haplotype lengths) at each 
    distance in the range [0, distance], with an increment of `increment`.

    Args:
        df: A pandas DataFrame containing the following columns:
            - "Chromosome": The chromosome where the haplotype block is located.
            - "Switch_Center": The position of the switch center.
            - "Distance": The distance between the current and previous switch centers.
            - "samples_name": The name of the sample where the haplotype block is located.
        distance: The maximum distance to consider (default: 50000).
        increment: The distance increment to use (default: 100).

    Returns:
        A pandas DataFrame containing the following columns:
        - "Distance": The distance from the previous switch center.
        - "Fraction": The fraction of gaps at that distance.

    Raises:
        None
    """
    total_gaps = len(df)
    df_gaps = pd.DataFrame()
    #find the number of gaps at each distance (add 1 to distance to include the last distance)
    for i in range(0, distance+1, increment):

        #print progress every 1000
        if i % 1000 == 0:
            logging.debug(f"Processing distance {i}")

        df_temp = df[df["Distance"] < i+increment].copy()
        distance = i
        fraction = len(df_temp)/total_gaps
        row = {"Distance": distance, "Fraction": fraction}
        df_gaps = pd.concat([df_gaps, pd.DataFrame(row, index=[0])], ignore_index=True)

    df_gaps = df_gaps[["Distance", "Fraction"]]

    return df_gaps

def getHaplotypeLengthDistributionVsExpected(
        switches_list: list[SwitchesTable],
        filter_set: list[str], 
        inverse: bool = False, 
        random_trials: int = 100,
        fraction_increment: int = 100
        ) -> pd.DataFrame:
    """
    Returns a DataFrame containing the fraction of haplotype blocks at each distance from the previous switch center.

    Args:
        switches_list: A list of Switches objects.
        filter_set: A set of sample names to include or exclude, depending on the value of `inverse`.
        inverse: If True, exclude samples in `filter_set` instead of including them (default: False).

    Returns:
        A pandas DataFrame containing the following columns:
        - "Distance_kb": The distance from the previous switch center, in kilobases.
        - "Fraction_real": The fraction of haplotype blocks at that distance in the real data.
        - "Fraction_random": The fraction of haplotype blocks at that distance in the random data.
    """
    #get the haplotype block lengths for real data
    master_df = getHaplotypeBlockLengths(switches_list, filter_set, random=False ,inverse=inverse)
    logging.info(f"Total number of gaps: {len(master_df)}")
    df_gaps = getFractionOfGapsAtDistance(master_df, increment=fraction_increment)
    df_gaps["Fraction_real"] = df_gaps["Fraction"]

    #here might need to create multiple random datasets and average them...
    from mayo.settings import chromosomes as chromosomes
    df_gaps_rand_master = pd.DataFrame()
    for i in range(random_trials):
        createDatasetNullHypothesis(switches_list, chromosomes)
        #get the haplotype block lengths for random data
        master_df_rand = getHaplotypeBlockLengths(switches_list, filter_set, random=True ,inverse=inverse)
        logging.info(f"Total number of random gaps trial{i}: {len(master_df_rand)}")
        df_gaps_rand = getFractionOfGapsAtDistance(master_df_rand, increment=fraction_increment)
        df_gaps_rand["trial"] = i
        df_gaps_rand_master = pd.concat([df_gaps_rand_master, df_gaps_rand], ignore_index=True)
    
    print(df_gaps_rand_master)
    #get the average of the Fraction column at each distance
    df_gaps_rand = df_gaps_rand_master.groupby(["Distance"])["Fraction"].mean().reset_index()
    df_gaps_rand["Fraction_random"] = df_gaps_rand["Fraction"]

    #combine the two dataframes
    df_gaps = pd.concat([df_gaps, df_gaps_rand], ignore_index=True)

    #divide Distance by 1000
    df_gaps["Distance_kb"] = df_gaps["Distance"]/1000

    #group by distance and create two columns: Fraction_real and Fraction_random with 
    df_gaps = df_gaps.groupby(["Distance_kb"])["Fraction"].apply(list).reset_index()
    #split the list into two columns
    df_gaps[["Fraction_real", "Fraction_random"]] = pd.DataFrame(df_gaps["Fraction"].tolist(), index=df_gaps.index)
    df_gaps.drop(columns=["Fraction"], inplace=True)
    df_gaps["Fraction_real_adjusted"] = df_gaps["Fraction_real"] - df_gaps["Fraction_random"]

    return df_gaps

def plotHaplotypeLengthDistributionVsExpected(
        df_gaps: pd.DataFrame, 
        figsize: tuple[int, int] = (12, 6),
        xticks: list[int] = np.arange(0, 51, 5),
        yticks: list[int] = np.arange(0, 0.51, 0.1),
        show_adjusted:bool = False,
        show:bool = False, 
        save_path: str|None = None,
        
        ) -> None:
    '''Plots the fraction of gaps at each distance between haplotypes, 
    and the expected fraction of gaps at each distance between haplotypes.
    The expected fraction is calculated by randomly selecting the same number of 
    haplotypes as in the real dataset, and calculating the fraction of gaps between them.

    Args:
        df_gaps (DataFrame): A DataFrame with columns: Distance_kb, Fraction_real, Fraction_random, Fraction_real_adjusted

    Returns:
        None
    '''
    sns.set_style("whitegrid")
    plt.figure(figsize=figsize)
    sns.set_context(context="talk", font_scale = 1, rc={"lines.linewidth": 3})
    sns.lineplot(data=df_gaps, x="Distance_kb", y="Fraction_real", color="tab:blue", label="Actual")
    sns.lineplot(data=df_gaps, x="Distance_kb", y="Fraction_random", color="tab:orange", label="Random")
    if show_adjusted:
        sns.lineplot(data=df_gaps, x="Distance_kb", y="Fraction_real_adjusted", color="tab:red", label="Real - Random")
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.xlim(0, 50)
    plt.xlabel("Haplotype block length (kb)", fontweight='bold')
    plt.ylabel("Fraction of all haplotype blocks", fontweight='bold')
    plt.title("Cumulative fraction of haplotype block lengths", fontweight='bold')
    plt.legend()
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
        
    plt.close()

def plotHaplotypeLengthDistributionHist(
        df_blocks: pd.DataFrame,
        figsize: tuple[int, int] = (16, 4),
        xticks: list[int] = np.arange(0, 51, 5),
        save_path: str|None = None,
        show: bool = False
        ) -> None:

    sns.set_style("whitegrid")
    sns.set_context(context="talk", font_scale = 1, rc={"lines.linewidth": 2})
    plt.figure(figsize=figsize)
    sns.histplot(data=df_blocks, x="Distance_kb", bins=200, alpha=0.8, color="tab:blue")
    plt.xticks(xticks)
    plt.title("Haplotype block length distribution", fontweight='bold')
    plt.xlabel("Haplotype block length (kb)", fontweight='bold')
    plt.ylabel("Number of haplotype blocks", fontweight='bold')
    plt.xlim(0, None)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()
    
    plt.close()

def save_payload_to_pickle(payload: list, save_path: str = "payload.pkl") -> None:
    '''Saves a payload to a pickle file.

    Args:
        payload (list): A list of objects to save.
        save_path (str): Path to save the payload to.

    Returns:
        None
    '''
    import pickle
    with open(save_path, "wb") as f:
        pickle.dump(payload, f)

    logging.info(f"Saved payload to {save_path}")

def load_payload_from_pickle(load_path: str = "payload.pkl") -> list:
    '''Loads a payload from a pickle file.

    Args:
        load_path (str): Path to load the payload from.

    Returns:
        list: A list of objects.
    '''
    import pickle
    with open(load_path, "rb") as f:
        payload = pickle.load(f)
    
    logging.info(f"Loaded payload from {load_path}")
    return payload

def smoothen_derive_col_data(
        df: pd.DataFrame,
        col_name: str, 
        window_length: int = 10,
        derivative = 0) -> pd.DataFrame:
    """
    Applies a Savitzky-Golay filter to the specified column of the input DataFrame, and optionally computes its derivative.

    Args:
        df (pd.DataFrame): The input DataFrame.
        col_name (str): The name of the column to process.
        window_length (int, optional): The length of the filter window. Defaults to 10.
        derivative (int, optional): The order of the derivative to compute. Defaults to 0 (no derivative).

    Returns:
        pd.DataFrame: The input DataFrame with a new column containing the smoothed and/or derived data.
    """

    from scipy.signal import savgol_filter
    df[f"{col_name}_smoothed"] = savgol_filter(df[col_name], window_length=window_length, polyorder=3, deriv=derivative)
    return df

def plot_col_series(
        df: pd.DataFrame, 
        col_name: str, 
        x_col: str = "Window_Size_kb",
        show: bool = True, 
        show_smoothed_derived: bool = True, 
        show_original: bool = True,
        save_path: str = None,
        title: str | None = None,
        label: str = "Association Rate",
        x_label: str = "Window size (kb)",
        y_label: str = "Rate"
        ) -> None:
    """
    Plots a line chart of the specified column of the input DataFrame, along with its smoothed and/or derived data.

    Args:
        df (pd.DataFrame): The input DataFrame.
        col_name (str): The name of the column to plot.
        show (bool, optional): Whether to show the plot. Defaults to True.
        show_smoothed_derived (bool, optional): Whether to show the smoothed and/or derived data. Defaults to True.
        show_original (bool, optional): Whether to show the original data. Defaults to True.
        save_path (str, optional): The path to save the plot image. Defaults to None.

    Returns:
        None
    """
    #set size
    plt.figure(figsize=(12, 6))
    sns.set_style("whitegrid")
    if show_original:
        sns.lineplot(data=df, x=x_col, y=col_name, color="black", label=col_name)
    if show_smoothed_derived:
        sns.lineplot(data=df, x=x_col, y=f"{col_name}_smoothed", color="purple", linewidth=4, label=label)
    
    # x axis ticks every 5 kb
    plt.xticks(np.arange(0, 51, 5))
    plt.xlabel(x_label, fontweight='bold', fontsize=20)
    plt.ylabel(y_label, fontweight='bold', fontsize=20)
    
    if title is not None:
        plt.title(title, fontweight='bold', fontsize=24)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()

def findNonRepresentativeSNPs(
    samples_parental_SNPs,
    fig_name_freq_all = "outputs/all_parental_snps_combined_frequency_count_distribution_all60_haploung1d.png",
    fig_name_freq_non_repr = "outputs/all_parental_snps_combined_frequency_count_distribution_first57_haploung1d.png",
    non_representative_SNPs_out_path = "outputs/unreliable_SNPs_freq_count.txt") -> None:
    """
    This function takes a list of samples containing SNPs and returns a list of SNPs that are not representative.
    A SNP is considered representative if it is present in 50 or more (out of 60) samples.
    
    Args:

        samples_parental_SNPs (list): a list of samples containing SNPs
        fig_name_freq_all (str): a path to save the figure with the distribution of the Frequency_count for all SNPs
        fig_name_freq_non_repr (str): a path to save the figure with the distribution of the Frequency_count for non-representative SNPs
        non_representative_SNPs_out_path (str): a path to save the list of non-representative SNPs

    Returns:

        None
    """
    from mayo.settings import haploid_samples_ung1d

    df_parental_SNPs_all = pd.DataFrame()

    for sample in samples_parental_SNPs:
        sample_name = sample.name

        if sample.heterozygous:
            continue

        df_parental_SNPs = sample.df.copy()
        df_parental_SNPs = df_parental_SNPs[["Chromosome", "Region", "Allele", "Frequency"]]
        df_parental_SNPs["sample_name"] = sample_name
        #concat to the df_parental_SNPs_all
        df_parental_SNPs_all = pd.concat([df_parental_SNPs_all, df_parental_SNPs])

        #save to file
    df_parental_SNPs_all.to_csv("outputs/all_parental_snps_combined_all_haploid.csv", index=False)

    logging.info(f"The number of haploid samples is: {len(haploid_samples_ung1d)}")

    df_parental_SNPs_all = pd.DataFrame()

    for sample in samples_parental_SNPs:
        sample_name = sample.name

        if sample_name in haploid_samples_ung1d:

            df_parental_SNPs = sample.df.copy()
            df_parental_SNPs = df_parental_SNPs[["Chromosome", "Region", "Allele", "Frequency"]]
            df_parental_SNPs["sample_name"] = sample_name
            #concat to the df_parental_SNPs_all
            df_parental_SNPs_all = pd.concat([df_parental_SNPs_all, df_parental_SNPs])

    #save to file
    df_parental_SNPs_all.to_csv("outputs/all_parental_snps_combined.csv", index=False)

    #group by chromosome and region and Allele, average the frequency, give count of samples
    df_parental_SNPs_all_grouped = df_parental_SNPs_all.groupby(["Chromosome", "Region"]).agg({"Frequency": ["mean", "count"]})

    #collapse the multiindex
    df_parental_SNPs_all_grouped.columns = ["_".join(x) for x in df_parental_SNPs_all_grouped.columns.ravel()]

    #reset the index
    df_parental_SNPs_all_grouped = df_parental_SNPs_all_grouped.reset_index()

    logging.info(df_parental_SNPs_all_grouped.describe())

    #plot the distribution of the Frequency_mean
    plt.figure(figsize=(12, 5))
    sns.distplot(df_parental_SNPs_all_grouped["Frequency_count"], kde=False, bins=60)
    #save to file
    plt.savefig(fig_name_freq_all, dpi=300)

    plt.figure(figsize=(12, 5))
    sns.distplot(df_parental_SNPs_all_grouped[df_parental_SNPs_all_grouped["Frequency_count"] < 58]["Frequency_count"], kde=False, bins=57)
    plt.savefig(fig_name_freq_non_repr, dpi=300)

    df_unreliable = (df_parental_SNPs_all_grouped[df_parental_SNPs_all_grouped["Frequency_count"] < 50])

    #save df_unreliable to a txt file for further analysis
    #iterate over the df_unreliable and save the SNPs to a file
    out_str = ''
    for index, row in df_unreliable.iterrows():
        chromosome = row["Chromosome"]
        region = row["Region"]
        out_str += f"['{chromosome}', {region}],\n"

    with open(non_representative_SNPs_out_path, "w") as f:
        f.write(out_str)

def qqplot(
    x,
    y, 
    quantiles=None, 
    interpolation='nearest', 
    ax=None, rug=False,
    rug_length=0.05, 
    rug_kwargs=None, 
    **kwargs):
    """Draw a quantile-quantile plot for `x` versus `y`.

    Parameters
    ----------
    x, y : array-like
        One-dimensional numeric arrays.

    ax : matplotlib.axes.Axes, optional
        Axes on which to plot. If not provided, the current axes will be used.

    quantiles : int or array-like, optional
        Quantiles to include in the plot. This can be an array of quantiles, in
        which case only the specified quantiles of `x` and `y` will be plotted.
        If this is an int `n`, then the quantiles will be `n` evenly spaced
        points between 0 and 1. If this is None, then `min(len(x), len(y))`
        evenly spaced quantiles between 0 and 1 will be computed.

    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        Specify the interpolation method used to find quantiles when `quantiles`
        is an int or None. See the documentation for numpy.quantile().

    rug : bool, optional
        If True, draw a rug plot representing both samples on the horizontal and
        vertical axes. If False, no rug plot is drawn.

    rug_length : float in [0, 1], optional
        Specifies the length of the rug plot lines as a fraction of the total
        vertical or horizontal length.

    rug_kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.axvline() and
        matplotlib.axes.Axes.axhline() when drawing rug plots.

    kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.scatter() when drawing
        the q-q plot.
    """
    import numbers
    import numpy as np
    # Get current axes if none are provided
    if ax is None:
        ax = plt.gca()

    if quantiles is None:
        quantiles = min(len(x), len(y))

    # Compute quantiles of the two samples
    if isinstance(quantiles, numbers.Integral):
        quantiles = np.linspace(start=0, stop=1, num=int(quantiles))
    else:
        quantiles = np.atleast_1d(np.sort(quantiles))
    x_quantiles = np.quantile(x, quantiles, interpolation=interpolation)
    y_quantiles = np.quantile(y, quantiles, interpolation=interpolation)

    # Draw the rug plots if requested
    if rug:
        # Default rug plot settings
        rug_x_params = dict(ymin=0, ymax=rug_length, c='gray', alpha=0.5)
        rug_y_params = dict(xmin=0, xmax=rug_length, c='gray', alpha=0.5)

        # Override default setting by any user-specified settings
        if rug_kwargs is not None:
            rug_x_params.update(rug_kwargs)
            rug_y_params.update(rug_kwargs)

        # Draw the rug plots
        for point in x:
            ax.axvline(point, **rug_x_params)
        for point in y:
            ax.axhline(point, **rug_y_params)

    # Draw the q-q plot
    ax.scatter(x_quantiles, y_quantiles, **kwargs)

    #add diagonal line
    ax.plot([0,1],[0,1], transform=ax.transAxes, color='red', linestyle='--')

def get_all_SNVs_df(
        frames_list: list,
        filter_sample_names: list = []
        ) -> pd.DataFrame:
    """
    Returns a pandas DataFrame containing all SNVs from a list of frames.

    Args:
        frames_list (list): A list of frames to extract SNVs from.
        filter_sample_names (list, optional): A list of sample names to filter by. Defaults to [].

    Returns:
        pd.DataFrame: A DataFrame containing all SNVs from the input frames.
    """
    #create a dataframe with all SNVs
    df_SNPs_all = pd.DataFrame()
    for i in frames_list:
        
        #only take samples from the list
        if i.name not in filter_sample_names:
            logging.info(f"{i.name} not in sample_names list. skipping...")
            continue
        
        #df_temp = i.df_SNPs.copy()
        df_temp = i.df.copy()
        #remove NaNs and print out the number of NaNs and cleaned dataframe
        len_before = len(df_temp)
        df_temp = df_temp[df_temp["Region"].notna()]
        len_after = len(df_temp)
        logging.info(f"Number of NaNs removed: {len_before - len_after}. Total number of SNVs: {len(df_temp)}")
        #print the number of Is_Clustered==True SNPs
        #print(f"Number of clustered SNPs: {len(df_temp[df_temp['Is_Clustered']==True])}")
        #remove NaNs
        df_temp = df_temp[df_temp["Region"].notna()]
        df_temp["sample_name"] = i.name
        df_temp["Event"] = "SNV"
        df_temp.rename(columns={"Region": "Event_Position"}, inplace=True)
        df_SNPs_all = pd.concat([df_SNPs_all, df_temp], ignore_index=True)
    #print out total number of SNVs
    logging.info(f"Total number of SNVs: {len(df_SNPs_all)}")

    return df_SNPs_all

def get_all_switches_df(
        frames_list: list, 
        filter_sample_names: list = []
        ) -> pd.DataFrame:
    """
    Returns a pandas DataFrame containing all switches from a list of frames.

    Args:
        frames_list (list): A list of frames to extract switches from.
        filter_sample_names (list, optional): A list of sample names to filter by. Defaults to [].

    Returns:
        pd.DataFrame: A DataFrame containing all switches from the input frames.
    """
    #create a dataframe with all switches
    df_switches_all = pd.DataFrame()
    for i in frames_list:

        #only take samples from the list
        if i.name not in filter_sample_names:
            continue

        df_temp = i.switchesTable.copy()
        df_temp["sample_name"] = i.name
        df_temp["Event"] = "Switch"
        df_temp.rename(columns={"Switch_Center": "Event_Position"}, inplace=True)
        df_switches_all = pd.concat([df_switches_all, df_temp], ignore_index=True)

    #print out total number of switches
    logging.info(f"Total number of switches: {len(df_switches_all)}")

    return df_switches_all

def createCombinedSNVsSwitchesDF(
        frames_list: list, 
        filter_sample_names: list = []
        ) -> pd.DataFrame:
    """
    Returns a pandas DataFrame containing all SNVs and switches from a list of frames, with parental SNPs added.

    Args:
        frames_list (list): A list of frames to extract SNVs and switches from.
        filter_sample_names (list, optional): A list of sample names to filter by. Defaults to [].

    Returns:
        pd.DataFrame: A DataFrame containing all SNVs and switches from the input frames, with parental SNPs added.
    """
    df_SNPs_all = get_all_SNVs_df(frames_list, filter_sample_names)
    df_switches_all = get_all_switches_df(frames_list, filter_sample_names)

    #merge the two dataframes
    df_combined = pd.concat([df_SNPs_all, df_switches_all], ignore_index=True)
    #add parental SNP column by merging with the parental SNP table
    parental_origin_db_df = pd.read_csv("outputs/all_parental_SNPs.tsv", sep="\t")
    parental_SNPs_df = parental_origin_db_df.copy()
    parental_SNPs_df = parental_SNPs_df[parental_SNPs_df["Parent"] == "parent_1"]
    parental_SNPs_df.rename(columns={"Region": "Event_Position"}, inplace=True)
    parental_SNPs_df = parental_SNPs_df[["Chromosome", "Event_Position"]]
    parental_SNPs_df["Event"] = "SNP"

    df_combined = pd.concat([df_combined, parental_SNPs_df], ignore_index=True)

    return df_combined

def correlateSNVsWithSwitches(
        frames_list: list, 
        snp_num: int,
        show: bool = True,
        ) -> None:
    '''Correlates SNVs and switches from a list of SampleTable objects representing 
    normal samples, and plots their distribution on each chromosome as a histogram.

    Args:
        frames_list_samples_normal_ref (list): A list of SampleTable objects 
            representing the normal samples to be analyzed.
        snp_num (int): The number of SNPs to used for the analysis.
        show (bool, optional): Whether to show the plots. Defaults to True.

    Returns:
        None
    '''
    import scipy.stats as stats
    from mayo.settings import haploid_samples_ung1d as haploid_samples_ung1d
    df_combined = createCombinedSNVsSwitchesDF(frames_list, haploid_samples_ung1d)
    
    #rename the Chromosomes to chr1, chr2, etc.
    df_combined = renameChrToRoman(df_combined, "Chromosome")
    x_var="Switch"
    y_var="SNV"
    #plot the distribution of SNVs and switches on each chromosome as a histogram (x-axis: chromosome position, y-axis: number of events)
    for i in df_combined["Chromosome"].unique():

        # if chromosome is mitochondrial, skip
        if i == "ref|NC_001224|":
            continue

        #perform a Kolmogorov-Smirnov test on the SNVs and switches
        print("Kolmogorov-Smirnov test for ", i)
        test = stats.ks_2samp(df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]==x_var)]["Event_Position"], df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]==y_var)]["Event_Position"])
        p_value = test[1]
        stat = test[0]
        print(f"p-value: {p_value}, stat: {stat}\n")
        
        #make the histogram of SNVs and switches, on the same plot, with different colors and separate y-axes
        sns.set_style("ticks")
        plt.rcParams['figure.figsize'] = [30, 5]
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        #plot histogram
        df_plot_x = df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]==x_var)]
        df_plot_y = df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]==y_var)]

        ax1.hist(df_plot_x["Event_Position"], bins=200, color="blue", alpha=0.5)
        ax2.hist(df_plot_y["Event_Position"], bins=200, color="red", alpha=0.5)

        #set the y-axis limits
        ax1.set_ylim(0, 50)
        ax2.set_ylim(0, 50)
        ax1.set_ylabel(x_var, color="blue")
        ax2.set_ylabel(y_var, color="red")
        ax1.set_xlabel("Chromosome position")
        plt.title(f'Chr {i}, SNVs and switches ({snp_num} SNPs) histograms; KS test p-value: {p_value:.3e}, stat: {stat:.3e}')
        plt.savefig(f"figures/SNVs_switches_histograms_{snp_num}snp/{i}.png", dpi=300)
        if show:
            plt.show()
        plt.close()

        sns.set_style("ticks")
        #also, draw a 2-sample Q-Q plot between the SNVs and switches (sm.qqplot_2samples)
        plt.figure(figsize=(12, 12))
        u = df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]==x_var)]["Event_Position"]
        v = df_combined[(df_combined["Chromosome"]==i) & (df_combined["Event"]==y_var)]["Event_Position"]

        #to find number of quantiles, take the length of the chrmosome and divide by 10,000
        qqplot(u, v, c='r', alpha=0.5, edgecolor='k', rug=True, rug_length=0.05, s=100)
        plt.xlabel(x_var)
        plt.ylabel(y_var)
        plt.title(f'Chr {i}, Q-Q plot for SNVs and switches ({snp_num} SNPs); KS test p-value: {p_value:.3e}, stat: {stat:.3e}')
        plt.savefig(f"figures/SNVs_switches_{snp_num}snp_QQplots/{i}.png", dpi=300)
        if show:
            plt.show()
        plt.close()

def countSNVsSwitchesInSlidingWindow(
        df_combined: pd.DataFrame,
        window_size: int = 5000, 
        step_size: int = 5000
        ) -> pd.DataFrame:
    """
    Counts the number of SNVs and switches in a sliding window of a specified size and step size.

    Args:
        df_combined (pd.DataFrame): A DataFrame containing SNVs and switches.
        window_size (int, optional): The size of the sliding window. Defaults to 5000.
        step_size (int, optional): The step size of the sliding window. Defaults to 5000.

    Returns:
        pd.DataFrame: A DataFrame containing the counts of SNVs and switches in each sliding window.
    """
    from mayo.settings import chromosomes as chromosomes
    df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])

    #count the number of SNVs and switches total
    num_SNVs = len(df_combined[df_combined["Event"] == "SNV"])
    num_switches = len(df_combined[df_combined["Event"] == "Switch"])
    logging.info(f"Total number of SNVs: {num_SNVs}, total number of switches: {num_switches}")

    #iterate through all chromosomes
    out_df = pd.DataFrame()
    for chromosome in df_chr["chromosome"].unique():

        #create a sliding window of W=5000 slide it along the chromosome and count the number of SNVs and switches in each window
        df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])
        df_chr = df_chr[df_chr["chromosome"] == chromosome]
        end_position = df_chr["end_position"].values[0]
        df_combined_chr = df_combined[df_combined["Chromosome"] == chromosome]
        df_combined_chr = df_combined_chr.sort_values(by=["Event_Position"])
        df_combined_chr.reset_index(inplace=True, drop=True)
        #create a window of size W=5000
        window_start = 0
        window_end = window_start + window_size
        df_window = df_combined_chr[(df_combined_chr["Event_Position"] >= window_start) & (df_combined_chr["Event_Position"] < window_end)]
        #create a dataframe to store the results
        df_window_counts = pd.DataFrame(columns=["Chromosome", "Window_Start", "Window_End", "SNVs", "Switches"])
        #slide the window along the chromosome
        while window_end <= end_position:
            #count the number of SNVs and switches in the window
            num_SNVs = len(df_window[df_window["Event"] == "SNV"])
            num_switches = len(df_window[df_window["Event"] == "Switch"])
            #add the counts to the dataframe thorugh concat
            df_window_counts = pd.concat([df_window_counts, pd.DataFrame({"Chromosome": chromosome, "Window_Start": window_start, "Window_End": window_end, "SNVs": num_SNVs, "Switches": num_switches}, index=[0])], ignore_index=True)
            #slide the window
            window_start += step_size
            window_end += step_size
            df_window = df_combined_chr[(df_combined_chr["Event_Position"] >= window_start) & (df_combined_chr["Event_Position"] < window_end)]
        #add the dataframe to the output dataframe
        out_df = pd.concat([out_df, df_window_counts], ignore_index=True)

    return out_df

def plot_smoothed_SNV_Switches_lineplot(
        out_df: pd.DataFrame,
        snp_num: int,
        show: bool = True,
        ) -> None:
    """
    Plots a smoothed line plot of the number of SNVs and switches in each window for each chromosome.

    Args:
        df_combined (pd.DataFrame): A DataFrame containing SNVs and switches.
        snp_num (int): The number of SNPs.
        save_path (str | None, optional): The path to save the plot to. Defaults to None.
        show (bool, optional): Whether to show the plot. Defaults to True.

    Returns:
        None
    """
    #smooth the data
    out_df = smoothen_derive_col_data(out_df, "SNVs", window_length=100)
    out_df = smoothen_derive_col_data(out_df, "Switches", window_length=100)

    SNVs_col_name = "SNVs_smoothed"
    switches_col_name = "Switches_smoothed"

    #rename the Chromosomes to chr1, chr2, etc.
    out_df = renameChrToRoman(out_df, "Chromosome")

    for i in out_df["Chromosome"].unique():
        #i = "ref|NC_001135|"
        #plot a line plot of the number of SNVs and switches in each window
        sns.set_style("ticks")
        #get max of the x-axis
        x_max = out_df[out_df["Chromosome"] == i]["Window_End"].max()
        x_max_unit = x_max / 10000

        plt.rcParams['figure.figsize'] = [1*x_max_unit, 5]
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        #plot a line plot where window_end is the x-axis and SNVs and switches are the y-axis
        df_plot = out_df[out_df["Chromosome"] == i]
        df_plot = df_plot.sort_values(by=["Window_End"])
        df_plot.reset_index(inplace=True, drop=True)
        ax1.plot(df_plot["Window_End"], df_plot[f"{SNVs_col_name}"], color="red")
        ax1.set_ylabel("SNVs", color="red")
        ax1.set_xlabel("Window End")
        ax2.plot(df_plot["Window_End"], df_plot[f"{switches_col_name}"], color="blue")
        ax2.set_ylabel("Switches", color="blue")
        #set the y-axis limits
        #ax1.set_ylim(0, 50)
        ax2.set_ylim(0, 30)
        ax1.set_ylim(0, 130)
        plt.title(f'Chr {i}, SNVs and switches ({snp_num} SNPs) line plot (smoothed)')
        plt.savefig(f"figures/SNVs_switches_lineplots_{snp_num}snp/{i}.png", dpi=300)

        if show:
            plt.show()
        
        plt.close()

def compute_and_add_JT_clusters(
        samples_SNVs, 
        pvals_to_test, 
        max_distance_between_mutations = 10000, 
        min_cluster_mutations = 3):
    
    from .nbinom import analyzeSampleForMutationClusters

    for pval in pvals_to_test:
        logging.info(f"Testing pval = {pval}")
        pval_master_df = pd.DataFrame()
        for sample in samples_SNVs:
            df_mutations = sample.df.copy()[["Chromosome", "Region"]]
            new_clusters = analyzeSampleForMutationClusters(
                df_mutations,
                min_p_value = pval,
                max_distance_between_mutations = max_distance_between_mutations,
                min_cluster_mutations = min_cluster_mutations)
            
            sample.cluster_dict_JT.update({pval: new_clusters})
            new_clusters["sample"] = sample.name
            new_clusters["all_mutations"] = len(sample.df)

            pval_master_df = pd.concat([pval_master_df, new_clusters])
            #reset the index so that it is unique
            pval_master_df.reset_index(drop=True, inplace=True)
            logging.info(f"{sample.name} has {len(sample.df):3} mutations and {len(new_clusters):3} JT clusters")

        #ID each row (cluster) by index into "Dataset_Cluster_ID" column
        pval_master_df["Dataset_Cluster_ID"] = pval_master_df.index
        #since the first index is 0, add 1 to each index to avoid confusion with 0
        pval_master_df["Dataset_Cluster_ID"] = pval_master_df["Dataset_Cluster_ID"] + 1

        pval_master_df.to_csv(f"data/new_clusters/mutation_clusters/mutation_clusters_pval_{str(pval)}.tsv", sep="\t", index=False)

def add_PMACD_clusters(
        samples_SNVs, 
        PMACD_pvals_to_add):
    
    from mayo.clusters import aggregate_anz5
    for pval in PMACD_pvals_to_add:
        #open the file
        PMACD_file_path = f"data/cluster_files/clusters_diff_thresholds/master_{pval}/1_master_para_bed_sorted_anz5.txt"
        PMACD_df = pd.read_csv(PMACD_file_path, sep="\t")
        for sample in samples_SNVs:

            #filter the dataframe to the sample
            PMACD_df_sample = PMACD_df[PMACD_df["Tumor_Sample_Barcode"] == sample.name].copy()

            # parse the dataframe (currently marked SNPs)
            PMACD_df_sample_clusters = aggregate_anz5(PMACD_df_sample)

            logging.info(f"There are {len(PMACD_df_sample_clusters)} PMACD clusters in {sample.name} for {pval} pval")

            sample.cluster_dict_PMACD.update({pval: PMACD_df_sample_clusters})

def find_commonly_heterozygous_chromosomes(samples_parental_SNPs: list[SampleTable]) -> tuple:
    """
    This function finds the most commonly heterozygous and aneuploid chromosomes from a list of samples.

    Args:
        samples_parental_SNPs (list[SampleTable]): A list of SampleTable objects. Each SampleTable object should have 
        'heterozygous_chromosomes' and 'aneuploid_chromosomes' attributes which are lists of chromosomes.

    Returns:
        tuple: A tuple containing two Counter objects. The first Counter object contains the count of each 
        heterozygous chromosome across all samples. The second Counter object contains the count of each aneuploid 
        chromosome across all samples.

    Raises:
        Exception: If a SampleTable object does not have an 'aneuploid_chromosomes' attribute, an exception is raised 
        and a message is printed.
    """
    from collections import Counter

    all_chr_list_heterozygous = []
    all_chr_list_aneouploid = []
    for sample in samples_parental_SNPs:
        all_chr_list_heterozygous.extend(sample.heterozygous_chromosomes)
        try:
            all_chr_list_aneouploid.extend(sample.aneuploid_chromosomes)
        except Exception as e:
            print(f"Sample {sample.name} has no aneuploid_chromosomes")

    #find how many times each chromosome is heterozygous
    heterozygous_chr_counter = Counter(all_chr_list_heterozygous)

    #find how many times each chromosome is aneuploid
    aneuploid_chr_counter = Counter(all_chr_list_aneouploid)

    #find most commonly heterozygous chromosomes
    most_common_heterozygous_chr = heterozygous_chr_counter.most_common(10)

    #find most commonly aneuploid chromosomes
    most_common_aneuploid_chr = aneuploid_chr_counter.most_common(10)

    print(f"out of {len(samples_parental_SNPs)} samples: Most commonly heterozygous chromosomes:")
    print(most_common_heterozygous_chr)
    print(f"out of {len(samples_parental_SNPs)} samples: Most commonly aneuploid chromosomes:")
    print(most_common_aneuploid_chr)

    return heterozygous_chr_counter, aneuploid_chr_counter

def getTotalssDNAinSamples_SNVs(
        samples_SNVs: List[OverlapSwitchesTable], 
        cluster_kind: Literal["JT", "PMACD"], 
        focus_set: List[str] | None = None):
    """
    Calculate the total ssDNA for a list of samples.

    Parameters:
    samples_SNVs (List[OverlapSwitchesTable]): A list of samples, each represented as an OverlapSwitchesTable.
    cluster_kind (Literal["JT", "PMACD"]): The kind of cluster to be used in the calculation. Can be either "JT" or "PMACD".
    focus_set (List[str] | None): A list of sample names to be included in the calculation. If None, all samples are included.

    Returns:
    df_total_ssDNA (pd.DataFrame): A DataFrame containing the total ssDNA for each sample. The DataFrame has columns "sample", "p_val", and "total_ssDNA".
    """    
    output_list_all = []
    for sample in samples_SNVs:
        if focus_set is not None and sample.name not in focus_set:
            continue
        else:
            output_list_sample = sample.calculateTotalssDNA(cluster_kind=cluster_kind)
            output_list_all.extend(output_list_sample)

    df_total_ssDNA = pd.DataFrame(output_list_all, columns=["sample", "p_val", "total_ssDNA"])

    return df_total_ssDNA

def showTotalssDNAinSamples_SNVs(df_total_ssDNA, p_val) -> None:
    """
    Display the average total amount of ssDNA in clusters with a specific p-value.

    Parameters:
    df_total_ssDNA (pd.DataFrame): A DataFrame containing the total ssDNA for each sample. 
                                   The DataFrame should have columns "sample", "p_val", and "total_ssDNA".
    p_val (float): The p-value to filter the DataFrame on.

    Returns:
    None
    """  
    df_total_ssDNA = df_total_ssDNA[df_total_ssDNA["p_val"] == p_val]

    #take the average of the total_ssDNA column
    print(f"Stats for p_val = {p_val}:")
    print(f"  Average total ssDNA: {df_total_ssDNA['total_ssDNA'].mean()}, median: {df_total_ssDNA['total_ssDNA'].median()}, min: {df_total_ssDNA['total_ssDNA'].min()}, max: {df_total_ssDNA['total_ssDNA'].max()}")

def find_GC_CO_from_switches(samples_SNVs: List[OverlapSwitchesTable], threshold_GC_max = 5000):
    for sample in samples_SNVs:
        sample.classify_GC_CO(threshold_GC_max = threshold_GC_max)

def plot_common_mutations_histogram(
        df_all_mutations: pd.DataFrame,
        max_duplicates: int,
        bins: int = 100,
        savename: str = "figures/common_mutations_histogram.png",
        show: bool = False
        ) -> None:
    """
    Plots a histogram of the number of samples each mutation appears in.

    Parameters:
    df_all_mutations (pd.DataFrame): A DataFrame containing all mutations.
    max_duplicates (int): The maximum number of samples a mutation can appear in to be kept.

    Returns:
    None
    """  
    #plot the distribution of mutations as a histogram
    plt.figure(figsize=(10,5))
    plt.title("Distribution of mutations")
    plt.xlabel("Number of samples (duplicates)")
    plt.ylabel("Count")
    plt.hist(df_all_mutations["mutation"].value_counts(), bins=bins, log=True)
    #draw a vertical line at the number of samples that we want to filter out, shade the area to the right of the line
    plt.axvline(x=max_duplicates, color='r', linestyle='--')
    max_x = plt.gca().get_xlim()[1]
    plt.axvspan(max_duplicates, max_x, alpha=0.3, color='red')
    #save to file
    plt.savefig(savename, dpi=300)
    if show:
        plt.show()
    else:
        plt.close()

def createCombinedSNVsTable(samples_SNVs: List[OverlapSwitchesTable]) -> pd.DataFrame:
    df_all_mutations = pd.DataFrame()
    for sample in samples_SNVs:
        df = sample.df.copy()
        df["sample_name"] = sample.name
        df_all_mutations = pd.concat([df_all_mutations, df], ignore_index=True)
    return df_all_mutations

def filterCommonMutations(
        samples_SNVs: List[OverlapSwitchesTable], 
        max_duplicates: int, 
        genotypes_dict: dict = None,
        specific_genotypes: list = None,
        max_genotype_duplicates: int = None,
        ) -> None:
    """
    Filters out common mutations that appear in more than a specified number of samples and optionally, 
    more than a specified number of genotypes. The function also generates a histogram plot showing the 
    distribution of mutations across samples. Additionally, the function creates a new column ("mutation") 
    in each sample describing the mutation.

    Parameters:
    samples_SNVs (List[OverlapSwitchesTable]): A list of samples, each containing a DataFrame of SNVs.
    max_duplicates (int): The maximum number of samples a mutation can appear in to be kept.
    genotypes_dict (dict, optional): A dictionary mapping genotypes to samples. Default is None.
    max_genotype_duplicates (int | None, optional): The maximum number of genotypes a mutation can appear in to be kept. Default is None.

    Returns:
    None. The function modifies the input samples_SNVs in-place.

    This function also generates a histogram plot showing the distribution of mutations across samples.
    Mutations appearing in more than 'max_duplicates' samples are highlighted.

    Note: This function uses logging to output the progress and results of the filtering process.
    """

    for sample in samples_SNVs:
        sample.df["mutation"] = sample.df["Chromosome"].astype(str) + "_" + sample.df["Region"].astype(str) + "_" + sample.df["Reference"] + "_" + sample.df["Allele"]


    df_all_mutations = createCombinedSNVsTable(samples_SNVs)
    plot_common_mutations_histogram(df_all_mutations, max_duplicates, savename="figures/common_mutations_histogram_before_filtering.png")

    all_raw_mutation_count = len(df_all_mutations)
    logging.info(f"There are {all_raw_mutation_count} raw mutations in total in all samples.")
    

    #iterate through each genotype and remove mutations that appear in more than max_duplicates samples
    if genotypes_dict is not None:
        genotypes_dict_by_genotype = invert_dict(genotypes_dict)

        for genotype in genotypes_dict_by_genotype.keys():

            if specific_genotypes is None:
                logging.info(f"Filtering genotype {genotype} for extra filtering.")

            elif genotype not in specific_genotypes:
                logging.info(f"Skipping genotype {genotype} for extra filtering because it is not in the list of specific genotypes.")
                continue

            samples_genotype = genotypes_dict_by_genotype[genotype]
            df_all_mutations_genotype = df_all_mutations[df_all_mutations["sample_name"].isin(samples_genotype)]
            df_dups_genotype = df_all_mutations_genotype["mutation"].value_counts()[df_all_mutations_genotype["mutation"].value_counts() > max_genotype_duplicates]
            
            for sample in samples_SNVs:

                if sample.name in samples_genotype:

                    len_before = len(sample.df)
                    sample.df = sample.df[~sample.df["mutation"].isin(df_dups_genotype.index)]
                    len_after = len(sample.df)
                    removed_sample = len_before - len_after
                    logging.info(f"Genotype {genotype} removed {removed_sample:3} mutations from {sample.name}. {len_after} mutations left.")

    df_all_mutations_filtered_genotype = createCombinedSNVsTable(samples_SNVs)
    plot_common_mutations_histogram(df_all_mutations_filtered_genotype, max_duplicates, savename="figures/common_mutations_histogram_after_genotype_filtering.png")

    genotype_removed_count = all_raw_mutation_count - len(df_all_mutations_filtered_genotype)
    logging.info(f"Total removed by genotype: {genotype_removed_count}. {len(df_all_mutations_filtered_genotype)} mutations ({(len(df_all_mutations_filtered_genotype) / all_raw_mutation_count * 100):.5f}%) left.")

    #remove mutations that appear in more than max_duplicates samples
    logging.info(f"Removing mutations that appear in more than {max_duplicates} samples.")
    total_removed = 0
    df_dups = df_all_mutations["mutation"].value_counts()[df_all_mutations["mutation"].value_counts() > max_duplicates]
    for sample in samples_SNVs:
        len_before = len(sample.df)
        sample.df = sample.df[~sample.df["mutation"].isin(df_dups.index)]
        len_after = len(sample.df)
        removed_sample = len_before - len_after
        total_removed += removed_sample
        try:
            genotype = genotypes_dict[sample.name]
        except KeyError:
            genotype = "unknown"
        logging.info(f"Removed {removed_sample:3} mutations from {sample.name} ({genotype}). {len_after} mutations left.")


    df_all_mutations_filtered = createCombinedSNVsTable(samples_SNVs)
    plot_common_mutations_histogram(df_all_mutations_filtered, max_duplicates, bins=max_duplicates, savename="figures/common_mutations_histogram_after_all_filters.png")

    #print top 30 value counts of df_all_mutations_filtered from most to least common
    print(df_all_mutations_filtered["mutation"].value_counts()[:30])

    grand_total_removed = all_raw_mutation_count - len(df_all_mutations_filtered)
    logging.info(f"Grand total removed: {grand_total_removed}. {len(df_all_mutations_filtered)} mutations ({(len(df_all_mutations_filtered) / all_raw_mutation_count * 100):.5f}%) left.")

def filterCommonMutations_sectors(
        df_all_mutations: pd.DataFrame, 
        specific_genotypes: list[list] = None,
        minimum_count_within_genotype: int | None = 1,
        max_outside_genotype_duplicates: int | None = 0
        ) -> None:
    """
    Filters a DataFrame of mutations based on specific genotypes and counts.

    Parameters:
    df_all_mutations (pd.DataFrame): The DataFrame containing all mutations.
    specific_genotypes (list[list], optional): A list of specific genotypes to filter by. Default is None.
    minimum_count_within_genotype (int, optional): The minimum count of a mutation within a genotype to be included. Default is 1.
    max_outside_genotype_duplicates (int, optional): The maximum count of a mutation outside a genotype to be included. Default is 0.

    Returns:
    pd.DataFrame: The filtered DataFrame of mutations.
    """
    df_all_mutations = df_all_mutations[df_all_mutations["genotype"].isin(specific_genotypes)].copy()

    df_all_mutations['count_within_spore'] = df_all_mutations.groupby(['mutation'])['mutation'].transform('count')

    # Group the DataFrame by genotype and mutation, and calculate the sum of count_within_genotype for each group
    df_all_mutations['count_within_colony'] = df_all_mutations.groupby(['genotype', 'mutation'])['mutation'].transform('count')

    df_all_mutations['common_between_colonies'] = df_all_mutations['count_within_spore'] - df_all_mutations['count_within_colony']


    #print the number of colonies common between colonies and the number of mutations within each colony
    #print(df_all_mutations['common_between_colonies'].value_counts())
    #print(df_all_mutations['count_within_colony'].value_counts())

    #filter the all_mutations_df to only include mutations


    if minimum_count_within_genotype is not None:
        df_all_mutations = df_all_mutations[(df_all_mutations["count_within_colony"] >= minimum_count_within_genotype)]
    if max_outside_genotype_duplicates is not None:
        df_all_mutations = df_all_mutations[(df_all_mutations["common_between_colonies"] <= max_outside_genotype_duplicates)]

    return(df_all_mutations)

def invert_dict(dictionary: dict) -> dict:
    """
    Inverts a dictionary, so that the keys become values and the values become keys.
    It keeps all the duplicated values, and puts them in a list.

    Parameters:
    dictionary (dict): Dictionary to invert.

    Returns:
    dict: Inverted dictionary.
    """
    
    inverted_dict = dict()

    for key, value in dictionary.items():
        inverted_dict.setdefault(value, list()).append(key)

    return inverted_dict

def get_samples_with_genotype(query_genotype: str, genotype_dict: dict) -> list:
    """
    Returns a list of samples with a given genotype.

    Parameters:
    query_genotype (str): Genotype/s to look for. Can be a single genotype or a combination of genotypes separated by a "+" sign.
    genotype_dict (dict): Dictionary containing sample names as keys and genotypes as values.

    Returns:
    list: A list of sample names with the given genotype.
    """
    genotype_dict_inv = invert_dict(genotype_dict)
    samples_with_genotype = sum([genotype_dict_inv[i] for i in query_genotype.split("+")], [])

    return samples_with_genotype

def calculate_GC_CO(samples_SNVs: List[OverlapSwitchesTable]) -> None:
    """
    Calculate the total number of GC and CO events for each sample.

    This function iterates over each sample in the provided list. For each sample, it calculates the total number of GC 
    (Gene Conversion) and CO (Crossover) events by summing up the lengths of the lists associated with each chromosome 
    in the sample's GC_dict and CO_dict respectively. The results are stored in the sample's total_GC and total_CO attributes.

    Parameters:
    samples_SNVs (List[OverlapSwitchesTable]): A list of samples. Each sample is expected to have GC_dict and CO_dict 
    attributes, which are dictionaries where the keys are chromosome names and the values are lists of events.

    Returns:
    None
    """
    for sample in samples_SNVs:
        
        total_GC = 0
        total_CO = 0

        for chr in sample.GC_dict.keys():
            total_GC += len(sample.GC_dict[chr])

        for chr in sample.CO_dict.keys():
            total_CO += len(sample.CO_dict[chr])
        
        sample.total_GC = total_GC
        sample.total_CO = total_CO

def pull_mutations_relative_to_position(df, chromosome, position, window_size=10000) -> pd.DataFrame:
    """
    Pulls mutations from a dataframe that are within a given window size of a given position.

    Parameters:
    df (DataFrame): DataFrame containing mutations.
    position (int): Position to find mutations relative to.
    window_size (int): Window size to look for mutations in. Default is 1000.

    Returns:
    DataFrame: A DataFrame containing mutations within the window size of the given position in a relative coordinate system.
    (-window_size; upstream, +window_size; downstream)
    """
    df_chr = df[df["Chromosome"] == chromosome]
    position_min = position - window_size
    position_max = position + window_size
    df_chr_pos = df_chr[(df_chr["Region"] >= position_min) & (df_chr["Region"] <= position_max)].copy()
    df_chr_pos["Original_Region"] = df_chr_pos["Region"]
    df_chr_pos["Original_Chromosome"] = chromosome
    df_chr_pos["Region"] = df_chr_pos["Region"] - position

    return df_chr_pos

def pull_mutations_from_sample_recev_dict(ev_dict: dict, df, window_size=10000) -> pd.DataFrame:
    """
    Pulls mutations from a sample's recombination event dictionary.

    Parameters:
    ev_dict (dict): Recombination event dictionary.

    Returns:
    DataFrame: A DataFrame containing mutations from a sample's recombination event dictionary.
    """
    df_out = pd.DataFrame()

    for chromosome in ev_dict.keys():
        for position in ev_dict[chromosome]:
            df_position = pull_mutations_relative_to_position(df, chromosome, position, window_size)
            #create event_id column
            event_ID = f"{chromosome}_{position}"
            df_position["RecEv_ID"] = event_ID
            df_out = pd.concat([df_out, df_position], ignore_index=True)

    return df_out

def cumulate_column(
        df: pd.DataFrame,
        column_name: str,
        col_range: tuple | None = None,
        window_size: int = 100, 
        window_slide: int = 100,
        ) -> pd.DataFrame:
    """
    This function applies a sliding window to a specified column of a DataFrame and counts the number of features in each window.

    Parameters:
    df (pd.DataFrame): The DataFrame to process.
    column_name (str): The name of the column to apply the sliding window to.
    col_range (tuple, optional): The range of values to consider in the column. If None, the function will use the min and max values in the column.
    window_size (int, optional): The size of the sliding window. Default is 100.
    window_slide (int, optional): The amount to slide the window for each step. Default is 100.

    Returns:
    pd.DataFrame: A new DataFrame with the original column values and their corresponding counts within the sliding window.
    """

    df_sliding_window = pd.DataFrame(columns=[column_name, "count"])

    if col_range is None:
        col_range = (df[column_name].min(), df[column_name].max())
    
    range_start, range_end = col_range


    for i in range(range_start, range_end, window_slide):

        df_window = df[(df[column_name] >= i) & (df[column_name] <= i+window_size)].copy()
        total_features = len(df_window)
        concat_df = pd.DataFrame([[i, total_features]], columns=[column_name, "count"])
        df_sliding_window = pd.concat([df_sliding_window, concat_df], ignore_index=True)

    #window_size smoothing shifts the data to the left, so we need to shift it back to the right
    df_sliding_window[column_name] = df_sliding_window[column_name] + (window_size / 2)

    return df_sliding_window

def getA3A_motifs(seq_record):
    """
    This function counts the occurrences of specific motifs (TCT, TCA, TCG, AGA, TGA, CGA) in a given sequence record.
    The sequence record is assumed to be a tuple where the second element is the sequence.

    Parameters:
    seq_record (tuple): A tuple where the second element is the sequence to search for motifs.

    Returns:
    tuple: A tuple containing the counts of each motif in the sequence. The order of counts is as follows:
           (TCT count, TCA count, TCG count, AGA count, TGA count, CGA count)
    """
    fasta_seq = seq_record[1]

    #forward strand:
    tct_count = fasta_seq.count("TCT")

    tca_count = fasta_seq.count("TCA")

    tcg_count = fasta_seq.count("TCG")

    #reverse strand:
    aga_count = fasta_seq.count("AGA")

    tga_count = fasta_seq.count("TGA")

    cga_count = fasta_seq.count("CGA")

    return tct_count, tca_count, tcg_count, aga_count, tga_count, cga_count

def map_genotype_to_sample(samples_SNVs, genotype_dict):
    """
    Maps the genotype of each sample in samples_SNVs to the corresponding genotype in genotype_dict.
    """
    for sample in samples_SNVs:
        try:
            sample.genotype = genotype_dict[sample.name]
        except KeyError:
            print(f"Genotype for {sample.name} not found in genotype_dict. Setting to 'unknown'")
            sample.genotype = "unknown"

def constrain_to_region(
        df: pd.DataFrame,
        region: tuple = (-12, 12)
        ) -> pd.DataFrame:
    """
    Constrains the DataFrame to a specified region.

    Parameters:
    df (pd.DataFrame): DataFrame to be constrained.
    region (tuple): Tuple specifying the lower and upper bounds of the region. Default is (-12, 12).

    Returns:
    pd.DataFrame: DataFrame constrained to the specified region.
    """
    return(df[(df["Region"] >= region[0]) & (df["Region"] <= region[1])].copy())

def plot_SNVs_around_rec_ev(
        master_relative_GC_CO_df: pd.DataFrame, 
        title: str = "",
        save_path: str | None = None,
        show: bool = False,
        mutation_type: Literal["clustered", "scattered", "all"] | None = None,
        window_size = 10000
        ) -> None:
    """
    Plots the distribution of Single Nucleotide Variants (SNVs) around recombination events.

    Parameters:
    master_relative_GC_CO_df (pd.DataFrame): DataFrame containing information about the SNVs and their relative positions to recombination events.
    title (str): Title for the plot. Default is an empty string.
    save_path (str | None): Path to save the plot. If None, the plot is not saved. Default is None.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is False.

    Returns:
    None
    """
    import numpy as np
    master_relative_GC_CO_df = master_relative_GC_CO_df.copy()
    #rename the C->T and G->A mutations to C→T and G→A
    master_relative_GC_CO_df["Reference"] = master_relative_GC_CO_df["Reference"].replace({
        "C->T": "C→T", 
        "G->A": "G→A",
        "C->T or G->A": "C→T or G→A"})
    pallette = {
        "C→T": "tab:red", 
        "G→A": "tab:blue",
        "C→T or G→A": "tab:green"}
    pallette = {key: pallette[key] for key in master_relative_GC_CO_df["Reference"].unique()}

    #plot the distribution of mutations relative to GC and CO events
    sns.set_context("poster")
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(12, 8))
    hue_order = ["C→T", "G→A", "C→T or G→A"]
    hue_order = [i for i in hue_order if i in master_relative_GC_CO_df["Reference"].unique()]
    alpha = 0.2
    if "C→T or G→A" in hue_order:
        alpha = 0.3
    sns.histplot(
        ax=ax,
        data=master_relative_GC_CO_df, 
        x="Region",
        hue="Reference", 
        palette=pallette, 
        kde=False, 
        bins=100, 
        kde_kws={'bw_adjust': .3},
        stat='density',
        hue_order=hue_order, 
        legend=False,
        alpha=alpha)
    
    #rename Reference to Mutation Type
    master_relative_GC_CO_df.rename(columns={"Reference": "Mutation Type"}, inplace=True)
    sns.kdeplot(
        ax=ax,
        data=master_relative_GC_CO_df, 
        x="Region", 
        hue="Mutation Type", 
        palette=pallette, 
        bw_adjust=0.2, 
        linewidth=4,   # Set the desired line thickness
        hue_order=hue_order,
        fill=False,
        common_norm=True,
        legend=True,
        #alpha=alpha
        )

    # Convert the y-axis to percentage
    ax = plt.gca()
    y_ticks = ax.get_yticks()
    ax.set_yticklabels([f'{int(y * 100)}' for y in y_ticks])

    plt.axvline(x=0, color="black", linestyle="--", linewidth=2) # draw a vertical line at 0

    if mutation_type is not None:
        plt.title(f"{mutation_type.capitalize()} mutations relative to recombination events", fontweight="bold", fontsize=28)
    else:
        plt.title(f"Mutations relative to recombination events", fontweight="bold", fontsize=28)

    #add a genotype label in the upper left corner
    if title:
        plt.text(0.05, 0.98, f"$\it{{{title}}}$", transform=ax.transAxes, fontsize=24, verticalalignment="top")
    
    plt.xlim(-window_size/1000, window_size/1000)
    plt.xticks(np.arange(-window_size/1000, window_size/1000+1, 2))
    plt.xlabel("Distance (kb) from DSB axis", fontweight="bold", fontsize=24)
    plt.ylabel("Mutation Density", fontweight="bold", fontsize=24)

    #replace legend title to Mutation Type instead of Reference
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles=handles[1:], labels=labels[1:], title="Mutation Type", title_fontsize=20, fontsize=20, loc="upper right")

    #remove y axis ticks
    ax.yaxis.set_major_locator(plt.NullLocator())

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

def plot_SNVs_around_rec_ev_meta(
        meta_master_relative_GC_CO_df, 
        title: str = "",
        save_path: str | None = None,
        show: bool = True
        ):
    import numpy as np
    """
    Plots the distribution of Single Nucleotide Variants (SNVs) around recombination events for a meta-analysis.

    Parameters:
    meta_master_relative_GC_CO_df (pd.DataFrame): DataFrame containing information about the SNVs and their relative positions to recombination events for a meta-analysis.
    title (str): Title for the plot. Default is an empty string.
    save_path (str | None): Path to save the plot. If None, the plot is not saved. Default is None.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is True.

    Returns:
    None
    """
    sns.set_context("talk")
    plt.figure(figsize=(16, 10))
    sns.set_style("whitegrid")
    #sns.histplot(data=meta_master_relative_GC_CO_df, x="Region", hue="genotype", kde=True, bins=100, kde_kws={'bw_adjust': .3})
    sns.kdeplot(data=meta_master_relative_GC_CO_df, x="Region", hue="genotype", fill=False, common_norm=False, bw_adjust=.5)
    plt.axvline(x=0, color="black", linestyle="--", linewidth=1)
    #cut the plot at -12 and 12
    plt.xlim(-12, 12)
    plt.title(f"Mutations relative to recombination events ({title})")
    plt.xticks(np.arange(-12, 13, 1))
    plt.xlabel("Distance from event (kb)")
    plt.ylabel("Frequency")

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

def calculateClusterStatistics(
        samples_SNVs: list[OverlapSwitchesTable],
        cluster_kind: Literal["JT", "PMACD"] = "JT",
        focus_set: list = None,
        pval_list = [0.01, 0.001, 0.0001, 0.00001, 0.000001]
        ) -> tuple[List, List]:
    """
    Calculates some cluster statistics for a given set of samples_SNVs.

    Parameters:
    samples_SNVs (list): List of samples.
    cluster_kind (str): Type of cluster to consider. Default is "JT". Other option is "PMACD".
    focus_set (list): List of sample names to focus on. If None, all samples are considered. Default is None.
    pval_list (list): List of p-values to consider for the analysis. Default is [0.01, 0.001, 0.0001, 0.00001, 0.000001].

    Returns:
    tuple: A tuple containing two lists. The first list contains the ratios of clustered to scattered mutations for each p-value. 
           The second list contains the total number of clustered mutations and total ssDNA for each p-value.
    """
    output_ratios_list_pvals = []
    clust_mut_vs_total_ssdna_list = []

    for pval in pval_list:

        output_ratios_list = []

        for sample in samples_SNVs:

            if focus_set is not None and sample.name not in focus_set:
                continue

            sample_total_mutation_count = len(sample.df)

            if cluster_kind == "JT":
                sample_total_clustered_mutations = sample.cluster_dict_JT[pval]["Mutations"].sum()
                total_ssDNA = sample.totalssDNA_JT_dict[pval]
            
            elif cluster_kind == "PMACD":
                #this is not tested yet
                sample_total_clustered_mutations = sample.cluster_dict_PMACD[pval]["Mutations"].sum()
                total_ssDNA = sample.totalssDNA_PMACD_dict[pval]

            sample_scattered_mutations = sample_total_mutation_count - sample_total_clustered_mutations
            clustered_scattered_ratio = (sample_total_clustered_mutations / (sample_scattered_mutations + sample_total_clustered_mutations)) * 100

            output_ratios_list.append(clustered_scattered_ratio)
            clust_mut_vs_total_ssdna_list.append([pval, sample_total_clustered_mutations, total_ssDNA])

            print(f"pval={pval}; {cluster_kind} clusters stats: {sample.name} mutations; total: {sample_total_mutation_count}, clustered: {sample_total_clustered_mutations}, scattered: {sample_scattered_mutations}, ratio: {clustered_scattered_ratio:.4f}, total ssDNA: {total_ssDNA}")

        output_ratios_list_pvals.append(output_ratios_list)

    return output_ratios_list_pvals, clust_mut_vs_total_ssdna_list

def get_proportion_of_A3A_mutations(samples_SNVs: list[OverlapSwitchesTable]) -> pd.DataFrame:
    """
    Calculates the proportion of A3A mutations (C->T and G->A) vs all mutations for each sample in samples_SNVs.

    Parameters:
    samples_SNVs (list): List of samples.

    Returns:
    DataFrame: A DataFrame containing the proportion of A3A mutations vs all mutations for each sample.
    """
    output_df = pd.DataFrame()
    for sample in samples_SNVs:
        total_mutations = len(sample.df)
        C_to_T = len(sample.df[(sample.df["Reference"] == "C") & (sample.df["Allele"] == "T")])
        G_to_A = len(sample.df[(sample.df["Reference"] == "G") & (sample.df["Allele"] == "A")])
        try:
            total_a3a = C_to_T + G_to_A
            a3a_fraction = total_a3a / total_mutations
        except ZeroDivisionError:
            total_a3a = 0
        #concatenate to output_df
        output_df = pd.concat([
            output_df, 
            pd.DataFrame({
                "Sample": [sample.name], 
                "Total mutations": [total_mutations], 
                "C->T": [C_to_T], 
                "G->A": [G_to_A], 
                "Total A3A": [total_a3a], 
                "A3A fraction": [a3a_fraction],
                "Genotype": [sample.genotype]
                })])
        
    output_df = output_df.reset_index(drop=True)
    return output_df

def plot_proportion_A3A_mutations(
        A3A_mutations_proportion_df: pd.DataFrame, 
        genotypes_include: list,
        save_path: str | None = "figures/A3A_fraction_boxplot_fig2.png",
        show: bool = True,
        fig_size: tuple = (12, 5)
        ) -> None:
    """
    Plots the proportion of A3A mutations as a box plot with a stripplot overlay.

    Parameters:
    A3A_mutations_proportion_df (pd.DataFrame): DataFrame containing the proportion of A3A mutations.
    genotypes_include (list): List of genotypes to include in the plot.
    save_path (str | None): Path to save the plot. If None, the plot is not saved. Default is "figures/A3A_fraction_boxplot_fig2.png".
    show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is True.

    Returns:
    None
    """
    df = A3A_mutations_proportion_df[A3A_mutations_proportion_df["Genotype"].isin(genotypes_include)].copy()
    replace_dict = {
        'ung1∆ non-selected': 'ung1∆\nnon-selected',
        'inc. Tetrad ung1∆ non-selected': 'inc. Tetrad\nung1∆\nnon-selected',
        'ung1∆ premeiotic non-selected': 'ung1∆\npremeiotic\nnon-selected',
        'spo13∆spo11∆ premeiotic': 'ung1∆\nspo13∆\nspo11∆\npremeiotic',
        'exo1-ndpol32∆': 'ung1∆\nexo1-nd\npol32∆',
        'exo1-ndsgs1∆C': 'ung1∆\nexo1-nd\nsgs1∆C',
        'spo13∆spo11∆': 'ung1∆\nspo13∆\nspo11∆',
        'exo1-nd': 'ung1∆\nexo1-nd',
        'sgs1∆C': 'ung1∆\nsgs1∆C',
        'pol32∆': 'ung1∆\npol32∆',
        'spo13∆': 'ung1∆\nspo13∆',
    }

    #rename genotypes for plotting so they fit on the x-axis
    df["Genotype"] = df["Genotype"].replace(replace_dict)

    #replace in genotype_include list as well
    genotypes_include = [replace_dict.get(i, i) for i in genotypes_include]


    #print out the A3A fraction (average) for each genotype
    print(df.groupby("Genotype")["A3A fraction"].mean())

    #plot the A3A fraction as a box plot with a stripplot overlay
    sns.set_context("talk")
    plt.figure(figsize=fig_size)
    sns.set_style("whitegrid")
    #set order according to genotypes_include
    sns.boxplot(data=df, x="Genotype", y="A3A fraction", color="lightgray", showfliers=False, linewidth=2, order=genotypes_include)
    import betterbeeswarm
    sns.swarmplot(data=df, overflow='random', x="Genotype", y="A3A fraction", color="white", size=6, edgecolor="black", order=genotypes_include, alpha=0.6, linewidth=1)
    # plt.setp(plt.gca().collections, edgecolor="black", linewidth=1)
    # for collection in plt.gca().collections:
    #     collection.set_edgecolor((0, 0, 0, 1))  # RGBA tuple where A is the alpha channel
    #     collection.set_facecolor((1, 1, 1, 0)) 
    #sns.stripplot(data=df, x="Genotype", y="A3A fraction", color="white", size=6, edgecolor="black", jitter=0.25, linewidth=0.5, order=genotypes_include, alpha=0.7)
    plt.title(f"A3A fraction", fontsize=20, fontweight="bold")
    plt.xlabel("Sample", fontsize=16, fontweight="bold")
    plt.ylabel("A3A fraction", fontsize=16, fontweight="bold")
    plt.xticks(rotation=0, ha="center", fontsize=14, fontstyle="italic", fontweight="bold")
    plt.yticks(fontsize=14, fontweight="bold")

    #save the figure
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    
    if show:
        plt.show()
    else:
        plt.close()

def bin_rec_relative_mutations(
        df: pd.DataFrame,
        region: tuple = (-12, 12),
        bins: int = 24,
        bin_size: int = 1,
        slide_size: float = 1
        ) -> pd.DataFrame:
    """
    Bins the mutations relative to recombination events into specified regions with sliding bins.

    Parameters:
    df (pd.DataFrame): DataFrame containing information about the mutations and their relative positions to recombination events.
    region (tuple): Tuple specifying the lower and upper bounds of the region. Default is (-12, 12).
    bin_size (int): Size of each bin in kb. Default is 1kb.
    slide_size (float): Size of the slide for each bin in kb. Default is 0.1kb.

    Returns:
    pd.DataFrame: DataFrame with the mutations binned into the specified regions.
    """
    bottom_region = region[0]
    top_region = region[1]
    bin_width = bin_size
    slide_width = slide_size
    print(f"Analysis region: {bottom_region} to {top_region} kb, bin width: {bin_width} kb, slide width: {slide_width} kb")

    df_trimmed = constrain_to_region(df, region)

    df_trimmed = df_trimmed[["genotype", "Region", "Reference", "norm_count"]]

    bins = []
    current_start = bottom_region
    while current_start + bin_width <= top_region:
        bins.append((current_start, current_start + bin_width))
        current_start += slide_width

    bin_labels = range(len(bins))
    df_trimmed["Region"] = pd.cut(df_trimmed["Region"], bins=[b[0] for b in bins] + [bins[-1][1]], labels=bin_labels, include_lowest=True)

    # Recalculate the Region values to be relative to the recombining end
    df_trimmed["Region"] = df_trimmed["Region"].apply(lambda x: bins[x][0] + bin_width / 2 if pd.notna(x) else x)

    return df_trimmed

def normalize_rec_relative_bins_to_wt(df, normalize_to="ung1∆"):
    """
    Normalizes the counts of mutations in each bin to a specified genotype.

    Parameters:
    df (pd.DataFrame): DataFrame containing information about the mutations, their genotypes, and their relative positions to recombination events.
    normalize_to (str): Genotype to normalize the counts to. Default is "ung1∆".

    Returns:
    pd.DataFrame: DataFrame with the counts of mutations in each bin normalized to the specified genotype.
    """
    df_counts = df.groupby(["genotype", "Region"]).sum() #sum up the norm_count values for each Region bin
    df_counts = df_counts.unstack(level=0) #transform the counts into a DataFrame, with the bins as columns
    df_counts = df_counts.fillna(0)

    df_norm = df_counts.div(df_counts["norm_count"][normalize_to], axis=0)
    df_norm.columns = df_norm.columns.droplevel()
    df_norm = df_norm.reset_index()
    df_norm = df_norm.rename(columns={"Region": "Region_bin"})

    return df_norm

def conduct_wilcoxon_signed_rank_test(
        meta_master_relative_GC_CO_df: pd.DataFrame,
        region_flank: int | float = 8,
        bins: int = 64
        ) -> None:
    """
    Conducts a Wilcoxon signed-rank test on the normalized counts of mutations in each bin.

    Parameters:
    meta_master_relative_GC_CO_df (pd.DataFrame): DataFrame containing information about the mutations and their relative positions to recombination events for a meta-analysis.
    region_flank (int | float): Number specifying the lower and upper bounds of the region. Default is 8.
    bins (int): Number of bins to divide the region into. Default is 64.

    Returns:
    None: The function prints the results of the Wilcoxon signed-rank test and does not return anything.
    """
    from scipy.stats import wilcoxon

    region = (-region_flank, region_flank)

    df_trimmed = bin_rec_relative_mutations(meta_master_relative_GC_CO_df, region, bins=bins)
    df_norm = normalize_rec_relative_bins_to_wt(df_trimmed, normalize_to="ung1∆")
    print("Dataframe with normalized counts to ung1∆ for wilcoxon signed rank test:")
    print(df_norm)

    for focus_genotype in df_norm.columns[1:].tolist():

        #conduct a wilcoxon singed rank test, between Region values that are smaller than 0 and those that are larger than 0
        #make sure that region_bin is a float type
        df_norm["Region_bin"] = df_norm["Region_bin"].astype(float)
        df_smaller_than_zero = df_norm[df_norm["Region_bin"] < 0]
        df_larger_than_zero = df_norm[df_norm["Region_bin"] > 0]
        
        try:
            stat, p = wilcoxon(df_smaller_than_zero[focus_genotype], df_larger_than_zero[focus_genotype])
            print(f"Wilcoxon signed rank test for {focus_genotype} left and right of 0: {p}, stat: {stat}")

        except ValueError as e:
            print(f"Wilcoxon signed rank test needs equal number of samples in both groups and they difference needs to not be all 0s. Skipping {focus_genotype}.")
            print(e)
      
def test_third_order_poly_fit_F_test(df_subset):
    """
    Conducts an F-test to compare the fit of a third order polynomial to a linear fit on a subset of data.

    Parameters:
    df_subset (pd.DataFrame): DataFrame containing the data to be fitted. It should have columns "Region_bin" and "norm_count".

    Returns:
    float: The p-value from the F-test. A small p-value (typically ≤ 0.05) indicates strong evidence that the third order polynomial is a better fit to the data than a linear fit.
    """
    import scipy.stats as stats
    import numpy as np

    coefficients = np.polyfit(df_subset["Region_bin"], df_subset["norm_count"], 3) # Fit a 3rd degree polynomial
    predicted = np.polyval(coefficients, df_subset["Region_bin"]) # Get the predicted values
    residuals = df_subset["norm_count"] - predicted # Calculate the residuals
    dof_residuals = len(df_subset["Region_bin"]) - 4 # Degrees of freedom of the residuals

    # Calculate the residuals for a linear hypothesis
    residuals_linear = df_subset["norm_count"] - np.polyval(np.polyfit(df_subset["Region_bin"], df_subset["norm_count"], 1), df_subset["Region_bin"])

    # Degrees of freedom of the residuals under the linear hypothesis
    dof_linear = len(df_subset["Region_bin"]) - 2

    # Perform the F-test
    f_stat = (np.sum(residuals_linear**2) - np.sum(residuals**2)) / (dof_linear - dof_residuals) / (np.sum(residuals**2) / dof_residuals)
    p_value = stats.f.sf(f_stat, dof_linear - dof_residuals, dof_residuals)

    return p_value

def get_GC_CO_numbers(
    samples_SNVs: List[OverlapSwitchesTable],
    query_genotypes: str | None = None,
    min_percent: int | float = 80,
    include_aneuploid_heterozygous: bool = False,
    ) -> tuple:
    """
    This function calculates and returns the GC and CO numbers for each haploid and homozygous sample in the provided list of genotypes.
    
    Parameters:
    samples_SNVs (List[OverlapSwitchesTable]): A list of sample objects to calculate GC and CO numbers for.
    query_genotypes (str | None): A string representing the genotypes to query. If None, all genotypes are considered. Default is "ung1∆+ung1∆NAT+exo1-nd+pol32∆+exo1-ndpol32∆".
    min_percent (int | float): The minimum percentage of the Parental SNPs detected for a sample to be considered. Default is 80.
    
    Returns:
    tuple: A tuple containing two lists. The first list contains the GC numbers for each sample, and the second list contains the CO numbers for each sample.
    
    Note:
    If a sample has fewer than `min_percent` of the genome covered by SNVs, a warning is logged and the sample is skipped.
    If a sample has fewer than 10 GCs or COs, a message is printed to the console.
    """
    from mayo.settings import config as cfg
    from mayo.settings import genotype_dict as genotype_dict_master
    genotype_dict = genotype_dict_master[cfg["genotype_dict"]]

    GC_numbers = []
    CO_numbers = []

    used_samples = 0

    for sample in samples_SNVs:

        if query_genotypes is not None and sample.name not in get_samples_with_genotype(query_genotypes, genotype_dict):
            continue

        if sample.percent <= min_percent:
            print(f"{sample.name} sample has fewer than {min_percent}% of the detected parental SNPs. GC and CO numbers can be inaccurate. Skipping in summary.")
            logging.warning(f"{sample.name} sample has fewer than 80% of the detected parental SNPs. GC and CO numbers can be inaccurate. Skipping in summary.")
            continue

        if not include_aneuploid_heterozygous:
            if sample.aneuploid == True or sample.heterozygous == True:
                logging.info(f"{sample.name} sample is aneuploid or heterozygous. Skipping in summary.")
                continue
        
        #check if total_GC and total_CO are already calculated
        if not hasattr(sample, "total_CO") and hasattr(sample, "total_CO"):
            print(f"Sample {sample.name} has no total_GC and total_CO attributes. Calculating... for entire dataset")
            calculate_GC_CO(samples_SNVs)
            
        total_GC = sample.total_GC
        total_CO = sample.total_CO

        GC_numbers.append(total_GC)
        CO_numbers.append(total_CO)

        used_samples += 1

        logging.info(f"{sample.name} sample has a total of {total_GC:3} GCs and {total_CO:3} COs")

        if total_GC < 10:
            print(f"{sample.name} sample has fewer than 10 detected GCs... Please verify.")
        if total_CO < 10:
            print(f"{sample.name} sample has fewer than 10 detected COs... Please verify.")
    
    print(f"Used {used_samples} samples for GC and CO numbers.")

    return GC_numbers, CO_numbers

def get_mutations_around_rec_ev(
        samples_SNVs: list[OverlapSwitchesTable], 
        genotype_dict: dict, 
        focus_genotype: str | None = None, 
        mutation_type: Literal["clustered", "scattered", "all"] = "all",
        homozygous_only: bool = False,
        event_type: Literal["GC", "CO", "both"] = "both",
        window_size: int = 10000
        ) -> pd.DataFrame:
    """
    Retrieves mutations around recombination events for a list of samples.

    Parameters:
    samples_SNVs (list): List of samples with SNVs.
    genotype_dict (dict): Dictionary mapping genotypes to samples.
    focus_genotype (str | None): Specific genotype to focus on. If None, all genotypes are considered. Default is None.
    mutation_type (str): Type of mutations to consider. Can be "clustered", "scattered", or "all". Default is "all".
    homozygous_only (bool): If True, only considers homozygous mutations. Default is False.
    event_type (Literal["GC", "CO", "both"]): Type of recombination event to consider. Can be "GC" (Gene Conversion), "CO" (Crossover), or "both". Default is "both".
    window_size (int): Size of the window around the recombination event to consider. Default is 10000.

    Returns:
    pd.DataFrame: DataFrame containing information about the mutations around recombination events.
    """
    from mayo.settings import rtg_list
    genotype_dict_inv = invert_dict(genotype_dict)

    master_relative_GC_df = pd.DataFrame()
    master_relative_CO_df = pd.DataFrame()

    for sample in samples_SNVs:
        #print(f"Processing {sample.name}...")

        GC_dict = sample.GC_dict.copy()
        CO_dict = sample.CO_dict.copy()

        #new, skip the sample if it's in the rtg_list
        if sample.name in rtg_list:
            print(f"Skipping {sample.name} as it's in the rtg_list.")
            logging.info(f"Skipping {sample.name} as it's in the rtg_list.")
            continue

        if focus_genotype is not None and sample.name not in sum([genotype_dict_inv[i] for i in focus_genotype.split("+")], []):
            continue

        if homozygous_only:
            if sample.heterozygous is True:
                
                #remove heterozygous chromosome keys from GC_dict and CO_dict
                GC_dict = {k: v for k, v in GC_dict.items() if k not in sample.heterozygous_chromosomes}
                CO_dict = {k: v for k, v in CO_dict.items() if k not in sample.heterozygous_chromosomes}
        
        if mutation_type == "clustered":
            assert sample.df_SNPs[sample.df_SNPs["Is_Clustered"] == False].empty, "There aren't clustered mutations in the sample's df_SNPs DataFrame. Please address this before running this script."
            relative_mut_GC_df = pull_mutations_from_sample_recev_dict(GC_dict, sample.df_SNPs, window_size)
            relative_mut_CO_df = pull_mutations_from_sample_recev_dict(CO_dict, sample.df_SNPs, window_size)

        elif mutation_type == "scattered":
            relative_mut_GC_df = pull_mutations_from_sample_recev_dict(GC_dict, sample.df_scattered, window_size)
            relative_mut_CO_df = pull_mutations_from_sample_recev_dict(CO_dict, sample.df_scattered, window_size)

        elif mutation_type == "all":
            relative_mut_GC_df = pull_mutations_from_sample_recev_dict(GC_dict, sample.df, window_size)
            relative_mut_CO_df = pull_mutations_from_sample_recev_dict(CO_dict, sample.df, window_size)

        #add sample_name information
        relative_mut_GC_df["Sample"] = sample.name
        relative_mut_CO_df["Sample"] = sample.name            
        
        #here could normalize by the number of A3A motifs available?

        master_relative_GC_df = pd.concat([master_relative_GC_df, relative_mut_GC_df], ignore_index=True)
        master_relative_CO_df = pd.concat([master_relative_CO_df, relative_mut_CO_df], ignore_index=True)
        
        total_relative_mutations = len(master_relative_GC_df) + len(master_relative_CO_df)

    if len(master_relative_GC_df) > 0:
        #master_relative_GC_df["norm_count"] = ((1 / sample.total_GC) / len(master_relative_GC_df)) * 1000
        #master_relative_GC_df["norm_count"] = ((1 / len(master_relative_GC_df)) * 100)
        master_relative_GC_df["norm_count"] = ((1 / total_relative_mutations) * 100)
    else:
        master_relative_GC_df["norm_count"] = 0
        
    if len(master_relative_CO_df) > 0:
        #master_relative_CO_df["norm_count"] = ((1 / sample.total_CO) / len(master_relative_CO_df)) * 1000
        #master_relative_CO_df["norm_count"] = ((1 / len(master_relative_CO_df)) * 100)
        master_relative_CO_df["norm_count"] = ((1 / total_relative_mutations) * 100)
    else:
        master_relative_CO_df["norm_count"] = 0

    master_relative_GC_df["event"] = "GC"
    master_relative_CO_df["event"] = "CO"

    master_relative_GC_CO_df = pd.concat([master_relative_GC_df, master_relative_CO_df], ignore_index=True)
    master_relative_GC_CO_df = master_relative_GC_CO_df[~master_relative_GC_CO_df["Reference"].isin(["A", "T"])]
    master_relative_GC_CO_df = master_relative_GC_CO_df[~master_relative_GC_CO_df["Allele"].isin(["C", "G"])]
    master_relative_GC_CO_df = master_relative_GC_CO_df.drop(columns=["Chromosome", "Allele", "Coverage", "Frequency"])

    master_relative_GC_CO_df["Reference"] = master_relative_GC_CO_df["Reference"].replace(["C"], "C->T")
    master_relative_GC_CO_df["Reference"] = master_relative_GC_CO_df["Reference"].replace(["G"], "G->A")

    master_relative_GC_CO_df["Region"] = master_relative_GC_CO_df["Region"] / 1000 #normalize x to 1kb

    if event_type != "both":
        master_relative_GC_CO_df = master_relative_GC_CO_df[master_relative_GC_CO_df["event"] == event_type]

    print(master_relative_GC_CO_df["Reference"].value_counts())
    print(master_relative_GC_CO_df["event"].value_counts())

    return master_relative_GC_CO_df

def plot_SNVs_around_rec_ev_kde_meta(
        meta_master_relative_GC_CO_df, 
        title: str = "",
        save_path: str | None = None,
        show: bool = True,
        mutation_type: Literal["clustered", "scattered", "all"] | None = None
        ):
    import numpy as np
    """
    Plots the distribution of Single Nucleotide Variants (SNVs) around recombination events for a meta-analysis.

    Parameters:
    meta_master_relative_GC_CO_df (pd.DataFrame): DataFrame containing information about the SNVs and their relative positions to recombination events for a meta-analysis.
    title (str): Title for the plot. Default is an empty string.
    save_path (str | None): Path to save the plot. If None, the plot is not saved. Default is None.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is True.
    mutation_type (str): Type of mutations to consider. Can be "clustered", "scattered", "all", or None. Default is None.

    Returns:
    None
    """
    sns.set_context("talk")
    plt.figure(figsize=(16, 10))
    sns.set_style("whitegrid")
    #sns.histplot(data=meta_master_relative_GC_CO_df, x="Region", hue="genotype", kde=True, bins=100, kde_kws={'bw_adjust': .3})
    plot = sns.kdeplot(
        data=meta_master_relative_GC_CO_df, 
        x="Region", 
        hue="genotype", 
        cut=0, 
        fill=False, 
        common_norm=False, 
        common_grid=True, 
        bw_adjust=.5)
    plt.axvline(x=0, color="black", linestyle="--", linewidth=1)
    #cut the plot at -12 and 12
    plt.xlim(-12, 12)

    if mutation_type is not None:
        plt.title(f"{mutation_type.capitalize()} mutations relative to recombination events ({title})", fontweight="bold", fontsize=30)
    else:
        plt.title(f"Mutations relative to recombination events ({title})", fontweight="bold", fontsize=30)

    plt.xticks(np.arange(-12, 13, 1))

    from matplotlib import lines as mlines
    # Extract unique genotypes and their corresponding colors
    genotypes = meta_master_relative_GC_CO_df["genotype"].unique()
    palette = sns.color_palette(n_colors=len(genotypes))

    # Create custom legend handles using lines
    handles = []
    labels = []

    for genotype, color in zip(genotypes, palette):
        handle = mlines.Line2D([], [], color=color, label=genotype, linestyle='-', linewidth=2)
        handles.append(handle)
        labels.append(genotype)

    # Add custom legend with italic font style
    plt.legend(handles=handles, labels=labels, title='', prop={'style': 'italic'})

    # Add custom legend with italic font style
    plt.legend(handles=handles, labels=labels, title='', prop={'style': 'italic'})

    plt.xlabel("Distance from event (kb)", fontweight="bold")
    plt.ylabel("Frequency", fontweight="bold")

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

def get_rec_ev_relative_A3A_SNVs_and_plot(
        samples_SNVs, 
        genotype_dict, 
        genotype_list, 
        mutation_type: Literal["clustered", "scattered", "all"] = "all",
        window_size = 8000,
        event_type: Literal["GC", "CO", "both"] = "both"
        ) -> pd.DataFrame:
    """
    This function calculates mutations around recombination events for each genotype in the genotype list, plots the results, and returns a DataFrame with the results for all genotypes.

    Parameters:
    samples_SNVs (DataFrame): A DataFrame containing information about the SNVs in the samples.
    genotype_dict (dict): A dictionary mapping genotypes to their corresponding samples.
    genotype_list (list): A list of genotypes to consider.
    mutation_type (str): Type of mutations to consider. Can be "clustered", "scattered", or "all". Default is "all".

    Returns:
    DataFrame: A DataFrame containing information about the mutations around recombination events for all genotypes in the genotype list.
    """
    meta_master_relative_GC_CO_df = pd.DataFrame()

    for focus_genotype in genotype_list:

        print(f"Calculating mutations around recombination events for {focus_genotype} samples...")

        master_relative_GC_CO_df = get_mutations_around_rec_ev(
            samples_SNVs,
            genotype_dict,
            focus_genotype=focus_genotype, #pol32∆ None ung1∆+ung1∆NAT+exo1-nd+pol32∆+exo1-ndpol32∆
            mutation_type=mutation_type,
            homozygous_only=True,
            event_type=event_type,
            window_size=window_size)
        
        plot_SNVs_around_rec_ev(
            master_relative_GC_CO_df, 
            title=focus_genotype, 
            save_path=f"figures/mutations_around_rec_events/mutations_around_rec_events_{focus_genotype}_{mutation_type}.png", 
            mutation_type=mutation_type,
            show=False)
        
        #if Reference == G->A, then multiply the Region by -1 and change the Reference to C->T
        master_relative_GC_CO_df.loc[master_relative_GC_CO_df["Reference"] == "G->A", "Region"] = master_relative_GC_CO_df["Region"] * -1
        master_relative_GC_CO_df["Reference"] = master_relative_GC_CO_df["Reference"].replace(["G->A"], "C->T")
        master_relative_GC_CO_df["Reference"] = master_relative_GC_CO_df["Reference"].replace(["C->T"], "C->T or G->A")
        master_relative_GC_CO_df["genotype"] = focus_genotype
        meta_master_relative_GC_CO_df = pd.concat([meta_master_relative_GC_CO_df, master_relative_GC_CO_df], ignore_index=True)

        plot_SNVs_around_rec_ev(
            master_relative_GC_CO_df, 
            title=focus_genotype, 
            save_path=f"figures/mutations_around_rec_events/mutations_around_rec_events_{focus_genotype}_flipped_{mutation_type}.png",
            mutation_type=mutation_type,
            show=False)

    plot_SNVs_around_rec_ev_kde_meta(
        meta_master_relative_GC_CO_df, 
        title="combined", 
        save_path=f"figures/mutations_around_rec_events/mutations_around_rec_events_combined_flipped_kde_{mutation_type}.png", 
        mutation_type=mutation_type,
        show=False)

    return meta_master_relative_GC_CO_df

def plot_GC_CO_distribution_by_genotype(
        samples_SNVs, 
        query_genotypes: str | None = None,
        save_name: str | None = "distribution_haploid",
        show: bool = True
        ) -> None:
    """
    This function plots the distribution of GCs (Gene Conversions) and COs (Crossover events) for the given genotypes. 
    It also calculates and prints the average, median, and standard deviation of GCs and COs.

    Parameters:
    samples_SNVs (DataFrame): A DataFrame containing information about the SNVs in the samples.
    query_genotypes (str | None): A string representing the genotypes to consider. If None, all genotypes are considered.
    save_name (str | None): The name to use when saving the plot. If None, the plot is not saved.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed.

    Returns:
    None
    """    
    import statistics

    GC_numbers, CO_numbers = get_GC_CO_numbers(samples_SNVs, query_genotypes=query_genotypes, min_percent=80)

    #combine the GC and CO numbers into one dataframe for plotting
    GC_CO_df = pd.DataFrame({"Event": ["GC"] * len(GC_numbers) + ["CO"] * len(CO_numbers), "Number": GC_numbers + CO_numbers})

    print(f"There are a total of {sum(GC_numbers)} GCs and {sum(CO_numbers)} COs in queried samples")

    for data, name in [
        (GC_numbers, "GCs"),
        (CO_numbers, "COs"),
        (GC_CO_df, "GCs and COs")]:
        
        #draw a distribution of GCs and COs as a violin plot with overlaid swarmplot
        sns.set_style('whitegrid')
        sns.set_context("notebook", font_scale=3)
        plt.figure(figsize=(6, 6))

        if name == "GCs and COs":
            sns.violinplot(x="Event", y="Number", data=data, color="whitesmoke", inner="quartile", density_norm="width", linewidth=3)
            sns.swarmplot(x="Event", y="Number", data=data, alpha=0.8, size=7, edgecolor="black")
        else:
            sns.violinplot(y=data, color="whitesmoke", inner="quartile", density_norm="width", linewidth=3)
            sns.swarmplot(y=data, alpha=0.9, size=8, edgecolor="black")

        plt.title(f"Distribution of {name}", fontsize=20, fontweight="bold")
        ymin, ymax = plt.gca().get_ylim()
        y_span = ymax - ymin
        y_span_percent_pad = y_span * 0.03
        plt.gca().set_ylim(bottom=0 - y_span_percent_pad)
        plt.gca().yaxis.set_major_locator(plt.MultipleLocator(10)) #make y ticks every 10
        plt.xticks(fontsize=18, fontweight="bold")
        plt.yticks(fontsize=18, fontweight="bold")
        plt.ylabel("Number of Events", fontsize=20, fontweight="bold")
        plt.xlabel("Event Type", fontsize=20, fontweight="bold")

        if save_name is not None:
            plt.savefig(f"figures/{name}_{save_name}_{query_genotypes}.png", dpi=300, bbox_inches='tight')
        
        if show:
            plt.show()
        else:
            plt.close()

    print(f"Average GCs: {sum(GC_numbers) / len(GC_numbers)}")
    print(f"Median GCs: {statistics.median(GC_numbers)}")
    print(f"Standard deviation GCs: {statistics.stdev(GC_numbers)}")

    print(f"Average COs: {sum(CO_numbers) / len(CO_numbers)}")
    print(f"Median COs: {statistics.median(CO_numbers)}")
    print(f"Standard deviation COs: {statistics.stdev(CO_numbers)}")

def find_chromosome_length(chromosome_name: str, chromosomes) -> int:
    """
    This function finds and returns the length of a specified chromosome.

    Parameters:
    chromosome_name (str): The name of the chromosome for which to find the length.
    chromosomes (list of tuples): A list of tuples, where each tuple contains a chromosome name and its end position, representing its length.

    Returns:
    int: The length of the specified chromosome.
    """
    df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])
    chromosome_length = df_chr[df_chr["chromosome"] == chromosome_name]["end_position"].values[0]
    return chromosome_length

def get_all_rec_events_df(
        samples_SNVs: list[OverlapSwitchesTable],
        genotype_dict: dict,
        query_genotypes: str = "ung1∆+ung1∆NAT+exo1-nd+pol32∆+exo1-ndpol32∆",
        skip_aneuploid_heterozygous_sample: bool = False,
        skip_aneuploid_heterozygous_chromosome: bool = True
        ) -> pd.DataFrame:
    """
    This function generates a DataFrame containing all recombination events (GCs and COs) for the specified genotypes.

    Parameters:
    samples_SNVs (list[OverlapSwitchesTable]): A list of OverlapSwitchesTable objects, each representing a sample.
    query_genotypes (str): A string representing the genotypes to consider. Default is "ung1∆+ung1∆NAT+exo1-nd+pol32∆+exo1-ndpol32∆".
    genotype_dict (dict): A dictionary mapping genotypes to samples.
    skip_aneuploid_heterozygous (bool): If True, aneuploid and heterozygous samples are skipped. If False, haploid chromosomes among
    aneuploid and heterozygous samples are considered in addition. Default is True.

    Returns:
    pd.DataFrame: A DataFrame containing information about all recombination events for the specified genotypes. The DataFrame has columns for Chromosome, Region, Type of event (GC or CO), and Sample name.
    """
    df_all_rec_events = pd.DataFrame()

    queried_samples = get_samples_with_genotype(query_genotypes, genotype_dict)

    pulled_samples_counter = 0

    for sample in samples_SNVs:

        if sample.name in queried_samples:
            if skip_aneuploid_heterozygous_sample:
                if sample.aneuploid == True or sample.heterozygous == True:
                    #print(f"Skipping {sample.name} because it is aneuploid or heterozygous")
                    continue
            
            pulled_samples_counter += 1
            GC_CO_positions = []

            #get all positions of GCs
            for chr in sample.GC_dict.keys():

                #check if the chromosome is aneuploid or heterozygous
                if skip_aneuploid_heterozygous_chromosome:
                    if sample.aneuploid == True or sample.heterozygous == True:
                        if chr in sample.heterozygous_chromosomes:
                            continue
                        if chr in sample.aneuploid_chromosomes:
                            continue

                positions = sample.GC_dict[chr]

                for position in positions:
                    GC_CO_positions.append([chr, position, "GC"])

            #get all positions of COs
            for chr in sample.CO_dict.keys():
                positions = sample.CO_dict[chr]

                for position in positions:
                    GC_CO_positions.append([chr, position, "CO"])

            df_all_rec_events_sample = pd.DataFrame(GC_CO_positions, columns=["Chromosome", "Region", "Type"])
            df_all_rec_events_sample["Sample"] = sample.name
            df_all_rec_events = pd.concat([df_all_rec_events, df_all_rec_events_sample], ignore_index=True)

    print(f"Pulled up from {pulled_samples_counter} samples")

    return df_all_rec_events

def smoothen_GC_CO_data(
        df: pd.DataFrame, 
        smoothing_range: tuple,
        event_type: Literal["GC", "CO"],
        window_size: int = 5000, 
        window_slide: int = 500,
        ) -> pd.DataFrame:
    """
    This function applies a sliding window approach to smoothen the data of GC (Gene Conversion) or CO (Crossover) events. 
    It counts the number of events in each window and shifts the data to the center of the window.

    Parameters:
    df (pd.DataFrame): A DataFrame containing information about the events.
    smoothing_range (tuple): A tuple representing the range over which to apply the smoothing.
    event_type (Literal["GC", "CO"]): The type of event to consider. Either "GC" for Gene Conversion or "CO" for Crossover.
    window_size (int): The size of the sliding window. Default is 5000.
    window_slide (int): The step size for the sliding window. Default is 500.

    Returns:
    pd.DataFrame: A DataFrame containing the smoothened data. The DataFrame has columns for Region, Events, and Type.
    """
    df_sliding_window = pd.DataFrame(columns=["Region", "Type", "Sample"])

    for i in range(smoothing_range[0], smoothing_range[1], window_slide):

        df_window = df[(df["Region"] >= i) & (df["Region"] <= i+window_size)].copy()
        total_events = len(df_window) #count the number of events in the window
        #add the window information to the df_sliding_window (#of events, position)
        df_sliding_window = pd.concat([df_sliding_window, pd.DataFrame([[i, total_events, event_type]], columns=["Region", "Events", "Type"])], ignore_index=True)
    
    #window_size smoothing shifts the data to the left, so we need to shift it back to the right
    df_sliding_window["Region"] = df_sliding_window["Region"] + (window_size / 2)

    return df_sliding_window

def plot_rec_events_for_all_chr(
        df_all_rec_events: pd.DataFrame, 
        SNP_NUM: int,
        event_types: list = ["GC", "CO"],
        save: bool = True,
        show: bool = False,
        lineplot_kde: bool = False
        ) -> None:
    """
    This function plots the distribution of recombination events (GCs and COs) for all chromosomes. 
    It applies a sliding window approach to smoothen the data and plots both the histogram of the events and the smoothened data.

    Parameters:
    df_all_rec_events (pd.DataFrame): A DataFrame containing information about all recombination events.
    SNP_NUM (int): The number of SNPs to consider.
    event_types (list): A list of event types to consider. Default is ["GC", "CO"].
    save (bool): If True, the plot is saved. Default is True.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is False.

    Returns:
    None
    """
    from mayo.settings import S288C_to_roman
    from mayo.settings import chromosomes as chromosomes
    import numpy as np

    #plot the distribution of COs on the chromosomes as a histogram
    for chromosome_to_plot in df_all_rec_events["Chromosome"].unique():

        if chromosome_to_plot == "ref|NC_001224|":
            continue

        df_all_rec_events_chr = df_all_rec_events[df_all_rec_events["Chromosome"] == chromosome_to_plot]
        df_both_slide = pd.DataFrame()

        for event_type in event_types:
            
            df_all_rec_events_chr_event = df_all_rec_events_chr[df_all_rec_events_chr["Type"] == event_type]

            print(f"total_{event_type}s on {chromosome_to_plot}: {len(df_all_rec_events_chr_event)}")
            print(f"Most common {event_type}s on {chromosome_to_plot}:")
            print(df_all_rec_events_chr_event["Region"].value_counts().head(5))
            print()

            chromosome_length = find_chromosome_length(chromosome_to_plot, chromosomes)
            df_sliding_window = smoothen_GC_CO_data(
                df=df_all_rec_events_chr_event, 
                smoothing_range=(0, chromosome_length), 
                event_type=event_type,
                window_size = 5000, 
                window_slide = 500)
            
            df_sliding_window = smoothen_derive_col_data(df_sliding_window, "Events", window_length=12)
            df_sliding_window.loc[df_sliding_window["Events_smoothed"] < 0, "Events_smoothed"] = 0 #if Events_smoothed < 0, set it to 0
            df_sliding_window["Type"] = event_type
            df_both_slide = pd.concat([df_both_slide, df_sliding_window], ignore_index=True)


        #plot a lineplot of the sliding window and the histogram of the GC and CO events
        chromosome_roman = S288C_to_roman[chromosome_to_plot]
        palette={"GC": "tab:blue", "CO": "tab:red"}
        fig, ax = plt.subplots(figsize=(20, 5))
        #create a stacked histogram of the GC and CO events
        hist = sns.histplot(ax=ax, data=df_all_rec_events_chr, x="Region", hue="Type", stat="count", bins=np.arange(0, chromosome_length, 2000), element="bars", palette=palette, multiple="stack", legend=True)
        hist.legend_.set_title("Event type")
        hist.legend_.get_title().set_fontsize('12')
        hist.legend_.get_title().set_fontweight('bold')
        for t in hist.legend_.texts:
            t.set_fontsize('12')

        if lineplot_kde:
            sns.lineplot(ax=ax, data=df_both_slide, x="Region", y="Events_smoothed", hue="Type", palette=palette, legend=False, alpha=0.8)
        plt.title(f"Recombination Events on Chromosome {chromosome_roman}", fontsize=24, fontweight="bold")
        plt.xlim(0, chromosome_length)
        plt.xlabel("Chromosome position", fontsize=20, fontweight="bold")
        plt.ylabel("Number of events", fontsize=20, fontweight="bold")
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tight_layout()

        if save:
            plt.savefig(f"figures/rec_events/rec_events_on_chr_{chromosome_roman}_{SNP_NUM}SNP.png", dpi=300)
        if show:
            plt.show()
        else:
            plt.close()
            
def plot_clust_scattered_ratio_violin(
        df_ratios: pd.DataFrame, 
        save_path: str | None = None,
        show: bool = False
        ) -> None:
    """
    This function plots the distribution of clustered/scattered mutation ratios for different p-values as a violin plot.
    Parameters:
    df_ratios (pd.DataFrame): A DataFrame containing the ratios of clustered to scattered mutations.
    """
    total_samples = len(df_ratios)
    plt.figure(figsize=(10, 6))
    sns.set_context("paper")
    sns.set_style("whitegrid")
    sns.violinplot(data=df_ratios, color="whitesmoke", inner="quartile", scale="width", edgecolor="gray", linewidth=2)
    sns.swarmplot(data=df_ratios, size=5, edgecolor=None, alpha=0.85)

    #get the max of the violinplot
    violin_max = df_ratios.max().max()
    plt.ylim(-0, violin_max + 8) #cut the plot at 0
    plt.title(f"% Clustered mutations for {total_samples} samples", fontsize=18, fontweight="bold")
    plt.xlabel("Cluster p-value", fontsize=14, fontweight="bold")
    plt.ylabel("% Clustered mutations", fontsize=14, fontweight="bold")
    plt.xticks(fontsize=12, fontweight="bold")
    plt.yticks(fontsize=12, fontweight="bold")

    if save_path is not None:
        plt.savefig(f"{save_path}", dpi=300, bbox_inches="tight")
        print(f"Plot saved as {save_path}")

    if show:
        plt.show()
    else:
        plt.close()

def plot_clust_scattered_ratio_kde(
        df_ratios: pd.DataFrame, 
        save_path: str | None = None,
        show: bool = False
        ) -> None:
    """
    This function plots the distribution of clustered/scattered mutation ratios for different p-values as a KDE (Kernel Density Estimation) line plot.

    Parameters:
    df_ratios (pd.DataFrame): A DataFrame containing the ratios of clustered to scattered mutations.
    save (bool): If True, the plot is saved. Default is True.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is False.

    Returns:
    None
    """
    total_samples = len(df_ratios)
    plt.figure(figsize=(10, 6))
    sns.set_context("poster")
    sns.set_style("whitegrid")
    sns.kdeplot(data=df_ratios)
    #cut plot at 0 and go to max of the violinplot
    plt.xlim(left=0)
    plt.title(f"% Clustered mutations \nfor {total_samples} ung1∆ samples by cluster p-value", fontsize=20, fontweight="bold")
    plt.xlabel("% Clustered mutations", fontsize=20, fontweight="bold")
    plt.ylabel("Frequency", fontsize=20, fontweight="bold")
    plt.xticks(fontsize=16, fontweight="bold")
    plt.yticks(fontsize=16, fontweight="bold")

    if save_path is not None:
        plt.savefig(f"{save_path}", dpi=300, bbox_inches="tight")
        print(f"Plot saved as {save_path}")
    
    if show:
        plt.show()
    else:
        plt.close()

def calculate_scattered_SNPs(samples_SNVs: list[OverlapSwitchesTable]) -> None:
    """
    Calculates the number of scattered SNVs for each sample.

    Parameters:
    samples_SNVs (list): List of samples with SNVs.

    Returns:
    None
    """
    for sample in samples_SNVs:
        df_all_snvs = sample.df.copy()
        df_clust_snvs = sample.df_SNPs.copy()

        total_clust_SNVs = len(df_clust_snvs)
        expected_scattered_SNVs = len(df_all_snvs) - total_clust_SNVs

        #remove all SNVs that are not in clusters from df_all_snvs by using both Chromosome and Region columns as a key
        df_all_snvs = df_all_snvs.set_index(["Chromosome", "Region"])
        df_clust_snvs = df_clust_snvs.set_index(["Chromosome", "Region"])
        df_scattered_snvs = df_all_snvs.drop(df_clust_snvs.index)
        df_scattered_snvs = df_scattered_snvs.reset_index()
        assert len(df_scattered_snvs) == expected_scattered_SNVs, f"df_scattered_snvs length: {len(df_scattered_snvs)}, expected_scattered_SNVs: {expected_scattered_SNVs}"

        sample.df_scattered = df_scattered_snvs

def expand_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands a DataFrame based on the interval length and value in each row.

    This function takes a DataFrame with columns 'position', 'interval_length', and 'Value'. 
    It creates a new DataFrame where for each row in the original DataFrame, 
    it creates 'interval_length' number of new rows in the new DataFrame. 
    The 'position' for these new rows starts from the 'position' in the original row and 
    increments by 1 for each new row. The 'Value' for these new rows is the same as the 'Value' 
    in the original row.

    After expanding the DataFrame, it adjusts the 'Position' column to be 1-indexed 
    (as opposed to 0-indexed) to match the 1-based indexing of the wig file format.

    Parameters:
    df (pandas.DataFrame): The original DataFrame to be expanded. 
                           It should have columns 'position', 'interval_length', and 'Value'.

    Returns:
    pandas.DataFrame: The expanded DataFrame with columns 'Position' and 'Value'.
    """
    positions = []
    values = []

    for _, row in df.iterrows():
        #chromosomes.extend([row.Chromosome] * row.interval_length)
        positions.extend(range(row.position, row.position + row.interval_length))
        values.extend([row.Value] * row.interval_length)

    df_new = pd.DataFrame({
        'Position': positions,
        'Value': values
    })

    #move the position by 1 (add 1 to each position) to match the 1-based indexing of the wig file format
    df_new['Position'] = df_new['Position'] + 1
    df_new['Position'] = df_new['Position'].astype(int)

    return df_new

def process_by_chromosome_bed(file_path: str) -> pd.DataFrame:
    """
    Processes a BED file by chromosome and converts it to a DataFrame.

    This function reads a BED file, which is a tab-delimited text file that defines data lines.
    Each line contains a chromosome, start, end, and value. The function processes the file by 
    chromosome, creating a new DataFrame with a single position column instead of start and end columns.
    It inserts new rows for each position between start and end and fills in the Value column with the same value.
    The function also renames the chromosome name to the S288C format using the roman_to_S288C dictionary.

    Parameters:
    file_path (str): The path to the BED file to be processed.

    Returns:
    pandas.DataFrame: The processed DataFrame with columns 'Chromosome', 'Position', and 'Value'.
    """
    #skip the first row, which is a header
    df = pd.read_csv(file_path, sep='\t', header=None, skiprows=1, names=['Chromosome', 'Start', 'End', 'Value'])
    df_out = pd.DataFrame(columns=['Chromosome', 'Position', 'Value'])

    for chromosome in df['Chromosome'].unique():
        df_chromosome = df[df['Chromosome'] == chromosome].copy()

        #instead of having start and end, we want to have a single position column
        #so we will insert new rows for each position between start and end and fill in the Value column with the same value
        df_chromosome['interval_length'] = df_chromosome['End'] - df_chromosome['Start']
        df_chromosome['position'] = df_chromosome['Start']
        df_chromosome['position'] = df_chromosome['position'].astype(int)
        df_chromosome['interval_length'] = df_chromosome['interval_length'].astype(int)

        df_new = expand_dataframe(df_chromosome)
        df_new['Chromosome'] = chromosome
        
        df_out = pd.concat([df_out, df_new])

    #rename the chromosome name to the S288C format by using the roman_to_S288C dictionary
    df_out['Chromosome'] = df_out['Chromosome'].str.replace('chr', '')

    from mayo.settings import roman_to_S288C
    df_out.replace({'Chromosome': roman_to_S288C}, inplace=True)

    df_out['Value'] = df_out['Value'].astype(float)

    return df_out

def check_for_duplicate_sample_names(samples_SNVs):
    #find duplicate names

    from collections import Counter
    names = [i.name for i in samples_SNVs]
    duplicates = [k for k,v in Counter(names).items() if v>1]
    print("Duplicates:", duplicates)

def filterCommonMutations_sectors_premeiotic(
        samples_SNVs, 
        specific_genotypes_2, 
        specific_genotypes_3, 
        specific_genotypes_4, 
        df_all_mutations, 
        genotype_dict,
        minimum_count_within_genotype=0,
        max_outside_genotype_duplicates=0):

    for specific_genotypes in (specific_genotypes_2 + specific_genotypes_3 + specific_genotypes_4):

        df_clean = filterCommonMutations_sectors(
            df_all_mutations=df_all_mutations,
            specific_genotypes=specific_genotypes,
            minimum_count_within_genotype=minimum_count_within_genotype,
            max_outside_genotype_duplicates=max_outside_genotype_duplicates)

        for sample in samples_SNVs:

            if sample.name not in genotype_dict.keys():
                continue

            initial_mutations = sample.df.shape[0]

            if genotype_dict[sample.name] in specific_genotypes:
                sample.df = sample.df[sample.df["mutation"].isin(df_clean["mutation"])]

                logging.info(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
                print(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
        
    # for specific_genotypes in specific_genotypes_2:

    #     df_clean = filterCommonMutations_sectors(
    #         df_all_mutations=df_all_mutations,
    #         specific_genotypes=specific_genotypes,
    #         minimum_count_within_genotype=minimum_count_within_genotype,
    #         max_outside_genotype_duplicates=max_outside_genotype_duplicates)
        
    #     for sample in samples_SNVs:

    #         if sample.name not in genotype_dict.keys():
    #             continue

    #         initial_mutations = sample.df.shape[0]

    #         if genotype_dict[sample.name] in specific_genotypes:
    #             sample.df = sample.df[sample.df["mutation"].isin(df_clean["mutation"])]

    #             logging.info(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
    #             print(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")

    # for specific_genotypes in specific_genotypes_3:

    #     df_clean = filterCommonMutations_sectors(
    #         df_all_mutations=df_all_mutations,
    #         specific_genotypes=specific_genotypes,
    #         minimum_count_within_genotype=minimum_count_within_genotype,
    #         max_outside_genotype_duplicates=max_outside_genotype_duplicates)
        
    #     for sample in samples_SNVs:

    #         if sample.name not in genotype_dict.keys():
    #             continue

    #         initial_mutations = sample.df.shape[0]

    #         if genotype_dict[sample.name] in specific_genotypes:
    #             sample.df = sample.df[sample.df["mutation"].isin(df_clean["mutation"])]

    #             logging.info(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
    #             print(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")

def showSpecificGenotypeClusters(
        samples_SNVs: list,
        specific_genotypes_2: list, 
        specific_genotypes_3: list,
        specific_genotypes_4: list,
        pval: float = 0.0001) -> None:
    
    specific_genotypes_all = specific_genotypes_2 + specific_genotypes_3 + specific_genotypes_4
    specific_genotypes_all = [item for sublist in (specific_genotypes_all) for item in sublist] #remove nests
    
    print("Tetrad/Dyad genotypes (Specific genotypes):")
    print(specific_genotypes_all)

    for sample in samples_SNVs:
        if sample.genotype in specific_genotypes_all:
            print(sample.name, sample.genotype)
            print(sample.cluster_dict_JT[pval])

def aggregateSectorsClusters(
        df_clusters: pd.DataFrame,
        specific_genotypes_2: list,
        specific_genotypes_3: list,
        specific_genotypes_4: list,
        genotype_dict: dict,
        remove_non_concordant: bool = True,
        cluster_outpath_individual: str = "outputs/cluster_table_with_missing_individuals_.tsv",
        cluster_outpath_totals: str = "outputs/cluster_table_with_missing_totals_.tsv"
        ) -> pd.DataFrame:

    genotypes_of_sector_samples = specific_genotypes_4 + specific_genotypes_2 + specific_genotypes_3
    genotypes_of_sector_samples = [item for sublist in genotypes_of_sector_samples for item in sublist]

    #get the sample_names that are in the genotypes_of_sector_samples
    reverse_genotype_dict = invert_dict(genotype_dict)
    sample_names = [reverse_genotype_dict[x] if x in reverse_genotype_dict.keys() else 'N/A' for x in genotypes_of_sector_samples]

    new_cluster_table = pd.DataFrame(columns=["Colony", "Chromosome", "Start_min", "End_max", "Length_max", "Mutations_max", "sample_count"])
    empty_cluster_table = pd.DataFrame(columns=["Colony", "Chromosome", "Start_min", "End_max", "Length_max", "Mutations_max", "sample_count"])

    for colony_samples in sample_names:

        if colony_samples == 'N/A':
            logging.warning(f"Some sectored samples are not detected in the genotype_dict and were returned as 'N/A'")
            logging.warning(f"Please check the genotypes_of_sector_samples list and the genotype_dict dictionary")
            logging.warning(f"genotypes_of_sector_samples: {genotypes_of_sector_samples}")
            continue

        colony_df = df_clusters[df_clusters["sample"].isin(colony_samples)].copy()
        colony_samples_string = "-".join(colony_samples)
        
        #we need to always add empty row in case we end up removing all clusters for a sample
        #not removing only when there is nothing there to begin with
        new_empty_row = {
            "Colony": colony_samples_string, 
            "Chromosome": None,
            "Start_min": None,
            "End_max": None,
            "Length_max": 0,
            "Mutations_max": 0,
            "sample_count": 0}
            
        empty_cluster_table = pd.concat([empty_cluster_table, pd.DataFrame([new_empty_row])], ignore_index=True)

        if colony_df.shape[0] == 0:
            logging.info(f"No clusters found in {colony_samples_string}")
            continue

        logging.info(f"Analyzing clusters in {colony_samples_string}")

        #group by Chromosome, Start, End, Length, Mutations and get count od sample and max of mutations
        colony_df = colony_df.groupby(["Chromosome", "Start", "End", "Length"]).agg({"sample": "count", "Mutations": "max"}).reset_index()
        colony_df.rename(columns={"Mutations": "max_mutations", "sample": "sample_count"}, inplace=True)
        logging.info(colony_df)

        #go through chromosome by chromosome in colony_df
        for chr in colony_df["Chromosome"].unique():
            colony_df_chr = colony_df[colony_df["Chromosome"] == chr].copy()

            #sort Start (ascending) and End (descending)
            colony_df_chr.sort_values(by=["Start", "End"], ascending=[True, False], inplace=True)

            #get the End position of the first cluster
            curr_cluster_start = colony_df_chr.iloc[0]["Start"]
            curr_cluster_end = colony_df_chr.iloc[0]["End"]
            max_mutations = colony_df_chr.iloc[0]["max_mutations"]
            total_samples = colony_df_chr.iloc[0]["sample_count"]
            logging.debug(f"chr:{chr} cluster: {curr_cluster_start} - {curr_cluster_end}")

            if colony_df_chr.shape[0] == 1:
                new_row = {
                    "Colony": colony_samples_string, 
                    "Chromosome": chr, 
                    "Start_min": curr_cluster_start, 
                    "End_max": curr_cluster_end, 
                    "Length_max": curr_cluster_end - curr_cluster_start, 
                    "Mutations_max": max_mutations,
                    "sample_count": total_samples}
                new_cluster_table = pd.concat([new_cluster_table, pd.DataFrame([new_row])], ignore_index=True)
                continue
            
            #go through the rest of the clusters
            for i in range(1, colony_df_chr.shape[0]):

                #if the next cluster is disjointed from the current cluster, save the current cluster
                if colony_df_chr.iloc[i]["Start"] > curr_cluster_end:
                    logging.debug(f"next chr:{chr} cluster is disjointed from {curr_cluster_start} - {curr_cluster_end}")
                    
                    #save old cluster by concatenating to new_cluster_table and concatenate to new_cluster_table
                    new_row = {
                        "Colony": colony_samples_string, 
                        "Chromosome": chr, 
                        "Start_min": curr_cluster_start, 
                        "End_max": curr_cluster_end, 
                        "Length_max": curr_cluster_end - curr_cluster_start, 
                        "Mutations_max": colony_df_chr.iloc[i-1]["max_mutations"],
                        "sample_count": total_samples}
                    new_cluster_table = pd.concat([new_cluster_table, pd.DataFrame([new_row])], ignore_index=True)

                    #start a new cluster
                    curr_cluster_start = colony_df_chr.iloc[i]["Start"]
                    curr_cluster_end = colony_df_chr.iloc[i]["End"]
                    max_mutations = colony_df_chr.iloc[i]["max_mutations"]
                    logging.debug(f"chr:{chr} cluster: {curr_cluster_start} - {curr_cluster_end}")
                
                #else next cluster has to to overlap with the current cluster
                else:

                    #add the sample count to the total_samples
                    current_samples = colony_df_chr.iloc[i]["sample_count"]
                    total_samples += current_samples

                    #extend the current cluster if end is higher than the current end
                    if colony_df_chr.iloc[i]["End"] > curr_cluster_end:
                        curr_cluster_end = colony_df_chr.iloc[i]["End"]
                        logging.debug(f"current chr:{chr} cluster is extended to: {curr_cluster_start} - {curr_cluster_end}")

                    #if the max_mutations is higher than the current max_mutations, update it
                    if colony_df_chr.iloc[i]["max_mutations"] > max_mutations:
                        logging.debug(f"max_mutations updated from {max_mutations} to {colony_df_chr.iloc[i]['max_mutations']}")
                        max_mutations = colony_df_chr.iloc[i]["max_mutations"]

            #save the last cluster in the chromosome
            new_row = {
                "Colony": colony_samples_string, 
                "Chromosome": chr, 
                "Start_min": curr_cluster_start, 
                "End_max": curr_cluster_end, 
                "Length_max": curr_cluster_end - curr_cluster_start, 
                "Mutations_max": max_mutations,
                "sample_count": total_samples}
            new_cluster_table = pd.concat([new_cluster_table, pd.DataFrame([new_row])], ignore_index=True)

    new_cluster_table["sectors_strategy"] = new_cluster_table["Colony"].apply(lambda x: len(x.split("-")))

    #remove samples with fewer than 3 sample count 
    if remove_non_concordant:
        new_cluster_table = new_cluster_table[
            (new_cluster_table["sample_count"] >= 2) & (new_cluster_table["sectors_strategy"] == 2) | #for samples with 2 sectors
            (new_cluster_table["sample_count"] >= 3) & (new_cluster_table["sectors_strategy"] != 2)   #for samples with 3 or more sectors
            ]

    #add empty clusters to new_cluster_table 
    cluster_table_with_missing = pd.concat([new_cluster_table, empty_cluster_table], ignore_index=True)
    
    logging.info(cluster_table_with_missing)
    print(cluster_table_with_missing)
    
    #save cluster_table_with_missing to file
    cluster_table_with_missing.to_csv(cluster_outpath_individual, sep="\t", index=False)

    #groupby Colony and sum the Length_max and Mutations_max
    cluster_table_with_missing_totals = cluster_table_with_missing.groupby("Colony").agg({"Length_max": "sum", "Mutations_max": "sum", "Start_min": "count"}).reset_index()

    #rename Length_max to Total_length and Mutations_max to Total_mutations
    cluster_table_with_missing_totals.rename(columns={
        "Length_max": "Total_length", 
        "Mutations_max": "Total_clust_mutations",
        "Start_min": "Total_clusters"
        }, inplace=True)

    #if total_length is 0, set total_clusters and total_clust_mutations to 0
    cluster_table_with_missing_totals["Total_clusters"] = cluster_table_with_missing_totals.apply(lambda x: 0 if x["Total_length"] == 0 else x["Total_clusters"], axis=1)
    cluster_table_with_missing_totals["Total_clust_mutations"] = cluster_table_with_missing_totals.apply(lambda x: 0 if x["Total_length"] == 0 else x["Total_clust_mutations"], axis=1)
    
    logging.info(cluster_table_with_missing_totals)
    print(cluster_table_with_missing_totals)
    
    #save cluster_table_with_missing_totals to file
    cluster_table_with_missing_totals.to_csv(cluster_outpath_totals, sep="\t", index=False)

def filterCommonMutations_sectors_postmeiotic(
        samples_SNVs, 
        specific_genotypes_2, 
        specific_genotypes_3, 
        specific_genotypes_4, 
        df_all_mutations, 
        genotype_dict,
        minimum_count_within_genotype_2=2,
        minimum_count_within_genotype_3=3,
        minimum_count_within_genotype_4=3,
        max_outside_genotype_duplicates=0):
        
    for specific_genotypes in specific_genotypes_2:

        df_clean = filterCommonMutations_sectors(
            df_all_mutations=df_all_mutations,
            specific_genotypes=specific_genotypes,
            minimum_count_within_genotype=minimum_count_within_genotype_2,
            max_outside_genotype_duplicates=max_outside_genotype_duplicates)
        
        for sample in samples_SNVs:

            if sample.name not in genotype_dict.keys():
                continue

            initial_mutations = sample.df.shape[0]

            if genotype_dict[sample.name] in specific_genotypes:
                sample.df = sample.df[sample.df["mutation"].isin(df_clean["mutation"])]

                logging.info(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
                print(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")

    for specific_genotypes in specific_genotypes_3:

        df_clean = filterCommonMutations_sectors(
            df_all_mutations=df_all_mutations,
            specific_genotypes=specific_genotypes,
            minimum_count_within_genotype=minimum_count_within_genotype_3,
            max_outside_genotype_duplicates=max_outside_genotype_duplicates)
        
        for sample in samples_SNVs:

            if sample.name not in genotype_dict.keys():
                continue

            initial_mutations = sample.df.shape[0]

            if genotype_dict[sample.name] in specific_genotypes:
                sample.df = sample.df[sample.df["mutation"].isin(df_clean["mutation"])]

                logging.info(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
                print(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")

    for specific_genotypes in specific_genotypes_4:

        df_clean = filterCommonMutations_sectors(
            df_all_mutations=df_all_mutations,
            specific_genotypes=specific_genotypes,
            minimum_count_within_genotype=minimum_count_within_genotype_4,
            max_outside_genotype_duplicates=max_outside_genotype_duplicates)

        for sample in samples_SNVs:

            if sample.name not in genotype_dict.keys():
                continue

            initial_mutations = sample.df.shape[0]

            if genotype_dict[sample.name] in specific_genotypes:
                sample.df = sample.df[sample.df["mutation"].isin(df_clean["mutation"])]

                logging.info(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")
                print(f"{sample.name}: Removed {initial_mutations - sample.df.shape[0]:3} non-meiotic mutations ({initial_mutations:3} -> {sample.df.shape[0]:3})")

def get_rec_ev_relative_A3A_SNVs_mod(
        samples_SNVs, 
        genotype_dict, 
        genotype_list, 
        mutation_type: Literal['clustered', 'scattered', 'all'] = "all", 
        window_size = 8000
        ) -> pd.DataFrame:
    """
    This function calculates mutations around recombination events for each genotype in the genotype list, plots the results, and returns a DataFrame with the results for all genotypes.

    Parameters:
    samples_SNVs (DataFrame): A DataFrame containing information about the SNVs in the samples.
    genotype_dict (dict): A dictionary mapping genotypes to their corresponding samples.
    genotype_list (list): A list of genotypes to consider.
    clustered_only (bool): If True, only consider clustered SNVs.

    Returns:
    DataFrame: A DataFrame containing information about the mutations around recombination events for all genotypes in the genotype list.
    """
    meta_master_relative_GC_CO_df = pd.DataFrame()

    for focus_genotype in genotype_list:

        print(f"Calculating mutations around recombination events for {focus_genotype} samples...")

        master_relative_GC_CO_df = get_mutations_around_rec_ev(
            samples_SNVs,
            genotype_dict,
            focus_genotype=focus_genotype, #pol32∆ None ung1∆+ung1∆NAT+exo1-nd+pol32∆+exo1-ndpol32∆
            mutation_type=mutation_type,
            homozygous_only=True,
            event_type="both",
            window_size=window_size)
        
        #if Reference == G->A, then multiply the Region by -1 and change the Reference to C->T
        # master_relative_GC_CO_df.loc[master_relative_GC_CO_df["Reference"] == "G->A", "Region"] = master_relative_GC_CO_df["Region"] * -1
        # master_relative_GC_CO_df["Reference"] = master_relative_GC_CO_df["Reference"].replace(["G->A"], "C->T")
        master_relative_GC_CO_df["genotype"] = focus_genotype
        meta_master_relative_GC_CO_df = pd.concat([meta_master_relative_GC_CO_df, master_relative_GC_CO_df], ignore_index=True)

        #calculate the number of C->T mutations in the negative region
        c_t_neg = len(master_relative_GC_CO_df[(master_relative_GC_CO_df['Region'] < 0) & (master_relative_GC_CO_df['Reference'] == 'C->T')])
        g_a_neg = len(master_relative_GC_CO_df[(master_relative_GC_CO_df['Region'] < 0) & (master_relative_GC_CO_df['Reference'] == 'G->A')])
        c_t_pos = len(master_relative_GC_CO_df[(master_relative_GC_CO_df['Region'] > 0) & (master_relative_GC_CO_df['Reference'] == 'C->T')])
        g_a_pos = len(master_relative_GC_CO_df[(master_relative_GC_CO_df['Region'] > 0) & (master_relative_GC_CO_df['Reference'] == 'G->A')])

        print(f"C->T mutations: negative_region={c_t_neg}, positive_region={c_t_pos}, negative_c_t/negative_total={c_t_neg / (c_t_neg + g_a_neg)}")
        print(f"G->A mutations: negative_region={g_a_neg}, positive_region={g_a_pos}, positive_g_a/positive_total={g_a_pos / (c_t_pos + g_a_pos)}")

        print(f"Ratio of C->T mutations in the negative region: {c_t_neg / (c_t_neg + g_a_neg)}")
        print(f"Ratio of G->A mutations in the positive region: {g_a_pos / (c_t_pos + g_a_pos)}")

        resection_pattern_mutations = c_t_neg + g_a_pos
        bir_pattern_mutations = g_a_neg + c_t_pos

        print(f"Resection pattern mutations: {resection_pattern_mutations}")
        print(f"BIR pattern mutations: {bir_pattern_mutations}")


    return meta_master_relative_GC_CO_df

def conduct_mann_whitney_U_test(
        meta_master_relative_GC_CO_df: pd.DataFrame,
        region_flank: int | float = 8,
        bins: int = 64
        ) -> None:
    """
    Conducts a Mann-Whitney U test for the distribution of SNVs around recombination events.

    Parameters:
    meta_master_relative_GC_CO_df (pd.DataFrame): DataFrame containing information about the SNVs and their relative positions to recombination events for a meta-analysis.
    region_flank (int | float): Size of the region around the recombination event to consider. Default is 8.
    bins (int): Number of bins to use for the distribution. Default is 64.

    Returns:
    None
    """
    from scipy.stats import mannwhitneyu

    region = (-region_flank, region_flank)

    df_trimmed = bin_rec_relative_mutations(meta_master_relative_GC_CO_df, region, bins=bins)
    df_norm = normalize_rec_relative_bins_to_wt(df_trimmed, normalize_to="ung1∆")

    print("Dataframe with normalized counts to ung1∆ for plotting:")
    print(df_norm)

    #if there are nan values, drop them
    df_norm = df_norm.dropna()

    #conduct test for 2 intervals separately, left and right of the recombination event (less than 0 and greater than 0)
    df_norm["Region_bin"] = df_norm["Region_bin"].astype(float)
    df_left_norm = df_norm[df_norm["Region_bin"] < 0]
    df_right_norm = df_norm[df_norm["Region_bin"] > 0]

    for focus_genotype in df_norm.columns[1:].tolist():

        #conduct Mann-Whitney U test
        print(f"Conducting Mann-Whitney U test for {focus_genotype}...")
        
        print(f"Left (symmetric component) of the recombination event (less than 0) for {focus_genotype}:")
        u_statistic, p_value = mannwhitneyu(df_left_norm["ung1∆"], df_left_norm[focus_genotype])
        print(f"U statistic: {u_statistic}, p-value: {p_value} for {focus_genotype}")

        print(f"Right (asymmetric component) of the recombination event (greater than 0) for {focus_genotype}:")
        u_statistic, p_value = mannwhitneyu(df_right_norm["ung1∆"], df_right_norm[focus_genotype])
        print(f"U statistic: {u_statistic}, p-value: {p_value} for {focus_genotype}")
         
def plot_GC_CO_distribution_by_genotype_box(
        samples_SNVs, 
        query_genotypes: str | None = None,
        save_name: str | None = "distribution_haploid",
        show: bool = True
        ) -> None:
    """
    This function plots the distribution of GCs (Gene Conversions) and COs (Crossover events) for the given genotypes. 
    It also calculates and prints the average, median, and standard deviation of GCs and COs.

    Parameters:
    samples_SNVs (DataFrame): A DataFrame containing information about the SNVs in the samples.
    query_genotypes (str | None): A string representing the genotypes to consider. If None, all genotypes are considered.
    save_name (str | None): The name to use when saving the plot. If None, the plot is not saved.
    show (bool): If True, the plot is displayed. If False, the plot is not displayed.

    Returns:
    None
    """    
    import statistics

    GC_numbers, CO_numbers = get_GC_CO_numbers(samples_SNVs, query_genotypes=query_genotypes, min_percent=80)

    print(f"There are a total of {sum(GC_numbers)} GCs and {sum(CO_numbers)} COs in {len(GC_numbers)}:{len(CO_numbers)} samples.")

    #combine the GC and CO numbers into one dataframe for plotting
    GC_CO_df = pd.DataFrame({"Event": ["GC"] * len(GC_numbers) + ["CO"] * len(CO_numbers), "Number": GC_numbers + CO_numbers})

    sns.set_context("poster")
    plt.figure(figsize=(6, 6))
    sns.boxplot(
        data=GC_CO_df,
        x="Event",
        y="Number",
        linewidth=2,
        color="silver",
        boxprops=dict(edgecolor="black"))   
    plt.title("Distribution of GCs and COs", fontsize=24, fontweight="bold")
    plt.ylim(0, max(GC_numbers + CO_numbers) + 5)
    plt.xticks(fontsize=18, fontweight="bold")
    plt.yticks(fontsize=18, fontweight="bold")
    plt.ylabel("Number of Events", fontsize=20, fontweight="bold")
    plt.xlabel("Event Type", fontsize=20, fontweight="bold")

    if save_name is not None:
        plt.savefig(f"figures/GC_and_CO_{save_name}_{query_genotypes}.png", dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()
    
    average_GC = sum(GC_numbers) / len(GC_numbers)
    median_GC = statistics.median(GC_numbers)
    stdev_GC = statistics.stdev(GC_numbers)

    average_CO = sum(CO_numbers) / len(CO_numbers)
    median_CO = statistics.median(CO_numbers)
    stdev_CO = statistics.stdev(CO_numbers)

    print(f"Query Genotypes: {query_genotypes}, GCs: Average={average_GC}, Median={median_GC}, Stdev={stdev_GC}")
    print(f"Query Genotypes: {query_genotypes}, COs: Average={average_CO}, Median={median_CO}, Stdev={stdev_CO}")

def prep_yeast_mine_data(path: str = ""):
    """Prepares yeast mine data for use in this script"""

    df_UTR = pd.read_csv(path)
    df_UTR = df_UTR[df_UTR['Gene > Transcripts > In _ Ypd'] == True]
    #df_UTR = df_UTR[df_UTR['Gene > Transcripts > In _ Gal'] == True]


    df_UTR["utr_type"] = df_UTR["Gene > Transcripts > UTRs > UTR > Primary DBID"].str.split("-").str[1]

    df_UTR["feture_direction"] = df_UTR["Gene > Chromosome Location > Strand"].map({"1": "+", "-1": "-"})

    #rename Gene > Transcripts > UTRs > UTR > Start Location to start
    #rename Gene > Transcripts > UTRs > UTR > End Location to end
    df_UTR = df_UTR.rename(columns={
        "Gene > Transcripts > UTRs > UTR > Start Location": "start", 
        "Gene > Transcripts > UTRs > UTR > End Location": "end",
        "Gene > Systematic Name": "gene_name"})

    df_UTR_3 = df_UTR[df_UTR['utr_type'] == "3prime"]
    df_UTR_5 = df_UTR[df_UTR['utr_type'] == "5prime"]

    return df_UTR, df_UTR_3, df_UTR_5

def pull_mutations_relative_to_feature(
        df_mutations: pd.DataFrame,
        chromosome: str,
        position: int,
        feature_direction: Literal["+", "-"],
        window_size = 10000):
    """
    Filters a DataFrame of mutations based on their relative position to a specific feature.

    Args:
        df_mutations (pd.DataFrame): The DataFrame containing mutation data.
        chromosome (str): The chromosome where the feature is located.
        position (int): The position of the feature on the chromosome.
        feature_direction (Literal["+", "-"]): The direction of the feature. "+" for forward, "-" for reverse.
        window_size (int, optional): The size of the window around the feature to consider. Defaults to 10000.

    Returns:
        pd.DataFrame: A DataFrame containing only the mutations within the specified window around the feature.
        If the feature direction is "-", the 'Region', 'Reference', and 'Allele' columns are inverted.
    """
    import BioAid as ba
    
    df_chr_pos = pull_mutations_relative_to_position(df_mutations, chromosome, position, window_size)

    #check that the mutations are in fact within -10000 and +10000 of the feature start
    # df_chr_pos = df_chr_pos[df_chr_pos['Region'] >= 0]
    # df_chr_pos = df_chr_pos[df_chr_pos['Region'] <= window_size]

    if feature_direction == '+':
        return df_chr_pos
    
    elif feature_direction == '-':
        df_chr_pos['Region'] = df_chr_pos['Region'] * -1

        #also need to flip the mutation type
        df_chr_pos['Reference'] = df_chr_pos['Reference'].apply(lambda x: ba.compl(x))
        df_chr_pos['Allele'] = df_chr_pos['Allele'].apply(lambda x: ba.compl(x))

        return df_chr_pos

def create_mutations_around_features_start_df(
        df_mutations, 
        features_path: str,
        query_position: Literal["Start", "End"] = "Start",
        window_size: int = 10000,
        just_tRNA_genes: bool = False,
        just_non_essential_genes: list | None = None,
        tss_and_tts_from_utrs: bool = True,
        UTR_data_path: str = "data/yeastmine_results_2024-05-24T10-03-55.csv"

        ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    This function creates mutations around the start of features in a DataFrame.

    Parameters:
    df_mutations (pd.DataFrame): DataFrame containing mutation data.
    features_path (str): Path to the features fasta file.
    query_position (Literal["Start", "End"], optional): Position to query for mutations. Defaults to "Start".
    window_size (int, optional): Size of the window to look for mutations around the feature start. Defaults to 10000.
    just_tRNA_genes (bool, optional): If True, only consider tRNA genes. Defaults to False.
    just_non_essential_genes (list, optional): List of non-essential genes to consider. Defaults to None.

    Returns:
    tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames. 
    The first DataFrame contains the mutations around the feature start. The second DataFrame contains the features data.
    """
    import BioAid as ba

    df_mutations = df_mutations.sort_values(by=["Chromosome", "Region"])
    df_features = pd.DataFrame()
    df_feature_mutations_master = pd.DataFrame()
    repetitive_DNA_list = ba.extractSeqFromFastaToList(features_path) #open the features fasta file with BioAid

    df_UTR, df_UTR_3, df_UTR_5 = prep_yeast_mine_data(path=UTR_data_path)

    analyzed_features = 0

    for seq_name, seq in repetitive_DNA_list:

        #extract name from seq_name
        feature_name = seq_name.split(" ")[0].strip(">").strip()
        chromosome_name = seq_name.split(",")[1].split(" ")[2].strip()

        if just_non_essential_genes is not None:
            if feature_name not in just_non_essential_genes:
                continue

        if just_tRNA_genes:
            if "tRNA_gene" not in seq_name: #tRNA_gene ARS
                continue

        # if "ARS" not in seq_name:
        #     continue
        # print(f"Evaluating {seq_name}")

        from mayo.settings import roman_to_S288C

        if chromosome_name not in roman_to_S288C.keys():
            logging.warning(f"Chromosome {chromosome_name} of {feature_name} is not in roman_to_S288C.keys(). This behavior is normal for Mitochondrial Sequences. Skipping...")
            continue

        #print(f"Evaluating {seq_name}")
        chromosome_name = roman_to_S288C[chromosome_name]
        seq_coords = seq_name.split(",")[1].split(" ")[-1].strip()
        seq_start = int(seq_coords.split("-")[0])
        seq_end = int(seq_coords.split("-")[1])

        if seq_start > seq_end:
            feature_direction = "-"

        else:
            feature_direction = "+"

        if tss_and_tts_from_utrs:

            if feature_name in df_UTR_5['gene_name'].values:
                #extract the 5' UTR start and end positions
                utr_seq_start = df_UTR_5[df_UTR_5['gene_name'] == feature_name]['start'].values[0]
                utr_seq_end = df_UTR_5[df_UTR_5['gene_name'] == feature_name]['end'].values[0]

                if feature_direction == "+":
                    seq_start = utr_seq_start
                if feature_direction == "-":
                    seq_start = utr_seq_end

                #print(f"Feature {feature_name} is in the 5' UTR of {chromosome_name} at {seq_start}-{seq_end}, *{feature_direction}*")

            if feature_name in df_UTR_3['gene_name'].values:

                #extract the 3' UTR start and end positions
                utr_seq_start = df_UTR_3[df_UTR_3['gene_name'] == feature_name]['start'].values[0]
                utr_seq_end = df_UTR_3[df_UTR_3['gene_name'] == feature_name]['end'].values[0]

                if feature_direction == "+":
                    seq_end = utr_seq_end
                if feature_direction == "-":
                    seq_end = utr_seq_start

                #print(f"Feature {feature_name} is in the 3' UTR of {chromosome_name} at {seq_start}-{seq_end}, *{feature_direction}*")

            #else:
                #print(f"Feature {feature_name} is not in the UTRs")


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

        #sort the genomic coordinates, repetitive DNA file has some coordinates in reverse order (transcription direction)
        seq_genomic_coord = [seq_start, seq_end]
        seq_genomic_coord.sort()

        if query_position == "Start":
            query_position_var = seq_start
        elif query_position == "End":
            query_position_var = seq_end
        else:
            raise ValueError(f"query_position {query_position} not recognized")

        df_feature_mutations = pull_mutations_relative_to_feature(
            df_mutations=df_mutations,
            chromosome=chromosome_name, 
            position=query_position_var, #or seq_start
            feature_direction=feature_direction,
            window_size=window_size)
        
        df_feature_mutations['Feature'] = feature_name

        df_features_append = pd.DataFrame([[
            feature_name, 
            chromosome_name, 
            seq_start, 
            seq_end, 
            feature_direction]], 
            columns=[
                "Feature", 
                "Chromosome", 
                "Start", 
                "End", 
                "Direction"])
        
        df_features = pd.concat([df_features, df_features_append], ignore_index=True)

        df_feature_mutations_master = pd.concat([df_feature_mutations_master, df_feature_mutations], ignore_index=True)

        analyzed_features += 1

    print(f"Evaluated {analyzed_features} features in {features_path}")

    return df_feature_mutations_master, df_features

def plot_mutations_around_feature(
        df_mutations: pd.DataFrame, 
        title: str = "Mutations around Feature Position", 
        window_size: int = 1000,
        show: bool = False
        ):
    """
    Plots a histogram of mutations around a feature start.

    Args:
        df_mutations (pd.DataFrame): The DataFrame containing mutation data.
        title (str, optional): The title of the plot. Defaults to "Mutations around Feature Position".
        window_size (int, optional): The size of the window around the feature to consider. Defaults to 2000.

    Returns:
        None: This function doesn't return anything; it shows a plot.
    """
    from matplotlib import ticker

    df_mutations = df_mutations[(df_mutations['Region'] <= window_size ) & (df_mutations['Region'] >= -window_size)]
    #print(f"Number of mutations around feature: {df_mutations.shape[0]}, sites={(window_size*2) + 1}")

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.set_context("poster")

    #set pallete for the residues
    pallette = {'C': 'tab:blue', 'G': 'tab:orange', 'T': 'tab:red', 'A': 'tab:green'}
    hue_order = ['C', 'G', 'T', 'A']

    sns.histplot(data=df_mutations, x="Region", hue="Reference", ax=ax, bins=200, palette=pallette, hue_order=hue_order, element="step", alpha=0.5, linewidth=0)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1) #draw a vertical line at 0
    ax.set_title(title, fontsize=24, weight='bold')
    ax.set_xlabel("Distance from Feature (bp)", fontsize=20, weight='bold')
    ax.set_ylabel("Number of Mutations", fontsize=20, weight='bold')
    ax.set_xlim(-window_size, window_size)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
    plt.tight_layout()

    plt.savefig(f"figures/transcription/{title}_{window_size}_flank.png", dpi=300)
    if show:
        plt.show()
    else:
        plt.close()

def countA3A_motifs_in_seq(seq: str) -> tuple:
    """
    Counts the number of A3A motifs in a given sequence.

    A3A motifs are defined as a 'C' preceded by a 'T' and not followed by a 'C' in the forward direction,
    or as a 'G' followed by an 'A' and not preceded by a 'G' in the reverse direction.

    Args:
        seq (str): The sequence to search for A3A motifs.

    Returns:
        tuple: A tuple containing two numpy arrays. The first array represents the positions of A3A motifs in the forward direction,
               and the second array represents the positions of A3A motifs in the reverse direction. The arrays have the same length as the input sequence,
               and positions with an A3A motif are marked with a 1, while all other positions are marked with a 0.
    """
    import numpy as np
    seq_len = len(seq)

    #create an array of 0s with the same length as the sequence
    seq_array_fw = np.zeros(seq_len, dtype=int)
    seq_array_rv = np.zeros(seq_len, dtype=int)

    for index, letter in enumerate(seq):
        if (index == 0) or (index == seq_len-1):
            continue

        if letter == "C":
            if seq[index-1] == "T":
                if seq[index+1] != "C":
                    seq_array_fw[index] = 1

        elif letter == "G":
            if seq[index+1] == "A":
                if seq[index-1] != "G":
                    seq_array_rv[index] = 1

    return seq_array_fw, seq_array_rv

def plot_A3A_motifs_relative_to_position(
        df: pd.DataFrame, 
        suffix:str = "",
        show: bool = False
        ) -> None:
    """
    Plots the number of A3A motifs relative to the feature start.

    Args:
        df (pd.DataFrame): The DataFrame containing the A3A motif data. It should have columns "Distance", "Reverse", and "Forward".

    Returns:
        None: This function doesn't return anything; it shows a plot.
    """
    from matplotlib import ticker
    sns.set_context("poster")

    fig, ax = plt.subplots(figsize=(12, 8))

    #for some reason matplotlib is giving warnings here
    sns.barplot(data=df, x="Distance", y="Reverse", color="tab:orange", ax=ax, alpha = 0.5, width=1, linewidth=0)
    sns.barplot(data=df, x="Distance", y="Forward", color="tab:blue", ax=ax, alpha = 0.5, width=1, linewidth=0)
    ax.set_title("A3A Motifs around Feature Position", fontsize=24, weight='bold')
    ax.set_xlabel("Distance from Feature Position (bp)", fontsize=20, weight='bold')
    ax.set_ylabel("Number of A3A Motifs", fontsize=20, weight='bold')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
    plt.tight_layout()
    plt.savefig(f"figures/transcription/A3A_motifs_{suffix}.png", dpi=300)

    if show:
        plt.show()
    else:
        plt.close()

def get_nucleotide_frequencies_around_feature(
        df_features: pd.DataFrame,
        context_flank: int = 100,
        query_position: Literal["Start", "End"] = "Start",
        fasta_file_path: str = "S288C_reference_sequencewoUra3wUra329wNAT.fa"
        ) -> pd.DataFrame:
    
    import BioAid as ba
    import numpy as np

    df = df_features.copy()

    list_of_positions = df[['Chromosome', query_position]].values.tolist()

    result = ba.pullGenomicContext(list_of_positions, fasta_file_path, context_flank)

    df[['ctx_L', 'ctx_C', 'ctx_R', 'ctx_chr']] = pd.DataFrame(result, index=df.index)
        
    df['seq'] = df['ctx_L'] + df['ctx_C'] + df['ctx_R']

    #make sure that df["seq"]'s length is equal to context_flank * 2 + 1, if not, check which features are causing the problem
    df['seq_len'] = df['seq'].apply(lambda x: len(x))

    #remove the features that are causing the problem
    print("removing features that are causing the problem:")
    print(df[df['seq_len'] != context_flank * 2 + 1][['Feature', 'seq_len']].values.tolist())
    df = df[df['seq_len'] == context_flank * 2 + 1]
    assert np.array_equal(df['seq'].apply(lambda x: len(x)).unique(), [context_flank * 2 + 1]), "Length of sequence is not equal to context_flank * 2 + 1"

    #if the feature is on the reverse strand, we need to reverse complement the sequence
    df['seq'] = df.apply(lambda x: ba.rev_compl(x['seq']) if x['Direction'] == '-' else x['seq'], axis=1)

    #create a new df to store nucleotide frequencies
    df_nucleotide_freq = pd.DataFrame()
    
    #create len(seq) columns that will store the value of the nucleotide at that position
    # Collect all new columns in a list
    new_columns = []

    for i in range(context_flank * 2 + 1):
        #df_nucleotide_freq[f"pos_{i}"] = df['seq'].apply(lambda x: x[i]) #old way, triggered a performance warning
        new_col = df['seq'].apply(lambda x: x[i])
        new_col.name = f"pos_{i}"
        new_columns.append(new_col)

    # Concatenate all new columns to the DataFrame at once
    df_nucleotide_freq = pd.concat([df_nucleotide_freq] + new_columns, axis=1)
    

    #count the number of each nucleotide at each position(column) calculate the frequency of each nucleotide and store it in a list
    #iterate over the columns of the new df
    nucleotide_freq_list = []
    for col in df_nucleotide_freq.columns:
        #count the number of each nucleotide at each position
        value_counts = df_nucleotide_freq[col].value_counts()
        #calculate the frequency of each nucleotide
        freq = value_counts / value_counts.sum()
        #convert value_counts to a dictionary
        freq = freq.to_dict()    
        #store the frequency of each nucleotide in a list
        nucleotide_freq_list.append(freq)


    #convert the list to a DataFrame, with nucleotides being rows and positions being columns
    df_nucleotide_freq = pd.DataFrame(nucleotide_freq_list)

    return df_nucleotide_freq

def countA3A_motifs_relative_to_position(
        df_features: pd.DataFrame,
        context_flank: int = 100,
        query_position: Literal["Start", "End"] = "Start",
        fasta_file_path: str = "S288C_reference_sequencewoUra3wUra329wNAT.fa"
        ) -> pd.DataFrame:
    """
    This function counts A3A motifs relative to a position in a given DataFrame. 
    It first copies the DataFrame and extracts the list of positions. 
    It then pulls the genomic context for these positions and adds this context to the DataFrame. 
    It also checks the length of the sequences and removes any features causing problems. 
    If the feature is on the reverse strand, it reverse complements the sequence. 
    It then counts the number of A3A motifs in the sequence and sums up these counts. 
    Finally, it converts the array of counts to a DataFrame and returns it.

    Parameters:
    df_features (pd.DataFrame): A DataFrame containing the features data. It should have the columns: 
        "Chromosome", "Start", "End", "Feature", "Direction", and "seq".
    context_flank (int, optional): The number of base pairs to include on either side of the feature. Defaults to 100.
    query_position (Literal["Start", "End"], optional): The position to query. Defaults to "Start".
    fasta_file_path (str, optional): The path to the fasta file. Defaults to a specific file path.

    Returns:
    df (pd.DataFrame): A DataFrame containing the counts of A3A motifs relative to the position. 
    It has the columns "Distance", "Forward", and "Reverse".
    """
    import BioAid as ba
    import numpy as np

    df = df_features.copy()

    list_of_positions = df[['Chromosome', query_position]].values.tolist()

    result = ba.pullGenomicContext(list_of_positions, fasta_file_path, context_flank)

    df[['ctx_L', 'ctx_C', 'ctx_R', 'ctx_chr']] = pd.DataFrame(result, index=df.index)
        
    df['seq'] = df['ctx_L'] + df['ctx_C'] + df['ctx_R']

    #make sure that df["seq"]'s length is equal to context_flank * 2 + 1, if not, check which features are causing the problem
    df['seq_len'] = df['seq'].apply(lambda x: len(x))

    #remove the features that are causing the problem
    print("removing unbound features (extending past reference genome):")
    print(df[df['seq_len'] != context_flank * 2 + 1][['Feature', 'seq_len']].values.tolist())
    df = df[df['seq_len'] == context_flank * 2 + 1]
    assert np.array_equal(df['seq'].apply(lambda x: len(x)).unique(), [context_flank * 2 + 1]), "Length of sequence is not equal to context_flank * 2 + 1"
    #print number of remaining features
    print(f"Number of features remaining: {len(df)}")


    #if the feature is on the reverse strand, we need to reverse complement the sequence
    df['seq'] = df.apply(lambda x: ba.rev_compl(x['seq']) if x['Direction'] == '-' else x['seq'], axis=1)

    #find the number of A3A motifs in the sequence
    df[['A3A_motifs_fw', 'A3A_motifs_rv']] = pd.DataFrame(df['seq'].apply(lambda x: countA3A_motifs_in_seq(x)).tolist(), index=df.index)

    array_fw = np.array(df['A3A_motifs_fw'].tolist())
    array_rv = np.array(df['A3A_motifs_rv'].tolist())

    #sum up both arrays vertically
    array_fw = np.sum(array_fw, axis=0)
    array_rv = np.sum(array_rv, axis=0)

    #convert the array to a pandas dataframe
    df = pd.DataFrame({"Distance": np.arange(-context_flank, context_flank+1), "Forward": array_fw, "Reverse": array_rv})

    return df

def normalize_mutations_around_feature(
        df_feature_mutations_master_coding, 
        df_features_mod, 
        feature_range=(-1000, 1000),
        sample_size_norm: int = 1,
        feature_count_norm: int = 1) -> pd.DataFrame:
    """
    This function normalizes the mutations around a feature. It first groups the data by "Region", 
    "Reference", and "Allele" and calculates the count of each group. 
    Then it joins this data with another DataFrame based on the "Region" and "Distance" columns. 
    It also creates a new column 'spectra' based on the 'Reference' and 'Allele' columns. 
    The function then normalizes the count based on the 'Forward' and 'Reverse' columns and the 'spectra' column. 
    Finally, it smooths the data for plotting purposes.

    Parameters:
    df_feature_mutations_master_coding (pd.DataFrame): A DataFrame containing the mutations data. It should have the 
    columns "Region", "Reference", and "Allele".
    df_features_mod (pd.DataFrame): A DataFrame containing the features data. It should have a column "Distance".
    feature_range (tuple): A tuple specifying the range of features to consider. Default is (-1000, 1000).
    sample_size_norm (int): The sample size to normalize the counts. Default is 1.
    feature_count_norm (int): The feature count to normalize the counts. Default is 1.

    Returns:
    df_join (pd.DataFrame): A DataFrame containing the normalized and smoothed mutations data. 
    It has the columns "Region", "Reference", "Allele", "count", "spectra", "norm_count", and "norm_count_smooth".
    """
    from numpy import NaN

    df_new = df_feature_mutations_master_coding[["Region", "Reference", "Allele"]]
    df_new = df_new.groupby(["Region", "Reference", "Allele"]).size().reset_index(name="count")

    #inner join the two dataframes; df_new on "Region" and df_features_mod on "Distance"
    df_join = pd.merge(df_new, df_features_mod, left_on="Region", right_on="Distance", how="inner")
    df_join = df_join.drop(columns=["Region"])
    df_join['spectra'] = df_join['Reference'] +"_to_"+ df_join['Allele']
    df_join = df_join[(df_join['spectra'] == "C_to_T") | (df_join['spectra'] == "G_to_A")]
    df_join["Forward_and_Reverse"] = df_join["Forward"] + df_join["Reverse"]

    #legacy code
    # df_join['norm_count'] = df_join.apply(
    #     lambda x: 
    #         0 if (x['Forward'] == 0) or (x['Reverse'] == 0)
    #         else (x['count'] / x['Forward']) if x['spectra'] == "C_to_T" 
    #         else (x['count'] / x['Reverse']) if x['spectra'] == "G_to_A" 
    #         else 0, axis=1)

    df_join['norm_count'] = df_join.apply(
        lambda x: 
            0 if (x['Forward'] == 0) or (x['Reverse'] == 0)
            else (x['count'] / (x['Forward'] / feature_count_norm) ) if x['spectra'] == "C_to_T" 
            else (x['count'] / (x['Reverse'] / feature_count_norm) ) if x['spectra'] == "G_to_A" 
            else 0, axis=1)

    df_join = df_join.reset_index(drop=True)

    #rename Distance to Region
    df_join = df_join.rename(columns={"Distance": "Region"})

    #combine the C_to_T and G_to_A spectra by adding the counts
    df_combined_spectra = df_join.copy()
    df_combined_spectra = df_combined_spectra.groupby(["Region"]).agg({"count": "sum", "Forward_and_Reverse": "max"}).reset_index()
    df_combined_spectra['spectra'] = "C_to_T_and_G_to_A"
    df_combined_spectra['norm_count'] = df_combined_spectra.apply(
        lambda x:
            0 if (x['Forward_and_Reverse'] == 0)
            else (x['count'] / (x['Forward_and_Reverse'] / feature_count_norm) ) if x['spectra'] == "C_to_T_and_G_to_A"
            else 0, axis=1)

    df_join_ct = df_join[df_join['spectra'] == "C_to_T"].copy()
    #make sure there is a row for each position in the region range (feature_range), if not, add a row with 0 count
    missing_regions = [i for i in range(feature_range[0], feature_range[1]+1) if i not in df_join_ct['Region'].values]
    new_rows = [{"Region": i, "count": 0, "Forward": 0, "Reverse": 0, "Forward_and_Reverse": 0, "spectra": "C_to_T", "norm_count": 0} for i in missing_regions]
    df_join_ct = pd.concat([df_join_ct, pd.DataFrame(new_rows)], ignore_index=True)

    df_join_ga = df_join[df_join['spectra'] == "G_to_A"].copy()
    #make sure there is a row for each position in the region range (feature_range), if not, add a row with 0 count
    missing_regions = [i for i in range(feature_range[0], feature_range[1]+1) if i not in df_join_ga['Region'].values]
    new_rows = [{"Region": i, "count": 0, "Forward": 0, "Reverse": 0, "Forward_and_Reverse": 0, "spectra": "G_to_A", "norm_count": 0} for i in missing_regions]
    df_join_ga = pd.concat([df_join_ga, pd.DataFrame(new_rows)], ignore_index=True)

    #add missing rows to the combined spectra dataframe
    missing_regions = [i for i in range(feature_range[0], feature_range[1]+1) if i not in df_combined_spectra['Region'].values]
    new_rows = [{"Region": i, "count": 0, "Forward_and_Reverse": 0, "spectra": "C_to_T_and_G_to_A", "norm_count": 0} for i in missing_regions]
    df_combined_spectra = pd.concat([df_combined_spectra, pd.DataFrame(new_rows)], ignore_index=True)    
    
    
    df_join = pd.concat([df_join_ct, df_join_ga], ignore_index=True)
    df_join = pd.concat([df_join, df_combined_spectra], ignore_index=True)

    df_join = df_join.reset_index(drop=True)

    df_join['norm_count'] = df_join['norm_count'] / sample_size_norm
    
    df_join['norm_count'] = df_join['norm_count'] / feature_count_norm

    return df_join

def describe_pd_series(series: pd.Series, sig_decimals: int = 10) -> None:
    """
    This function prints the description of a pandas Series with a given number of significant decimals.

    Args:
        series (pd.Series): The pandas Series to describe.
        sig_decimals (int, optional): The number of significant decimals to use. Defaults to 10.

    Returns:
        None: This function doesn't return anything; it prints the description of the Series.
    """
    print(f"Count: {series.count()}")
    print(f"Mean: {series.mean():.{sig_decimals}f}")
    print(f"Std: {series.std():.{sig_decimals}f}")
    print(f"Min: {series.min()}")
    print(f"25%: {series.quantile(0.25)}")
    print(f"50%: {series.median()}")
    print(f"75%: {series.quantile(0.75)}")
    print(f"Max: {series.max()}")

    

def conduct_mann_u_for_normalized_mutations_around_feature(
        df_join, 
        region_test: tuple = (-400, 0),
        region_control: tuple = (-1200, -800)
    ) -> list[list[str]] | None:
    """
    Conducts a Mann-Whitney U test for normalized mutations around a feature.

    The test is conducted separately for C_to_T and G_to_A mutations.
    The test and control regions are defined by the user.

    Args:
        df_join (pd.DataFrame): The DataFrame containing the mutation data. It should have columns 'Region', 'spectra', and 'norm_count'.
        region_test (tuple, optional): The range of the test region. Defaults to (-400, 0).
        region_control (tuple, optional): The range of the control region. Defaults to (-1200, -800).

    Returns:
        list[list[str]] | None: A list of lists containing the results of the Mann-Whitney U test for each mutation type.
        Each inner list contains the following elements: mutation type, test region, control region, p-value, mean count in test region, mean count in control region
    """
    from scipy.stats import mannwhitneyu

    region_test_lower, region_test_upper = region_test
    region_control_lower, region_control_upper = region_control

    df_join_1 = df_join[(df_join['Region'] >= region_test_lower) & (df_join['Region'] < region_test_upper)]
    df_join_2 = df_join[(df_join['Region'] >= region_control_lower) & (df_join['Region'] < region_control_upper)]

    df_join_1_ct = df_join_1[df_join_1['spectra'] == "C_to_T"].copy()
    df_join_1_ga = df_join_1[df_join_1['spectra'] == "G_to_A"].copy()
    df_join_1_both = df_join_1[df_join_1['spectra'] == "C_to_T_and_G_to_A"].copy()
    df_join_2_ct = df_join_2[df_join_2['spectra'] == "C_to_T"].copy()
    df_join_2_ga = df_join_2[df_join_2['spectra'] == "G_to_A"].copy()
    df_join_2_both = df_join_2[df_join_2['spectra'] == "C_to_T_and_G_to_A"].copy()
    
    try:
        test_region_string = f"{region_test_lower}_to_{region_test_upper}"
        control_region_string = f"{region_control_lower}_to_{region_control_upper}"

        print("***C_to_T***")
        print(f"Test region {test_region_string}:")
        describe_pd_series(df_join_1_ct['norm_count'])
        print(f"Control region {control_region_string}:")
        describe_pd_series(df_join_2_ct['norm_count'])

        mean_ct_test = df_join_1_ct['norm_count'].mean()
        mean_ct_ctrl = df_join_2_ct['norm_count'].mean()

        stat, p_1 = mannwhitneyu(df_join_1_ct['norm_count'], df_join_2_ct['norm_count'])
        print(f'Statistics={stat}, p={p_1}')

        print("***G_to_A***")
        print(f"Test region {test_region_string}:")
        describe_pd_series(df_join_1_ga['norm_count'])
        print(f"Control region {control_region_string}:")
        describe_pd_series(df_join_2_ga['norm_count'])

        mean_ga_test = df_join_1_ga['norm_count'].mean()
        mean_ga_ctrl = df_join_2_ga['norm_count'].mean()

        stat, p_2 = mannwhitneyu(df_join_1_ga['norm_count'], df_join_2_ga['norm_count'])
        print(f'Statistics={stat}, p={p_2}')

        print("***C_to_T_and_G_to_A***")
        print(f"Test region {test_region_string}:")
        describe_pd_series(df_join_1_both['norm_count'])
        print(f"Control region {control_region_string}:")
        describe_pd_series(df_join_2_both['norm_count'])

        mean_both_test = df_join_1['norm_count'].mean()
        mean_both_ctrl = df_join_2['norm_count'].mean()

        stat, p_3 = mannwhitneyu(df_join_1['norm_count'], df_join_2['norm_count'])
        print(f'Statistics={stat}, p={p_3}')

        print("***C_to_T_vs_G_to_A***")
        print(f"Test region 1 {test_region_string}:")
        describe_pd_series(df_join_1_ct['norm_count'])
        print(f"Test region 2 {test_region_string}:")
        describe_pd_series(df_join_1_ga['norm_count'])

        mean_vs_test_1 = df_join_1_ct['norm_count'].mean()
        mean_vs_test_2 = df_join_1_ga['norm_count'].mean()

        stat, p_4 = mannwhitneyu(df_join_1_ct['norm_count'], df_join_1_ga['norm_count'])
        print(f'Statistics={stat}, p={p_4}')

        return [
            ["C_to_T",            test_region_string, control_region_string, p_1, mean_ct_test,   mean_ct_ctrl],
            ["G_to_A",            test_region_string, control_region_string, p_2, mean_ga_test,   mean_ga_ctrl], 
            ["C_to_T_and_G_to_A", test_region_string, control_region_string, p_3, mean_both_test, mean_both_ctrl],
            ["C_to_T_vs_G_to_A",  test_region_string, test_region_string,    p_4, mean_vs_test_1, mean_vs_test_2]
        ]


    except Exception as e:
        print(e)
        return None

def smoothen_normalized_mutations_around_feature(
        df_norm: pd.DataFrame, 
        rolling_window: int = 100,
        method: Literal["rolling", "savgol"] = "rolling"
        ) -> pd.DataFrame:
    
    from scipy.signal import savgol_filter

    df_norm = df_norm.copy()
    
    df_norm = df_norm.sort_values(by="Region")

    df_join_ct = df_norm[df_norm['spectra'] == "C_to_T"].copy()
    df_join_ga = df_norm[df_norm['spectra'] == "G_to_A"].copy()
    df_join_both = df_norm[df_norm['spectra'] == "C_to_T_and_G_to_A"].copy()

    print(f"Using method {method} for smoothing.")
    if method == "rolling":

        df_join_ct['norm_count_smooth_rolling'] = df_join_ct['norm_count'].rolling(window=rolling_window).mean()
        df_join_ga['norm_count_smooth_rolling'] = df_join_ga['norm_count'].rolling(window=rolling_window).mean()
        df_join_both['norm_count_smooth_rolling'] = df_join_both['norm_count'].rolling(window=rolling_window).mean()

    elif method == "savgol":

        try:
            df_join_ct['norm_count_smooth'] = savgol_filter(df_join_ct['norm_count'], 100, 3)   
        except:
            print("Error with C_to_T. Using raw data.")
            df_join_ct['norm_count_smooth'] = df_join_ct['norm_count']

        try:
            
            df_join_ga['norm_count_smooth'] = savgol_filter(df_join_ga['norm_count'], 100, 3)
        except:
            print("Error with G_to_A. Using raw data.")
            df_join_ga['norm_count_smooth'] = df_join_ga['norm_count']
            
        try:
            
            df_join_both['norm_count_smooth'] = savgol_filter(df_join_both['norm_count'], 100, 3)
        except:
            print("Error with C_to_T_and_G_to_A. Using raw data.")
            df_join_both['norm_count_smooth'] = df_join_both['norm_count']

    # print("C_to_T data:")
    # print(df_join_ct)
    # print("G_to_A data:")
    # print(df_join_ga)
    # print("C_to_T_and_G_to_A data:")
    # print(df_join_both)

    df_nomr_line = pd.concat([df_join_ct, df_join_ga, df_join_both], ignore_index=True)

    return df_nomr_line

def plot_normalized_mutations_around_feature(
        df_norm: pd.DataFrame,
        smoothing_method: Literal["rolling", "savgol"] = "rolling",
        rolling_window: int = 100,
        y_lim: Union[tuple, bool, None] = None,
        suffix: str = "",
        show: bool = False,
        y_lim_pad: float = 1.1,
        include_scatter: bool = True,
        ) -> None:
    """
    This function plots the normalized mutations around a feature start. It creates 
    a scatterplot and lineplot with the normalized counts.

    Parameters:
    df_norm (pd.DataFrame): A DataFrame containing the normalized counts of mutations. 
    It should have the columns "Region", "norm_count", "norm_count_smooth", and "spectra".
    smoothing_method (Literal["rolling", "savgol"], optional): The smoothing method to use. Defaults to "rolling".
    rolling_window (int, optional): The rolling window for smoothing. Defaults to 100.
    y_lim (Union[tuple, bool, None], optional): The y-axis limits. Defaults to None.
    suffix (str, optional): The suffix to add to the plot filename. Defaults to "".
    show (bool, optional): Whether to show the plot. Defaults to False.
    y_lim_pad (float, optional): The padding for the y-axis limits. Defaults to 1.1.

    Returns:
    None: This function doesn't return anything, it just plots the data.

    """
    from matplotlib import ticker

    #create a scatterplot +lineplot with the normalized counts
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.set_context("poster")

    #set the y-axis labels to scientific notation
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    pallette = {'C_to_T': 'tab:blue', 'G_to_A': 'tab:orange', 'C_to_T_and_G_to_A': 'tab:green'}
    hue_order = ['C_to_T', 'G_to_A', 'C_to_T_and_G_to_A']

    if include_scatter:
        sns.scatterplot(data=df_norm, x="Region", y="norm_count", hue="spectra", ax=ax, palette=pallette, hue_order=hue_order, alpha=0.5, s=10)
    
    #smooth the data (for line-plotting purposes)
    df_norm_line = smoothen_normalized_mutations_around_feature(df_norm, rolling_window=rolling_window, method=smoothing_method)

    print(df_norm_line)

    if smoothing_method == "savgol":
        sns.lineplot(data=df_norm_line, x="Region", y="norm_count_smooth", hue="spectra", ax=ax, palette=pallette, hue_order=hue_order, alpha=1, linewidth=2)

    elif smoothing_method == "rolling":
        df_norm_line['Region'] = df_norm_line['Region'] - (rolling_window // 2)
        sns.lineplot(data=df_norm_line, x="Region", y="norm_count_smooth_rolling", hue="spectra", ax=ax, palette=pallette, hue_order=hue_order, alpha=1, linewidth=2)
    
    if y_lim is None:
        #set y-lim to 110% of the max value int the lineplot data
        y_lim = df_norm_line["norm_count_smooth_rolling"].max() * y_lim_pad
        print(f"Setting y limit to {y_lim}")
        ax.set_ylim(0, y_lim)
    elif y_lim is False:
        print("No y limit set. Plotting full range.")
        pass
    else:
        ax.set_ylim(y_lim)

    ax.axvline(x=0, color='black', linestyle='--', linewidth=1) #draw a vertical line at 0
    ax.set_title("Normalized A3A Mutation Frequencies around Feature", fontsize=24, weight='bold')
    ax.set_xlabel("Distance from Feature Position (bp)", fontsize=20, weight='bold')
    ax.set_ylabel("Normalized A3A Mutation Frequency", fontsize=20, weight='bold')
    #ax.set_xlim(df_norm['Region'].min(), df_norm['Region'].max()) #cut x at lowest and highest value
    ax.xaxis.set_major_locator(ticker.MultipleLocator(200)) #set the x ticks every 100 bp
    #set x axis limit to the same as the region
    ax.set_xlim(df_norm['Region'].min(), df_norm['Region'].max())

    #relabel the legend C_to_T and G_to_A to C-to-T and G-to-A
    handles, labels = ax.get_legend_handles_labels()
    relabel_dict = {
        'C_to_T': 'Non-Transcribed (C→T)', 
        'G_to_A': 'Transcribed (G→A)',
        'C_to_T_and_G_to_A': 'Both Strands',
        'spectra': 'Mutated Strand'}
    labels = [relabel_dict[label] for label in labels]
    ax.legend(handles=handles, labels=labels, title="Mutated Strand", title_fontsize=20, fontsize=16)

    plt.tight_layout()
    plt.savefig(f"figures/transcription/Normalized_Mutations_around_Feature_Position_{suffix}.png", dpi=300)

    if show:
        plt.show()
    else:
        plt.close()

def saveClusteredSNVsToFile(
        samples_SNVs: list[OverlapSwitchesTable], 
        path: str = "data/new_clusters/clustered_mutations/clustered_SNVs.tsv"):
    
    # df_clustered_SNVs = pd.DataFrame()
    # for sample in samples_SNVs:
    #     df_swap = sample.df_SNPs.copy()
    #     df_swap["Sample"] = sample.name

    #     df_clustered_SNVs = pd.concat([df_clustered_SNVs, df_swap], ignore_index=True)
    # df_clustered_SNVs.to_csv(path, sep="\t", header=True, index=False)

    cluster_counter = 0
    df_clustered_SNVs = pd.DataFrame()
    for sample in samples_SNVs:
        df_swap = sample.df_SNPs.copy()
        df_swap["Sample"] = sample.name
        df_swap["Dataset_Cluster_ID"] = df_swap["Cluster_ID"] + cluster_counter
        df_swap["genotype"] = sample.genotype
        cluster_counter += len(sample.df_SNPs["Cluster_ID"].unique())
        df_clustered_SNVs = pd.concat([df_clustered_SNVs, df_swap], ignore_index=True)

    df_clustered_SNVs.to_csv(path, sep="\t", header=True, index=False, encoding='utf-16')

def saveScatteredSNVsToFile(
        samples_SNVs: list[OverlapSwitchesTable], 
        path: str = "data/new_clusters/scattered_mutations/scattered_SNVs.tsv"):

    df_scattered_SNVs = pd.DataFrame()
    for sample in samples_SNVs:
        df_swap = sample.df_scattered.copy()
        df_swap["Sample"] = sample.name
        df_swap["genotype"] = sample.genotype
        df_scattered_SNVs = pd.concat([df_scattered_SNVs, df_swap], ignore_index=True)

    df_scattered_SNVs.to_csv(path, sep="\t", header=True, index=False, encoding='utf-16')
    

def showRecombinationHotspots(
        switches_path: str = "outputs/all_switches_2snps.tsv", 
        chromosome: str = "ref|NC_001135|",
        show: bool = True, 
        genotype_dict: dict | None = None,
        genotype_list: list | None = None,
        show_most_common: int | None = None,
        **kwargs):
    """
    This function visualizes recombination hotspots (switches) in a given chromosome.

    Parameters:
    switches_path (str): The path to the file containing switch data. Default is "outputs/all_switches_2snps.tsv".
    chromosome (str): The chromosome to visualize. Default is "ref|NC_001135|".
    show (bool): Whether to show the plot or not. Default is True.
    genotype_dict (dict): A dictionary mapping sample names to genotypes. Default is None.
    genotype_list (list): A list of genotypes to include in the visualization. Default is None.
    show_most_common (int): The number of most common recombination hotspots to print. Default is None.
    **kwargs: Additional keyword arguments for saving the plot. If 'save_path' is provided, the plot will be saved to this location.

    Returns:
    None
    """
    #load df_switches_all from a file
    df_switches_all = pd.read_csv(switches_path, sep="\t") #all_switches.tsv

    if genotype_dict is not None:
        df_switches_all["Genotype"] = df_switches_all["Sample_Name"].map(genotype_dict)

        if genotype_list is not None:
            df_switches_all = df_switches_all[df_switches_all["Genotype"].isin(genotype_list)]

    temp = df_switches_all[df_switches_all["Percent_Parental"] > 85]
    temp = temp[temp["Chromosome"] != "ref|NC_001224|"]
    temp = temp[temp["Chromosome"].isin([chromosome])]

    if show_most_common is not None:
        n=show_most_common
        print(f"Top {n} recombination hotspots:")
        print(temp.groupby(["Position_1", "Position_2", "Chromosome"]).size().reset_index().rename(columns={0: "count"}).sort_values(by="count", ascending=False).head(n))
        print(f"Top {n} position 1:")
        print(temp.groupby(["Position_1", "Chromosome"]).size().reset_index().rename(columns={0: "count"}).sort_values(by="count", ascending=False).head(n))
        print(f"Top {n} position 2:")
        print(temp.groupby(["Position_2", "Chromosome"]).size().reset_index().rename(columns={0: "count"}).sort_values(by="count", ascending=False).head(n))

    sns.set_context(context="talk", font_scale = 1, rc={"lines.linewidth": 2})
    sns.set_style("white")
    sns.displot(data=temp, x="Switch_Center", col="Chromosome", col_wrap=2, bins=500, aspect=2.5, height=8)

    plt.xlabel("Genomic position of a recombination event")
    plt.ylabel("Number of recombination events")

    #if user specified path to save the plot, save it
    if "save_path" in kwargs:
        plt.savefig(kwargs["save_path"], dpi=300, bbox_inches='tight')

    #if user specified to show the plot, show it. Otherwise, close it
    if show:
        plt.show()
    else:
        plt.close()

def create_project_directory_tree():
    """
    This function creates a directory tree for the project.
    """
    #create directory tree
    import os

    directories = [
        'figures',
            'figures/cluster_length_violin_by_pval',
                'figures/cluster_length_violin_by_pval/JT', 
            'figures/cluster_types_by_pval',
                'figures/cluster_types_by_pval/JT',
            'figures/clusters',
            'figures/combined',
            'figures/contextFigs',
            'figures/coverage',
            'figures/mutations_around_rec_events',
            'figures/null_hypothesis',
            'figures/other',
            'figures/rec_events',
            'figures/recombination_maps',
            #'figures/resources',
            'figures/SNP_parental_frequencies',
            'figures/SNVs_A3A_all',
            'figures/SNVs_A3A_cluster',
            'figures/SNVs_switches_1snp_QQplots',
            'figures/SNVs_switches_2snp_QQplots',
            'figures/SNVs_switches_histograms_1snp',
            'figures/SNVs_switches_histograms_2snp',
            'figures/SNVs_switches_lineplots_1snp',
            'figures/SNVs_switches_lineplots_2snp',
            'figures/switches',
            'figures/tests',
            'figures/transcription',
        'outputs',
            'outputs/associations',
                'outputs/associations/JT',
            'outputs/clusters',
                'outputs/clusters/JT',
            'outputs/high_coverage_SNP_positions',
            'outputs/logs',
            'outputs/low_coverage_SNP_positions',
            'outputs/partial_pickles',
            'outputs/simulations',
            'outputs/switches',
        ]
    
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
            logging.info(f"Created directory: {directory}")

def plot_rec_events_jointly(
        df_all_rec_events: pd.DataFrame, 
        event_types: list = ["GC", "CO"],
        save: bool = True,
        show: bool = False,
        lineplot_kde: bool = False,
        output_path: str = "figures/rec_events/rec_events_across_all_chr_SNP.png"
        ) -> None:
    """
    This function plots the distribution of recombination events (GCs and COs) for all chromosomes on separate subplots.
    It applies a sliding window approach to smoothen the data and plots both the histogram of the events and the smoothened data on a joint y-axis.

    Parameters:
    - df_all_rec_events (pd.DataFrame): A DataFrame containing information about all recombination events.
    - event_types (list): A list of event types to consider. Default is ["GC", "CO"].
    - save (bool): If True, the plot is saved. Default is True.
    - show (bool): If True, the plot is displayed. If False, the plot is not displayed. Default is False.
    - output_path (str): The path where the plot is saved. Default is "figures/rec_events/rec_events_across_all_chr_SNP.png".

    Returns:
    None
    """
    from mayo.settings import S288C_to_roman
    from mayo.settings import roman_to_S288C
    from mayo.settings import chromosomes as chromosomes
    import numpy as np

    sns.reset_defaults()
    df_both_slide = pd.DataFrame()
    
    fig, axes = plt.subplots(nrows=16, ncols=1, figsize=(20, 30), sharex=True)
    plt.subplots_adjust(hspace=0.1) #set little space between subplots

    df_all_rec_events = df_all_rec_events[df_all_rec_events["Chromosome"] != "ref|NC_001224|"]

    #sort df by chromosome name descending
    df_all_rec_events = df_all_rec_events.sort_values(by="Chromosome")
    df_all_rec_events["Chromosome"] = df_all_rec_events["Chromosome"].map(S288C_to_roman)
    df_all_rec_events["Region_kb"] = df_all_rec_events["Region"] / 1000
    longest_chromosome = max([chromosome[1] for chromosome in chromosomes])
   
    for idx, chromosome_to_plot in enumerate(df_all_rec_events["Chromosome"].unique()):

        ax = axes[idx]
        df_all_rec_events_chr = df_all_rec_events[df_all_rec_events["Chromosome"] == chromosome_to_plot]
        print("Chromosome to plot:", chromosome_to_plot)
        chromosome_to_plot_S288C = roman_to_S288C[chromosome_to_plot]
        chromosome_length = find_chromosome_length(chromosome_to_plot_S288C, chromosomes)
        
        bin_size = 5000
        n_bins = int(chromosome_length / bin_size)

        for event_type in event_types:
            df_all_rec_events_chr_event = df_all_rec_events_chr[df_all_rec_events_chr["Type"] == event_type]
            df_sliding_window = smoothen_GC_CO_data(
                df=df_all_rec_events_chr_event, 
                smoothing_range=(0, chromosome_length), 
                event_type=event_type,
                window_size=5000, 
                window_slide=500)

            df_sliding_window = smoothen_derive_col_data(df_sliding_window, "Events", window_length=12)
            df_sliding_window.loc[df_sliding_window["Events_smoothed"] < 0, "Events_smoothed"] = 0
            df_sliding_window["Type"] = event_type
            df_sliding_window["Chromosome"] = chromosome_to_plot
            df_both_slide = pd.concat([df_both_slide, df_sliding_window], ignore_index=True)

        df_chr = df_both_slide[df_both_slide["Chromosome"] == chromosome_to_plot].copy()
        df_chr["Region"].astype(float)
        df_chr["Region_kb"] = df_chr["Region"] / 1000

        # Plotting
        palette = {"GC": "tab:blue", "CO": "tab:red"}

        #add chromosome name to the plot in left corner of the subplot
        ax.text(0.01, 0.9, f"Chromosome {chromosome_to_plot}", transform=ax.transAxes, fontsize=16, fontweight="bold")
        ax.set_ylabel("Number of events", fontsize=14)
        ax.set_ylim(0, 40)
        ax.set_xlim(0, longest_chromosome / 1000)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if idx != 15:
            ax.set_xticks([])
        else:
            #set x axis ticks every 50kb
            ax.set_xticks(np.arange(0, longest_chromosome / 1000, 50))

        sns.histplot(ax=ax, data=df_all_rec_events_chr, x="Region_kb", hue="Type", stat="count", element="bars", palette=palette, multiple="stack", bins=n_bins, legend=False)
        if lineplot_kde:
            sns.lineplot(ax=ax, data=df_chr, x="Region_kb", y="Events_smoothed", hue="Type", palette=palette, legend=False, alpha=0.8)

    #add x axis label to the plot
    fig.text(0.5, -0.01, "Chromosome position (kb)", ha='center', fontsize=16, fontweight="bold")
    plt.tight_layout()
    
    #add a legend for a histogram
    handles = [plt.Rectangle((0,0),1,1, color=palette[event_type]) for event_type in event_types]
    fig.legend(handles, event_types, fontsize=14)

    if save:
        plt.savefig(output_path, dpi=300)
    if show:
        plt.show()
    else:
        plt.close()

def createAllSNVsTable(samples_SNVs: list, genotypes_keep: list | None = None) -> pd.DataFrame:
    """Creates a DataFrame with all SNVs from the provided samples.

    Args:
        samples_SNVs (list): List of samples containing SNV data.
        genotypes_keep (list, optional): List of genotypes to keep. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame containing all SNVs.
    """

    df_all_mutations = pd.DataFrame()

    for sample in samples_SNVs:
        df = sample.df.copy()
        df["sample_name"] = sample.name
        df["is_heterozygous"] = sample.heterozygous
        df["is_aneuploid"] = sample.aneuploid
        df["genotype"] = sample.genotype
        df_all_mutations = pd.concat([df_all_mutations, df], ignore_index=True)

    if genotypes_keep is not None:
        df_all_mutations = df_all_mutations[df_all_mutations["genotype"].isin(genotypes_keep)]

    columns_to_drop=[
        "Forward/reverse balance", 
        "# unique start positions",
        "# unique end positions",
        'BaseQRankSum',
        'Read position test probability',
        'Read direction test probability',
        'Hyper-allelic',
        'Homopolymer',
        'Homopolymer length',
        'Control count', 
        'Control coverage',
        'Control frequency',
        "Linkage",
        "Genotype"]
    
    #if column is not in DataFrame, it will throw an error, so we need to check first
    columns_to_drop = [i for i in columns_to_drop if i in df_all_mutations.columns]
    df_all_mutations = df_all_mutations.drop(columns=columns_to_drop)

    return df_all_mutations

def createMasterSNVDataFrames(samples_SNVs, sample_names, just_heterozygous=False, just_aneuploid=False):

    df_master_scattered = pd.DataFrame()
    df_master_all = pd.DataFrame()
    df_master_clusters = pd.DataFrame()

    sample_counter = 0

    for sample in samples_SNVs:

        if sample.name not in sample_names:
            continue

        if just_heterozygous and sample.heterozygous == False:
            continue
        
        if just_aneuploid and sample.aneuploid == False:
            continue

        df_master_all = pd.concat([df_master_all, sample.df.copy()], ignore_index=True)
        df_master_scattered = pd.concat([df_master_scattered, sample.df_scattered.copy()], ignore_index=True)
        
        sample_counter += 1
        
    df_master_clusters = df_master_all.drop(df_master_scattered.index)

    print(f"Number of samples in the analysis: {sample_counter}")
    print(f"df_master_all length: {len(df_master_all)}")
    print(f"df_master_scattered length: {len(df_master_scattered)}")
    print(f"df_master_clusters length: {len(df_master_clusters)}")

    return df_master_all, df_master_scattered, df_master_clusters, sample_counter

def addGenomeReferenceSystemOffesets(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Adjusts the positions of features based on the genome reference system offsets.

    Args:
    df (pd.DataFrame): DataFrame containing genomic data with columns 'Chromosome' and 'Position'.

    Returns:
    pd.DataFrame: DataFrame with adjusted positions.
    '''

    # add the length of the reporter to the Position column if the feature is on chromosome 3 
    # and the start position is greater than the reporter start position
    reporter_start_pos_chr3 = 211813
    offset_reporter = 2723
    df["Position"] = df.apply(
        lambda row: row["Position"] + offset_reporter 
        if (row["Chromosome"] == "ref|NC_001135|") & (row["Position"] >= reporter_start_pos_chr3) 
        else row["Position"], axis=1)

    # add the length of the ura3 deletion to the Position column if the feature is on chromosome 5 
    # and the start position is greater than the ura3 end position
    # first remove the positions that are in the range of the ura3 deletion, then adjust the positions
    ura3_start_pos_chr5 = 116166
    ura3_end_pos_chr5 = 116970
    ura3_del_len = 805
    df = df[~((df["Chromosome"] == "ref|NC_001137|") & (df["Position"].between(ura3_start_pos_chr5, ura3_end_pos_chr5)))].copy()
    df["Position"] = df.apply(
        lambda row: row["Position"] - ura3_del_len 
        if (row["Chromosome"] == "ref|NC_001137|") & (row["Position"] >= ura3_end_pos_chr5) 
        else row["Position"], axis=1)

    return df

def plot_normalized_rec_A3A(
        meta_master_relative_GC_CO_df, 
        region = (-12, 12), #in kb, -8 to 8 or -12 to 12
        bins = 24, # 16 or 24
        poly_order = 3, # 1 or 3
        genotypes_include = ["ung1∆", "pol32∆", "exo1-nd", "exo1-ndpol32∆"],
        save_path: str | None = "figures/normalized_rec_A3A.png",
        show: bool = True,
        normalize_to: str = "ung1∆",
        mutation_type: Literal["clustered", "scattered", "all"] | None = None,
        palette: dict | None = None
        ) -> None:
    
    from scipy.stats import linregress
    import numpy as np

    meta_master_relative_GC_CO_df_filtered = bin_rec_relative_mutations(meta_master_relative_GC_CO_df, region, bins=bins)
    df_norm = normalize_rec_relative_bins_to_wt(meta_master_relative_GC_CO_df_filtered, normalize_to=normalize_to)

    print("Dataframe with normalized counts to ung1∆ for plotting:")
    print(df_norm)

    df_norm = df_norm[["Region_bin"] + genotypes_include] #only keep the genotypes we want to plot

    #first melt the dataframe into a long format, all genotypes in one column
    df_counts_melted = pd.melt(
        df_norm, 
        id_vars=["Region_bin"], 
        value_vars=genotypes_include, 
        var_name="genotype", 
        value_name="norm_count")
    
    #apply log2 to the normalized counts
    df_counts_melted["norm_count"] = np.log2(df_counts_melted["norm_count"])
    #if there are values that are inf or nan, remove them
    df_counts_melted = df_counts_melted.replace([np.inf, -np.inf], np.nan).dropna()


    if palette is None:
        palette = {
            "ung1∆": "tab:blue",
            "ung1∆_clustered": "tab:blue",
            "ung1∆_scattered": "tab:cyan",

            "ung1∆NAT": "tab:orange",
            "ung1∆NAT_clustered": "tab:orange",
            "ung1∆NAT_scattered": "#FFA343",

            "pol32∆": "tab:green", 
            "pol32∆_clustered": "tab:green",
            "pol32∆_scattered": "tab:olive",

            "exo1-nd": "tab:red",
            "exo1-nd_clustered": "tab:red",
            "exo1-nd_scattered": "tab:pink",

            "exo1-ndpol32∆": "tab:purple",
            "exo1-ndpol32∆_clustered": "tab:purple",
            "exo1-ndpol32∆_scattered": "tab:gray",

            "sgs1∆C": "tab:pink",
            "exo1-ndsgs1∆C": "tab:gray",
            
        }

    sns.set_context("poster")
    plt.figure(figsize=(16, 8))
    sns.set_style("whitegrid")
    sns.lineplot(data=df_counts_melted, x="Region_bin", y="norm_count", hue="genotype", legend=False, alpha=0.3, palette=palette)

    for i, genotype in enumerate(df_counts_melted["genotype"].unique().tolist()):

        if genotype == "ung1∆NAT": # or genotype == "exo1-nd"
            continue

        df_subset = df_counts_melted[df_counts_melted['genotype'] == genotype]
        
        #overlay a regression line on the scatterplot
        #sns.regplot(data=df_subset, x="Region_bin", y="norm_count", color=palette[genotype], scatter=True, ci=95, order=poly_order)
        
        #or just plot the scatterplot
        sns.scatterplot(data=df_subset, x="Region_bin", y="norm_count", color=palette[genotype], legend=False)

        slope, intercept, r_value, p_value, std_err = linregress(df_subset["Region_bin"], df_subset["norm_count"])
        print(f"{genotype} regression line: slope={slope}, intercept={intercept}, r_value={r_value}, p_value={p_value}, std_err={std_err}")

        p_value_F = test_third_order_poly_fit_F_test(df_subset)
        print(f'{genotype} poly line p-value: {p_value_F}')

    if mutation_type is not None:
        plt.title(f"{mutation_type.capitalize()} mutations relative to recombining end", fontsize=30, fontweight='bold')
    else:
        plt.title(f"C->T mutations relative to recombining end", fontsize=30, fontweight='bold')
        
    plt.ylabel("Normalized mutation count (log2)", fontsize=20, fontweight='bold')
    #plt.ylabel("Normalized mutation count")
    plt.xlabel("Relative position (kb)", fontsize=20, fontweight='bold')

    from matplotlib import lines as mlines
    # Create custom legend handles with alpha set to 0
    handles = []
    labels = []

    for genotype in df_counts_melted["genotype"].unique():
        color = palette[genotype]
        handle = mlines.Line2D([], [], color=color, label=genotype, marker='o', linestyle='None', markersize=10)
        handles.append(handle)
        labels.append(genotype)

    # Add custom legend with italic font style
    legend = plt.legend(handles=handles, labels=labels, loc='upper right', fontsize=12, prop={'style': 'italic'})

    # Set the alpha for the markers in the legend
    for legend_handle in legend.legendHandles:
        legend_handle.set_alpha(1)

    # Add custom legend with italic font style
    plt.legend(handles=handles, labels=labels, loc='upper right', fontsize=12, prop={'style': 'italic'})

    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def find_chr_mutation_rates(df_master, sample_size):
    """
    Calculate and print the mutation rates for each chromosome.

    Parameters:
    df_master (pd.DataFrame): DataFrame containing all mutation data.
    sample_size (int): Number of samples analyzed.

    Returns:
    None
    """
    from mayo.settings import S288C_to_roman
    from mayo.settings import chromosomes as chromosomes
    df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])

    df_test = df_master.copy()
    df_test = df_test[df_test["Chromosome"] != "ref|NC_001224|"]

    for chromosome in df_test["Chromosome"].unique():
        chr_length = df_chr[df_chr['chromosome'] == chromosome]['end_position'].values[0]
        total_chr_mutations = len(df_test[df_test['Chromosome'] == chromosome])
        mutation_rate_pre_norm = total_chr_mutations / (chr_length * 1000)
        mutation_rate_norm = mutation_rate_pre_norm / sample_size
        chr_in_roman = S288C_to_roman[chromosome]
        logging.info(f"Chr {chr_in_roman} ({chromosome}) mutation rate: {mutation_rate_norm} SNVs/kb/sample.")

def plot_nucleotide_frequencies_around_feature(
        df_nucleotide_freq,
        plot_limits: tuple = (-500, 500),
        save_name: str | None = "nucleotide_frequencies_around_feature.png",
        title: str = "Nucleotide Frequencies around Feature Position",
        show: bool = False
        ):

    #prepare the dataframe for plotting
    df_plot = df_nucleotide_freq.copy()
    df_length = len(df_plot)
    flank_length = (df_length -1) / 2
    if not flank_length.is_integer(): #make sure that flank_length is an integer, if it is a fraction, throw an error
        raise ValueError("The length of the dataframe is not an odd number, cannot calculate the flank length")
    df_plot.reset_index(inplace=True) #move index to a column
    df_plot.rename(columns={"index": "relative_position"}, inplace=True)
    df_plot["relative_position"] = df_plot["relative_position"] - flank_length
    df_plot = df_plot[(df_plot["relative_position"] >= plot_limits[0]) & (df_plot["relative_position"] <= plot_limits[1])]
    df_plot.set_index("relative_position", inplace=True)


    #draw the plot
    fig, ax = plt.subplots(figsize=(24, 12))
    sns.set_context("poster")
    df_plot.plot(ax=ax)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=2)
    ax.set_title(title, fontsize=40, weight='bold')
    ax.set_xlabel("Position", fontsize=32, weight='bold')
    ax.set_ylabel("Frequency", fontsize=32, weight='bold')
    plt.tight_layout()

    if show:
        plt.show()

    if save_name is not None:
        plt.savefig(save_name, dpi=300)

    plt.close()

def calculateFeatureMutationRates(
        df_mutations: pd.DataFrame,
        sample_size_norm: int,
        set_name: str = "",
        query_position: Literal["Start", "End"] = "Start",
        pk_features_path = "data/orf_coding.fasta",
        rna_features_path = "data/rna_coding.fasta",
        just_non_essential_genes: list | None = None,
        tss_and_tts_from_utrs: bool = True,
        UTR_data_path: str = "data/yeastmine_results_2024-05-24T10-03-55.csv",
        smoothing_method: Literal["rolling", "savgol"] = "rolling",
        rolling_window: int = 100,
        region_control: tuple[int, int] = (600, 1200),
    ):
    
    df_feature_mutations_master_coding, df_features = create_mutations_around_features_start_df(
        df_mutations=df_mutations, 
        features_path=pk_features_path, #, "data/orf_coding.fasta", #"data/other_features_genomic.fasta", #,
        query_position=query_position,
        window_size = 1500,
        just_non_essential_genes=just_non_essential_genes,
        tss_and_tts_from_utrs=tss_and_tts_from_utrs, #True,
        UTR_data_path=UTR_data_path #"data/yeastmine_results_2024-05-24T10-03-55.csv"
        )
    pk_df_nucleotide_freq = get_nucleotide_frequencies_around_feature(
        df_features, 
        context_flank=1200, 
        query_position=query_position,
        fasta_file_path='S288C_reference_sequencewoUra3wUra329wNAT.fa'
        )
    
    plot_nucleotide_frequencies_around_feature(pk_df_nucleotide_freq, save_name=f"figures/nucleotide_frequencies_around_feature_{query_position}_protein_coding.png", title=f"Nucleotide Frequencies around Feature {query_position} (protein coding)")

    y_lim = None # None (0, 0.000075) (0, 0.001)

    plot_mutations_around_feature(df_feature_mutations_master_coding, title=f"Mutations around Feature {query_position} (protein coding)", window_size=1200)
    plot_mutations_around_feature(df_feature_mutations_master_coding, title=f"Mutations around Feature {query_position} (protein coding)", window_size=500)
    plot_A3A_motifs_relative_to_position(countA3A_motifs_relative_to_position(df_features, context_flank=200, query_position=query_position), suffix=f"orf_coding_200_{query_position}_{set_name}")
    df_features_mod = countA3A_motifs_relative_to_position(df_features, context_flank=1200, query_position=query_position)
    plot_A3A_motifs_relative_to_position(df_features_mod, suffix=f"orf_coding_1200_{query_position}_{set_name}")
    df_norm = normalize_mutations_around_feature(df_feature_mutations_master_coding, df_features_mod, feature_range=(-1200, 1200), sample_size_norm=sample_size_norm, feature_count_norm=4410)
    plot_normalized_mutations_around_feature(df_norm, suffix=f"orf_coding_1200_{query_position}_{set_name}", smoothing_method=smoothing_method, rolling_window=rolling_window, y_lim=y_lim)
    plot_normalized_mutations_around_feature(df_norm[(df_norm['Region'] > -500) & (df_norm['Region'] < 500)], suffix=f"orf_coding_500_{query_position}_{set_name}", smoothing_method=smoothing_method, rolling_window=rolling_window, y_lim=y_lim)
    
    result_list_pk = []
    result_list_pk.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(-25, 25), region_control=(-1200, -900)))
    result_list_pk.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(-25, 25), region_control=region_control))
    result_list_pk.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(-50, 50), region_control=region_control))
    result_list_pk.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(-100, 100), region_control=region_control))
    result_list_pk.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(-150, 150), region_control=region_control))

    ### RNA Coding
    df_feature_mutations_master_rna, df_features_rna = create_mutations_around_features_start_df(
        df_mutations = df_mutations, 
        features_path = rna_features_path, #"data/rna_coding.fasta",
        query_position = query_position,
        window_size = 1500,
        just_tRNA_genes = True,
        tss_and_tts_from_utrs=tss_and_tts_from_utrs, #True,
        UTR_data_path=UTR_data_path #"data/yeastmine_results_2024-05-24T10-03-55.csv"
        )
    trna_df_nucleotide_freq = get_nucleotide_frequencies_around_feature(
        df_features_rna, 
        context_flank=1200, 
        query_position=query_position,
        fasta_file_path='S288C_reference_sequencewoUra3wUra329wNAT.fa')
    
    plot_nucleotide_frequencies_around_feature(trna_df_nucleotide_freq, save_name=f"figures/nucleotide_frequencies_around_feature_{query_position}_trna.png", title=f"Nucleotide Frequencies around Feature {query_position} (tRNA)")

    y_lim = None # None (0, .6)

    plot_mutations_around_feature(df_feature_mutations_master_rna, title=f"Mutations around Feature {query_position} (rna coding)", window_size=1200)
    plot_mutations_around_feature(df_feature_mutations_master_rna, title=f"Mutations around Feature {query_position} (rna coding)", window_size=500)
    plot_A3A_motifs_relative_to_position(countA3A_motifs_relative_to_position(df_features_rna, context_flank=200, query_position=query_position), suffix=f"rna_coding_200_{query_position}_{set_name}")
    df_features_mod_rna = countA3A_motifs_relative_to_position(df_features_rna, context_flank=1200, query_position=query_position)
    plot_A3A_motifs_relative_to_position(df_features_mod_rna, suffix=f"rna_coding_1200_{query_position}_{set_name}")
    df_norm = normalize_mutations_around_feature(df_feature_mutations_master_rna, df_features_mod_rna, feature_range=(-1200, 1200), sample_size_norm=sample_size_norm, feature_count_norm=275)
    plot_normalized_mutations_around_feature(df_norm, suffix=f"rna_coding_1200_{query_position}_{set_name}", smoothing_method=smoothing_method, rolling_window=rolling_window, y_lim=y_lim)
    plot_normalized_mutations_around_feature(df_norm[(df_norm['Region'] > -500) & (df_norm['Region'] < 500)], suffix=f"rna_coding_500_{query_position}_{set_name}", smoothing_method=smoothing_method, rolling_window=rolling_window, y_lim=y_lim)
    
    result_list_trna = []
    result_list_trna.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(0, 100), region_control=region_control))
    result_list_trna.append(conduct_mann_u_for_normalized_mutations_around_feature(df_norm, region_test=(0, 100), region_control=(-1200, -900)))

    return result_list_pk, result_list_trna

def RecEventsTable_from_GC_CO_dict(GC_dict: dict, CO_dict: dict) -> pd.DataFrame:
    
    # Create lists to store the data
    gc_data = [{"Chromosome": chr, "Region": pos, "Type": "GC"} for chr in GC_dict.keys() for pos in GC_dict[chr]]
    co_data = [{"Chromosome": chr, "Region": pos, "Type": "CO"} for chr in CO_dict.keys() for pos in CO_dict[chr]]

    # Convert lists to DataFrames
    df_recEvents_GC = pd.DataFrame(gc_data)
    df_recEvents_CO = pd.DataFrame(co_data)

    # Concatenate the DataFrames
    df_recEvents = pd.concat([df_recEvents_GC, df_recEvents_CO], ignore_index=True)

    return df_recEvents

def add_GC_CO_legend_to_ax(ax, cmap_dict: dict) -> None:
    import matplotlib.lines as mlines
    
    legend_handles = []
    for rec_ev in cmap_dict:
        color = cmap_dict[rec_ev]
        marker = mlines.Line2D([], [], color=color, marker=11, linestyle='None', markersize=8, label=rec_ev)
        legend_handles.append(marker)
    
    new_legend = ax.legend(handles=legend_handles, title='Rec. Event', loc='center right', fontsize=10, title_fontsize=12, bbox_to_anchor=(0.9, .5))
    ax.add_artist(new_legend)

def calculate_similarity_score(df1, df2):
    # Step 1: Inner join on the "Position" column to find matching positions

    merged_df = pd.merge(df1, df2, on=['Chromosome', 'Region'], suffixes=('_df1', '_df2'))
    
    # Step 2: Compare the "Value" columns for matching positions
    matches = merged_df['Parent_df1'] == merged_df['Parent_df2']
    
    # Step 3: Calculate the similarity score
    similarity_score = matches.sum() / len(matches)
    
    return similarity_score

def create_similarity_matrix(sample_list):
    import numpy as np
    #create a similarity matrix for all parental samples, store in a numpy array
    similarity_matrix = np.zeros((len(sample_list), len(sample_list)))

    for iA in range(len(sample_list)):
        for iB in range(iA, len(sample_list)):
            if iA != iB:
                df1 = sample_list[iA].output_df.copy()
                df2 = sample_list[iB].output_df.copy()
                similarity_score = calculate_similarity_score(df1, df2)
                similarity_matrix[iA, iB] = similarity_score
                similarity_matrix[iB, iA] = similarity_score
            else:
                similarity_matrix[iA, iB] = 1.0  # Similarity with itself is 1

    print("finished calculating similarity matrix")
    return similarity_matrix


