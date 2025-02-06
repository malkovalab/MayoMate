# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
This is an independent script that can be used to generate a SNP density map from
a CSV file containing SNP data. It can also be used to plot the fraction of each 
chromosome covered by SNP-dense regions as a function of the maximum distance between SNPs.
'''

import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s:.%(funcName)s: %(message)s')
logger = logging.getLogger(__name__)
#logger = logging.getLogger('__main__.' + __name__)

def get_dens_regions(SNP_df: pd.DataFrame, max_distance: int = 500) -> dict:
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
    from settings.chromosome_lengths import chromosomes
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

def plot_chr_SNP_coverage_by_SNP_distance(
        df: pd.DataFrame, 
        distance_max: int = 2000, 
        distace_step: int = 10, 
        inverse: bool = False, 
        save_path: str = 'outputs/SNP_density.tsv', 
        fig_path: str = 'figures/SNP_density.png'
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
    from settings.chromosome_lengths import chromosomes

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
    
    #rename chr columns to roman numerals
    from settings import S288C_to_roman
    df_density['chr'] = df_density['chr'].map(S288C_to_roman)

    #plot the curve max_distance (x) vs density (y) for each chromosome
    plt.figure(figsize=(10, 6))
    sns.set_style("darkgrid")
    sns.lineplot(data=df_density, x='max_distance', y='density', hue='chr')
    plt.title("SNP coverage by distance", fontsize=16)
    plt.xlabel("Max distance between SNPs", fontsize=12)

    plt.ylabel("Fraction of chromosome covered by SNP-dense regions", fontsize=12)
    if inverse:
        plt.ylabel("Fraction of chromosome covered by SNP-poor regions", fontsize=12)

    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.subplots_adjust(right=0.80) #make room for the legend
    plt.savefig(fig_path, dpi=300)
    logging.info(f"Saved SNP density plot to {fig_path}")
    plt.show()

    #save the data
    df_density.to_csv(save_path, sep='\t', index=False)
    logging.info(f"Saved SNP density data to {save_path}")

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

if __name__ == "__main__":

    lowDensityMap = loadDensityMap(lowDensity=True, max_distance=500)

    data_path = 'outputs/all_parental_SNPs.tsv'
    df = pd.read_csv(data_path, sep='\t', header=0)
    df = df[df['Parent'] == "parent_1"]

    print(df['Chromosome'].value_counts())
    print(df['Zygosity'].value_counts())
    print(df['Parent'].value_counts())
    print(df[df['Zygosity'] != "Homozygous"])

    plot_chr_SNP_coverage_by_SNP_distance(
        df, 
        distance_max=5000, 
        distace_step=5,
        inverse=True,
        save_path='outputs/SNP_density_poor.tsv',
        fig_path='figures/SNP_density_poor.png')
