import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_chromosome_roman_name(chr_name: str) -> str:
    """
    This function takes a chromosome name in the format 'ref|NC_001133|' and returns its corresponding Roman numeral representation.
    
    Parameters:
    chr_name (str): The name of the chromosome in the format 'ref|NC_001133|'.
    
    Returns:
    str: The Roman numeral representation of the chromosome name.
    """
    
    chr_name_dict = {
        'ref|NC_001133|': 'I',
        'ref|NC_001134|': 'II',
        'ref|NC_001135|': 'III',
        'ref|NC_001136|': 'IV',
        'ref|NC_001137|': 'V',
        'ref|NC_001138|': 'VI',
        'ref|NC_001139|': 'VII',
        'ref|NC_001140|': 'VIII',
        'ref|NC_001141|': 'IX',
        'ref|NC_001142|': 'X',
        'ref|NC_001143|': 'XI',
        'ref|NC_001144|': 'XII',
        'ref|NC_001145|': 'XIII',
        'ref|NC_001146|': 'XIV',
        'ref|NC_001147|': 'XV',
        'ref|NC_001148|': 'XVI',
        'ref|NC_001224|': 'Mito'
    }

    return chr_name_dict[chr_name]

def merge_sample_SNPs_with_parental_origin_df(df_reference: pd.DataFrame, df_progeny: pd.DataFrame) -> pd.DataFrame:
    """
    Merges a dataframe of progeny SNPs with a reference dataframe based on Chromosome, Region, and Reference columns.

    Parameters:
    df_reference (pandas.DataFrame): The reference dataframe, containing Chromosome, Region, Reference, and Parent columns.
    df_progeny (pandas.DataFrame): The progeny dataframe, containing Chromosome, Region, and Reference columns.

    Returns:
    pandas.DataFrame: The merged dataframe, containing only the Chromosome, Region, Reference, and Parent columns.
    """
    #merge the two dataframes on the Chromosome, Region, and Reference columns
    df_merged = pd.merge(df_progeny, df_reference, on=['Chromosome', 'Region', 'Reference'])

    #keep only the Chromosome, Region, Reference columns
    df_merged = df_merged[['Chromosome', 'Region', 'Reference', 'Parent']]

    return df_merged

def get_snp_frequency_by_parent_df(df_merged: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the frequency of each parental SNP at each position in the merged dataframe.

    Parameters:
    df_merged (pandas.DataFrame): The merged dataframe, containing Chromosome, Region, Reference, and Parent columns.

    Returns:
    pandas.DataFrame: The grouped dataframe, containing Chromosome, Region, Reference, Parent, Count, Position_Count, and Frequency columns.
    """
    #group by Chromosome, Region, and Reference columns and count the number of times each combination occurs
    df_grouped = df_merged.groupby(['Chromosome', 'Region', 'Reference', 'Parent']).size().reset_index(name='Count')

    #get the number of times each Chromosome, Region combination occurs and store in a new column Position_Count
    df_grouped['Position_Count'] = df_grouped.groupby(['Chromosome', 'Region'])['Count'].transform('sum')

    #calculate the frequency of each parental snp at each position
    df_grouped['Frequency'] = df_grouped['Count']/df_grouped['Position_Count']

    return df_grouped

def plot_distribution_of_snp_frequencies(df_grouped: pd.DataFrame, show: bool = False) -> None:
    """
    Plots the distribution of the frequencies of each parental SNP.
    """
    #plot the distribution of occurences of each frequency
    sns.set(style="whitegrid")
    sns.distplot(df_grouped['Frequency'], kde=False, rug=True)
    plt.title('Distribution of Frequencies')
    if show:
        plt.show()
    plt.close()

def get_imbalanced_snps(df_grouped: pd.DataFrame) -> pd.DataFrame:
    """
    This function takes a DataFrame of grouped SNPs (Single Nucleotide Polymorphisms) and returns a DataFrame of imbalanced SNPs.
    An imbalanced SNP is defined as one where the frequency is either greater than 0.8 or less than 0.2.

    Parameters:
    df_grouped (pd.DataFrame): A DataFrame containing grouped SNPs. It is expected to have a 'Frequency' column.

    Returns:
    pd.DataFrame: A DataFrame containing imbalanced SNPs, sorted by 'Frequency' in descending order.
    """
    df_imbalanced = df_grouped.loc[(df_grouped['Frequency'] > .8) | (df_grouped['Frequency'] < .2)]
    df_imbalanced_top = df_imbalanced[df_imbalanced['Frequency'] > .8]
    df_imbalanced_bottom = df_imbalanced[df_imbalanced['Frequency'] < .2]
    df_imbalanced = pd.concat([df_imbalanced_top, df_imbalanced_bottom])
    df_imbalanced.sort_values(by=['Frequency'], ascending=False, inplace=True)

    return df_grouped

def plot_snp_frequencies_by_chromosome(
        df_grouped: pd.DataFrame, 
        save: bool = False, 
        save_path_root: str = 'SNP_parental_frequencies',
        show: bool = False
        ) -> None:
    """
    This function takes a DataFrame of grouped SNPs and returns a DataFrame of imbalanced SNPs.
    An imbalanced SNP is defined as one where the frequency is either greater than 0.8 or less than 0.2.

    Parameters:
    df_grouped (pd.DataFrame): A DataFrame containing grouped SNPs. It is expected to have a 'Frequency' column.

    Returns:
    pd.DataFrame: A DataFrame containing imbalanced SNPs, sorted by 'Frequency' in descending order.
    """
    #plot the frequency of each parental snp at each position on a scatter plot

    for chromosome in df_grouped['Chromosome'].unique():
        df_chromosome = df_grouped.loc[df_grouped['Chromosome'] == chromosome]

        #sort df_chromosome by frequency, from highest to lowest
        print(f'Chromosome {chromosome} ({get_chromosome_roman_name(chromosome)})')
        df_chromosome_sorted = df_chromosome.sort_values(by=['Frequency'], ascending=False)
        print(df_chromosome_sorted.head(20))

        #set figure size
        plt.figure(figsize=(12, 5))
        sns.set(style="whitegrid")
        #draw a horizontal line at .5
        plt.axhline(y=.5, color='black', linestyle='--')
        sns.scatterplot(x="Region", y="Frequency", hue="Parent", data=df_chromosome, palette=['red', 'blue'], linewidth=0, alpha=0.5)
        plt.ylim(0, 1)
        chr_roman = get_chromosome_roman_name(chromosome)
        plt.title(f'Chromosome {chr_roman} ({chromosome})')

        #save the plot to a file
        if save == True:
            plt.savefig(f'{save_path_root}_{chr_roman}.png', dpi=300)
            print(f'Plot saved to {save_path_root}_{chr_roman}.png')
        if show:
            plt.show()
        plt.close()

def find_imbalanced_snps(
        genotype_dict: dict,
        df_reference_path: str,
        df_progeny_path: str,
        genotype_list: list,
        by_genotype: bool = False,
        show: bool = False,
        save: bool = False,
        df_save_path: str | None = 'outputs/imbalanced_snps.csv'
        ) -> pd.DataFrame:
    """
    This function finds imbalanced SNPs in a given dataset. Imbalanced SNPs are defined as SNPs where the frequency of one parental allele is either greater than 0.8 or less than 0.2.

    Parameters:
    genotype_dict (dict): A dictionary mapping sample names to genotypes.
    df_reference_path (str): The path to the reference dataframe file.
    df_progeny_path (str): The path to the progeny dataframe file.
    genotype_list (list): A list of genotypes to consider.
    by_genotype (bool, optional): If True, the function will find imbalanced SNPs for each genotype in the genotype_list. If False, the function will find imbalanced SNPs for all genotypes together. Defaults to False.
    show (bool, optional): If True, the function will display plots of SNP frequencies. Defaults to False.
    save (bool, optional): If True, the function will save plots of SNP frequencies. Defaults to False.

    Returns:
    pd.DataFrame: A dataframe containing information about the imbalanced SNPs.
    """
    
    df_reference = pd.read_csv(df_reference_path, sep='\t', header=0)
    df_progeny = pd.read_csv(df_progeny_path, sep=',', header=0)

    df_progeny.rename(columns={'Allele': 'Reference'}, inplace=True)
    df_progeny["Genotype"] = df_progeny["sample_name"].map(genotype_dict)

    df_imbalanced_all = pd.DataFrame()

    if by_genotype:
        for genotype in genotype_list:

            df_progeny_genotype = df_progeny.copy()
            df_progeny_genotype = df_progeny_genotype[df_progeny_genotype['Genotype'] == genotype]

            df_merged = merge_sample_SNPs_with_parental_origin_df(df_reference, df_progeny_genotype)
            df_grouped = get_snp_frequency_by_parent_df(df_merged)

            df_imbalanced = get_imbalanced_snps(df_grouped)
            df_imbalanced["Genotype"] = genotype
            df_imbalanced_all = pd.concat([df_imbalanced_all, df_imbalanced])

            #plot_distribution_of_snp_frequencies(df_grouped, show=show)
            plot_snp_frequencies_by_chromosome(df_grouped, save=save, show=show, save_path_root=f'figures/SNP_parental_frequencies/{genotype}_SNP_parental_frequencies')

    if not by_genotype:
            
        df_merged = merge_sample_SNPs_with_parental_origin_df(df_reference, df_progeny)
        df_grouped = get_snp_frequency_by_parent_df(df_merged)

        df_imbalanced = get_imbalanced_snps(df_grouped)
        df_imbalanced_all = pd.concat([df_imbalanced_all, df_imbalanced])

        #plot_distribution_of_snp_frequencies(df_grouped, show=show)
        plot_snp_frequencies_by_chromosome(df_grouped, save=save, show=show, save_path_root=f'figures/SNP_parental_frequencies/SNP_parental_frequencies')

    #save the imbalanced SNPs to a file
    if df_save_path != None:
        df_imbalanced_all.to_csv(df_save_path, index=False)

    return df_imbalanced_all


if __name__ == '__main__':

    from settings import config as cfg
    from settings import genotype_dict as genotype_dict_master
    genotype_dict = genotype_dict_master[cfg["genotype_dict"]]

    genotype_list = [
    # "ung1∆",
    # "UNG1",
    # "ung1∆ non-selected",
    # "exo1-nd",
    # "pol32∆",
    # "exo1-ndpol32∆",
    "ung1∆NAT",
    # "ung1∆ EV",
    ]

    df_imbalanced_all = find_imbalanced_snps(
        genotype_dict=genotype_dict,
        df_reference_path='outputs/all_parental_SNPs.tsv',
        df_progeny_path='outputs/all_parental_snps_combined.csv',
        genotype_list=genotype_list,
        by_genotype=True,
        show=True,
        save=False
    )

    #save with utf-16 encoding
    #df_imbalanced_all.to_csv('outputs/imbalanced_snps.csv', index=False, encoding='utf-16', sep='\t')
    