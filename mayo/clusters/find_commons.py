# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
This file contains functions used to clean and prepare mutation data (CLC) for downstream analysis in P-MACD.
'''
import os
import pandas as pd

def load_data(path: str) -> pd.DataFrame:
    """
    Load all CSV files in the given directory and concatenate them into a single DataFrame.
    
    Args:
        path (str): The path to the directory containing the CSV files.
        
    Returns:
        pd.DataFrame: A DataFrame containing the concatenated data from all CSV files.
    """
    master_df = pd.DataFrame()
    for file in os.listdir(path):
        if file.endswith(".csv"):
            df = pd.read_csv(os.path.join(path, file))
            df["file_name"] = file[:8]
            master_df = pd.concat([master_df, df], ignore_index=True)
    return master_df

def find_commons(master_df: pd.DataFrame) -> pd.DataFrame:
    """
    Find all rows in the given DataFrame where Chromosome, Region, and Allele are the same.
    Count how many times each combination of values is duplicated and sort the results by count.
    
    Args:
        master_df (pd.DataFrame): The DataFrame to search for duplicates.
        
    Returns:
        pd.DataFrame: A DataFrame containing the unique combinations of Chromosome, Region, and Allele,
        sorted by the number of times each combination appears in the original DataFrame.
    """  
    df_duplicates = master_df[master_df.duplicated(["Chromosome", "Region", "Allele"], keep=False)]
    df_duplicates = df_duplicates.sort_values(by=["Chromosome", "Region", "Allele"])
    df_duplicates["count"] = df_duplicates.groupby(["Chromosome", "Region", "Allele"])["Chromosome"].transform("count")
    df_duplicates = df_duplicates.drop_duplicates(subset=["Chromosome", "Region", "Allele"], keep="first")
    df_duplicates = df_duplicates.sort_values(by=["count"], ascending=False)

    return df_duplicates

def plot_commons(df_duplicates: pd.DataFrame, x_cutoff: int) -> None:
    """
    Plot a histogram of the number of times each unique combination of Chromosome, Region, and Allele
    appears in the given DataFrame. Highlight the values that appear more than x_cutoff times.
    
    Args:
        df_duplicates (pd.DataFrame): The DataFrame containing the unique combinations of Chromosome,
            Region, and Allele, sorted by the number of times each combination appears in the original DataFrame.
        x_cutoff (int): The minimum number of times a combination must appear to be highlighted in the plot.
        
    Returns:
        None
    """
    import matplotlib.pyplot as plt
    df_duplicates["count"].plot.hist(bins=100)
    plt.xlabel("Number of Duplicates")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.axvline(x=x_cutoff, color="red")
    plt.axvspan(xmin=x_cutoff, xmax=df_duplicates["count"].max(), color="red", alpha=0.2)
    plt.show()

def remove_common_df_rows(df_in:pd.DataFrame, df_duplicates: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows from the input DataFrame that have the same Chromosome, Region, and Allele values as
    any row in the duplicates DataFrame.
    
    Args:
        df_in (pd.DataFrame): The input DataFrame to remove rows from.
        df_duplicates (pd.DataFrame): The DataFrame containing the unique combinations of Chromosome,
            Region, and Allele, sorted by the number of times each combination appears in the original DataFrame.
        
    Returns:
        pd.DataFrame: A new DataFrame with the rows that have the same Chromosome, Region, and Allele values
        as any row in the duplicates DataFrame removed.
    """
    df_in["is_duplicate"] = df_in.apply(lambda row: row["Chromosome"] in df_duplicates["Chromosome"] and row["Region"] in df_duplicates["Region"] and row["Allele"] in df_duplicates["Allele"], axis=1)
    df_in = df_in[df_in["is_duplicate"] == False]
    df_in = df_in.drop(columns=["is_duplicate"])
    return df_in

def renameChrCol(row: pd.Series) -> str:
    """
    Renames chromosome names in a pandas DataFrame.

    Args:
        row (pd.Series): A row of a pandas DataFrame containing a 'Chromosome' column.

    Returns:
        str: The renamed chromosome name.
    """
    chromosomes = {'ref|NC_001133|': 'chr1',
                'ref|NC_001134|':'chr2',
                'ref|NC_001135|':'chr3',
                'ref|NC_001136|':'chr4',
                'ref|NC_001137|':'chr5',
                'ref|NC_001138|':'chr6',
                'ref|NC_001139|':'chr7',
                'ref|NC_001140|':'chr8',
                'ref|NC_001141|':'chr9',
                'ref|NC_001142|':'chr10',
                'ref|NC_001143|':'chr11',
                'ref|NC_001144|':'chr12',
                'ref|NC_001145|':'chr13',
                'ref|NC_001146|':'chr14',
                'ref|NC_001147|':'chr15',
                'ref|NC_001148|':'chr16'}

    chromosome=row['Chromosome']
    return(chromosomes[chromosome])

def modFrame(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Modifies the input DataFrame to prepare it for downstream analysis (P-MACD).

    Args:
        frame (pd.DataFrame): A pandas DataFrame containing mutation data.

    Returns:
        pd.DataFrame: A modified pandas DataFrame with the following columns:
            - Tumor_Sample_Barcode
            - Chromosome
            - Start_position
            - Reference_Allele
            - Tumor_Seq_Allele2
            - Variant_Type

    Raises:
        None
    """
    frame=frame[['Chromosome', 'Region', 'Type', 'Reference', 'Allele', 'file_name']]
    frame = frame.rename(columns={"file_name": "Tumor_Sample_Barcode"})
    frame['Mutation_Set'] = 'CombinedSet'
    frame = frame.rename(columns={"Region": "Start_position", "Type": "Variant_Type",'Reference':'Reference_Allele','Allele':'Tumor_Seq_Allele2'})[['Tumor_Sample_Barcode','Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type']]
    frame = frame[frame.Variant_Type =='SNV']
    frame = frame[frame.Chromosome !='ref|NC_001224|']
    frame['Variant_Type']='SNP'
    frame['Chromosome'] = frame.apply(lambda row: renameChrCol(row), axis=1)

    return(frame)

def create_clean_pmacd_input(path: str, max_duplicates: int, save_master_file: str, save_duplicates_file: str) -> None:
    """
    Cleans and prepares a mutation data file for downstream analysis.

    Args:
        path (str): The path to the input mutation data file.
        max_duplicates (int): The maximum number of duplicates allowed for a mutation to be considered common.
        save_master_file (str): The path to save the cleaned mutation data file.
        save_duplicates_file (str): The path to save the file containing common mutations.

    Returns:
        None

    Raises:
        None
    """
    master_df = load_data(path)
    master_df = master_df[master_df["Type"] == "SNV"]
    master_df = master_df[master_df["Chromosome"] != "ref|NC_001224|"]
    df_duplicates = find_commons(master_df)
    plot_commons(df_duplicates, x_cutoff=max_duplicates)
    len_before = len(master_df)
    df_duplicates = df_duplicates[df_duplicates["count"] >= max_duplicates]
    len_after = len(df_duplicates)
    print(f"Removed {len_before - len_after} rows by removing duplicates with count >= {max_duplicates}")
    print(df_duplicates)

    master_df = remove_common_df_rows(master_df, df_duplicates)

    #remove poor quality data? What is poor quality data?

    mod_master_df = modFrame(master_df)
    mod_master_df.to_csv(f"{save_master_file}", index=False, header=True, sep='\t')

    #savethe df_duplicates as a CSV file
    df_duplicates.to_csv(f"{save_duplicates_file}", index=False)
