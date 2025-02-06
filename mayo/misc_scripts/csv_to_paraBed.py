# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
'''
This script processes CLC SNV output tables (csv files) so that they can be used in the P-MACD program.
The input is a directory containing csv files and the output is a single txt file containing all the data from the csv files.
The txt file is formatted so that it can be used as input for the P-MACD program.
The script:
    renames the chromosome names to integers.
   filters the data so that only SNVs are included in the output file.
   filters out the mitochondrial chromosome (ref|NC_001224|) from the output file.
   adds a 'Variant_Type' column to the output file and sets it to 'SNP' for all rows.
   adds a 'Tumor_Sample_Barcode' column to the output file and sets it to the sample name for all rows.
   (currently off) adds a 'Mutation_Set' column to the output file and sets it to the sample name for all rows.
   adds a 'Start_position' and 'End_position' columns to the output file and sets it to the 'Region' column for all rows.
   adds a 'Reference_Allele' column to the output file and sets it to the 'Reference' column for all rows.
   adds a 'Tumor_Seq_Allele2' column to the output file and sets it to the 'Allele' column for all rows.
   saves the output file as a tab-delimited txt file.

The script can be run from the command line using the following command:
   python csv_to_paraBed.py /path/to/csv/files/ output_file_name.txt

IMPORTANT:
   The script assumes that the sample name is the first 7 characters of the csv file name.
   If this is not the case, you will need to modify the script to change the STRING_LENGTH variable to the correct number of characters below.
'''

STRING_LENGTH=7

import os
import argparse
import logging
import pandas as pd

def modFrame(
        frame: pd.DataFrame, 
        sampleName: str = 'sample_name', 
        barcode: str = 'barcode'
        ) -> pd.DataFrame:
    """
    Modifies the input DataFrame by selecting specific columns, renaming them, and applying filters.

    Parameters:
    frame (pd.DataFrame): The input DataFrame to be modified.
    sampleName (str, optional): The sample name. Defaults to 'sample_name'.
    barcode (str, optional): The barcode. Defaults to 'barcode'.

    Returns:
    pd.DataFrame: The modified DataFrame with selected and renamed columns, and applied filters.
    """
    frame = frame[['Chromosome', 'Region', 'Type', 'Reference', 'Allele']].copy()
    frame['Tumor_Sample_Barcode'] = barcode
    #frame['Mutation_Set'] = sampleName
    frame = frame.rename(columns={"Region": "Start_position", "Type": "Variant_Type",'Reference':'Reference_Allele','Allele':'Tumor_Seq_Allele2'})[['Tumor_Sample_Barcode','Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type']]
    frame = frame[frame.Variant_Type =='SNV']
    frame = frame[frame.Chromosome !='ref|NC_001224|']
    frame['Variant_Type'] = 'SNP'
    frame['Chromosome'] = frame.apply(lambda row: renameChrCol(row), axis=1)

    return(frame)

def renameChrCol(row: pd.Series) -> int:
    """
    Renames the 'Chromosome' column in the input Series by mapping the chromosome names to integers.

    Parameters:
    row (pd.Series): The input Series with a 'Chromosome' column to be renamed.

    Returns:
    int: The integer representation of the chromosome name.
    """
    chromosomes_dic = { 'ref|NC_001133|' : 1,
                        'ref|NC_001134|' : 2,
                        'ref|NC_001135|' : 3,
                        'ref|NC_001136|' : 4,
                        'ref|NC_001137|' : 5,
                        'ref|NC_001138|' : 6,
                        'ref|NC_001139|' : 7,
                        'ref|NC_001140|' : 8,
                        'ref|NC_001141|' : 9,
                        'ref|NC_001142|' : 10,
                        'ref|NC_001143|' : 11,
                        'ref|NC_001144|' : 12,
                        'ref|NC_001145|' : 13,
                        'ref|NC_001146|' : 14,
                        'ref|NC_001147|' : 15,
                        'ref|NC_001148|' : 16}

    chromosome_string = row['Chromosome']
    chromosome_number = chromosomes_dic[chromosome_string]
    return(chromosome_number)

def csvDataFrameImport(directory: str, string_length = 7) -> tuple:
    """
    Imports all CSV files from a given directory, filters them based on a condition, and returns a tuple of dataframes and sample names.

    Parameters:
    directory (str): The directory path where the CSV files are located.
    string_length (int, optional): The number of characters to use from the beginning of the file name to identify the sample name. Defaults to 7.

    Returns:
    tuple: A tuple containing a list of dataframes (each dataframe corresponding to a CSV file) and a list of sample names.
    """
    frames_list_samples = []
    sample_names = []

    for counter, filename in enumerate(os.listdir(directory), start=1):
        if filename.endswith(".csv"):
            sample_name = filename[:string_length].strip("-")
            sample_names.append(sample_name)
            path = os.path.join(directory, filename)
            sample = pd.read_csv(path)
            logging.info(f"Processing sample_{counter}: {sample_name}")
            logging.info(f'dataframe {sample_name} has {len(sample)} total rows')
            sample = sample[sample.Type == "SNV"]
            frames_list_samples.append(sample)

    logging.info(f"found {len(frames_list_samples)} samples in {directory}")
    logging.info(f"{len(sample_names)} sample_names: {sample_names}")
    logging.info("Applied SNV only filter to all rows")

    return frames_list_samples, sample_names

def main(input_list) -> None:
    """
    Main function to import CSV files, modify the dataframes, and write the results to an output file.

    Parameters:
    input_list (list): List containing the directory of CSV files and the output file name.

    Returns:
    None
    """
    csv_files_dir = input_list[0]
    output_file_name = input_list[1]
    
    logging.info(f"Searching for csv files in {csv_files_dir}")

    logging.info(f"The String length is set to {STRING_LENGTH}. This is the number of characters that will be used from the beginning of the file name to identify the sample name. You can change this in the script.")
    frames_list_samples, sample_names = csvDataFrameImport(csv_files_dir, string_length=STRING_LENGTH)

    all_frames = []
    for i, frame in enumerate(frames_list_samples):
        logging.debug(frame)
        logging.info(f"{sample_names[i]} has {len(frame)} SNVs")
        frame = modFrame(frame, sampleName=sample_names[i], barcode=sample_names[i])
        all_frames.append(frame)

    result = pd.concat(all_frames)
    result.to_csv(os.path.join(csv_files_dir, output_file_name), sep="\t", mode='w', header=True, index=False)
    logging.info(f"Output file saved as {output_file_name}")

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    description = "This script processes CLC SNV output tables (csv files) so that they can be used in the P-MACD program."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_list', nargs='+', help='Input list containing csv files directory and output file name as command line arguments. Example: python csv_to_paraBed.py /path/to/csv/files/ output_file_name.txt')
    args = parser.parse_args()

    logging.info("Starting csv_to_paraBed.py")
    main(args.input_list)
    logging.info("Finished csv_to_paraBed.py")