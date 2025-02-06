# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import logging
import pandas as pd

logger = logging.getLogger('__main__.' + __name__)

def findSubstituteBase(row: pd.Series) -> str:
    '''Finds a base that can be used to substitute the reference allele in the 
    original reference, allowing for calling of all SNP positions.

    Args:
        row (pandas.Series): A pandas Series object representing a row of a 
            DataFrame containing information about a SNP.

    Returns:
        str: A single character string representing the base that can be used 
            to substitute the reference allele in the original reference.
    '''
    bases=['A','T','G','C']
    allele_p1=row['Allele_x']
    allele_p2=row['Allele_y']
    allele_ref=row['Reference']
    x=0
    while (bases[x] == allele_p1) | (bases[x] == allele_p2) | (bases[x] == allele_ref):
        x+=1
    return(bases[x])
    
def extractSeqFromFastaToRefList(path: str) -> list[list[str]]:
    '''Extracts sequence information from a FASTA file and returns it as a 
    list of lists.

    Args:
        path (str): A string representing the path to the FASTA file.

    Returns:
        list[list[str]]: A list of lists, where each inner list contains two 
            elements: the header line (starting with ">") and the sequence 
            string (without any line breaks).
    '''

    f = open(path,"r")
    contents = f.readlines()
    f.close()

    ref_list=[["",""]]
    chr=''
    seq=''
    for i in contents:
        if '>' in i:
            ref_list.append([i,''])
        else:
            ref_list[-1][1]=ref_list[-1][1]+i.strip()
    logging.info(f"Extraction of Sequence information from {path} finished.")
    return ref_list

def substituteBasesInRefList(ref_list: list[list[str]], new_parent: pd.DataFrame) -> list[list[str]]:
    '''Substitutes bases in a list of reference sequences based on information 
    in a DataFrame.

    Args:
        ref_list (list[list[str]]): A list of lists, where each inner list 
            contains two elements: the header line (starting with ">") and the 
            sequence string (without any line breaks).
        new_parent (pandas.DataFrame): A pandas DataFrame containing 
            information about the substitutions to be made. It must have the 
            following columns: "Chromosome", "Region", and "Substitute".

    Returns:
        list[list[str]]: A list of lists, where each inner list contains two 
            elements: the header line (starting with ">") and the sequence 
            string (without any line breaks), with the substitutions made.
    '''

    for index, row in new_parent.iterrows():
        chromosome = row['Chromosome']
        position = int(row['Region'])
        new_base = row['Substitute']

        for i in range(len(ref_list)):
            
            if chromosome in ref_list[i][0]:
                base_list=list(ref_list[i][1])
                base_list[position-1] = new_base
                new_seq = ''.join(base_list)
                ref_list[i][1]=new_seq
                continue

    return ref_list

def createFastaFromRefList(ref_list: list[list[str]], filename: str) -> None:
    '''Creates a FASTA file from a list of reference sequences.

    Args:
        ref_list (list[list[str]]): A list of lists, where each inner list 
            contains two elements: the header line (starting with ">") and the 
            sequence string (without any line breaks).
        filename (str): A string representing the name of the output file.

    Returns:
        None
    '''

    with open(filename, 'w') as file_out:
        for i in ref_list:
            file_out.write(i[0])
            file_out.write(i[1]+"\n")

def createParentalReferenceGenome(parent1: pd.DataFrame, parent2: pd.DataFrame, original_ref_path: str, new_ref_path: str) -> None:
    '''Creates a new reference genome by combining information from two parental 
    genomes and substituting bases in the original reference genome.

    Args:
        parent1 (pandas.DataFrame): A pandas DataFrame containing information 
            about the SNPs called in the first parental genome. It must have the 
            following columns: "Chromosome", "Region", "Reference", and "Call".
        parent2 (pandas.DataFrame): A pandas DataFrame containing information 
            about the SNPs called in the second parental genome. It must have the 
            following columns: "Chromosome", "Region", "Reference", and "Call".
        original_ref_path (str): A string representing the path to the original 
            reference genome in FASTA format.
        new_ref_path (str): A string representing the path to the new reference 
            genome to be created in FASTA format.

    Returns:
        None
    '''
    #Combine parental calls to create a table of positions to change
    new_parent = pd.merge(parent1, parent2, how="outer", on=['Chromosome','Region','Reference'])

    #find appropriate substitute base to ensure that all SNPs are calld in a new reference
    new_parent['Substitute']=new_parent.apply(lambda row: findSubstituteBase(row), axis=1)
    new_parent = new_parent[['Chromosome','Region','Substitute','Reference']]
    new_parent = new_parent.astype({'Region': int})
    new_parent = new_parent.sort_values(by=["Chromosome","Region"],ascending=(True,True),ignore_index=True)
    logging.info(f"New Parental: {new_parent}")

    #create a list of pairs [[name, seq]] of the original reference
    ref_list = extractSeqFromFastaToRefList(original_ref_path)

    #substitute bases in the original fasta (for yeast takes ~15min) based on
    #parental calls. (needs better implementation for large genomes)
    ref_list = substituteBasesInRefList(ref_list, new_parent)

    #save into a new fasta reference file
    createFastaFromRefList(ref_list, new_ref_path)