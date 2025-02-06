# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import logging

logger = logging.getLogger('__main__.' + __name__)

#rename to TableSNVs?
class TableSNP:
    """
    A superclass that represents a table of SNPs (Single Nucleotide Polymorphisms).

    Attributes:
    -----------
    df : pandas.DataFrame
        The table of SNPs.
    name : str
        The name of the table.
    usable : bool
        A flag that indicates whether the sample is usable or not.

    Methods:
    --------
    __init__(self, df: pandas.DataFrame, name: str) -> None:
        Initializes a new instance of the TableSNP class.

    __lt__(self, other: TableSNP) -> bool:
        Compares two instances of the TableSNP class based on their names.

    __repr__(self) -> str:
        Returns a string representation of the TableSNP instance.

    setUnusable(self) -> None:
        Sets the usable flag to False.

    printRows(self) -> None:
        Prints the number of rows in the table.

    extractSNVs(self, show: bool = False, log: bool = True) -> None:
        Filters the table to contain only SNVs (Single Nucleotide Variants).

    qualitySNVFilter(self, stringent: bool = False, log: bool = False, diploid: bool = False, 
                     params: dict = {"min_frequency": 70, "min_coverage": 10, 
                                     "max_coverage_stringent": 120, "min_frequency_stringent": 90, 
                                     "min_coverage_stringent": 15}) -> None:
        Filters the table based on quality criteria.

    regionAsInt(self) -> None:
        Converts the "Region" column to integer type.
    """

    def __init__(self, df: pd.DataFrame, name: str) -> None:
        '''Initializes a TableSNP object with a pandas DataFrame and a name.

        Args:
            df (pandas.DataFrame): A pandas DataFrame containing SNP information, 
                including the chromosome, region, reference allele, and alternate allele.
            name (str): A string representing the name of the sample.

        Returns:
            None
        '''
        if isinstance(df, pd.DataFrame):
            self.df = df
            self.name = name
        else:
            raise ValueError("Make sure the input table is a Pandas DataFrame!")
        
        # Usable is a boolean that is set to False if the sample is found to be faulty
        #e.g. bad quality, no parental data, etc. True by default
        self.usable = True

    def __lt__(self, other):
        '''Compares two TableSNP objects based on their name attribute.

        Args:
            other (TableSNP): Another TableSNP object to compare to.

        Returns:
            bool: True if the name of this TableSNP object is less than the name of 
                the other TableSNP object, False otherwise.
        '''
        return self.name < other.name
    
    def __repr__(self) -> str:
        '''Returns a string representation of a TableSNP object.

        Returns:
            str: A string representation of the TableSNP object, including the name 
                attribute, the length of the pandas DataFrame, and the usable attribute.
        '''
        return f"TableSNP({self.name}), df len: {len(self.df)}, usable: {self.usable}"
    
    def setUnusable(self) -> None:
        '''Sets the usable attribute of a TableSNP object to False.

        Returns:
            None
        '''
        self.usable = False

    def printRows(self) -> None:
        '''Prints the number of rows in the pandas DataFrame of a TableSNP object.

        Returns:
            None
        '''
        logging.info(f'Dataframe {self.name} has {len(self.df)} total rows')

    def extractSNVs(self, show: bool = False, log: bool = True) -> None:
        '''Filters the pandas DataFrame of a TableSNP object to contain only SNVs in 
        the "Type" column, and logs the number of SNVs.

        Args:
            show (bool, optional): A boolean indicating whether to print the filtered 
                DataFrame. Defaults to False.
            log (bool, optional): A boolean indicating whether to log the number of 
                SNVs. Defaults to True.

        Returns:
            None
        '''
        self.df = self.df[self.df.Type == "SNV"]
        if show:
            print(self.df)
        logging.info(f'Dataframe {self.name} has {len(self.df)} SNVs')

    def qualitySNVFilter(
        self, stringent: bool = False, 
        diploid: bool = False, 
        params: dict[str, int|float]= {
            "min_frequency": 35,
            "min_coverage": 10,
            "max_coverage_stringent": 400,
            "min_frequency_stringent": 90,
            "min_coverage_stringent": 15,
            "min_QUAL_stringent": 200
        } ) -> None:
        '''Filters SNVs in the pandas DataFrame of a TableSNP object based on various 
        quality parameters, including frequency, coverage, and zygosity.

        Args:
            stringent (bool, optional): A boolean indicating whether to apply stringent 
                filtering parameters. Defaults to False.
            diploid (bool, optional): A boolean indicating whether to filter for diploid 
                SNVs only. Defaults to False.
            params (dict[str, int|float], optional): A dictionary containing the quality 
                filtering parameters. Defaults to {
                    "min_frequency": 70,
                    "min_coverage": 10,
                    "max_coverage_stringent": 400,
                    "min_frequency_stringent": 90,
                    "min_coverage_stringent": 15,
                    "min_QUAL_stringent": 200
                }.

        Returns:
            None
        '''
        initial_len=len(self.df)

        self.df = self.df[self.df.Type =='SNV']

        if not diploid:
            logging.info(f"Sample {self.name} is considered haploid; removing heterozygous calls") 
            self.df = self.df[self.df.Zygosity == "Homozygous"]

        if stringent:
            # Stringent Parameters are max_coverage=120, min_frequency = 90 (%)
            # and are applied on the top of default parameters
            logging.info("Applying stringent SNV filtering")
            self.df = self.df[self.df.Coverage >= params["min_coverage_stringent"]]
            self.df = self.df[self.df.Frequency >= params["min_frequency_stringent"]]
            self.df = self.df[self.df.Coverage <= params["max_coverage_stringent"]]
            self.df = self.df[self.df.QUAL >= params["min_QUAL_stringent"]]

        else:
            self.df = self.df[self.df.Coverage >= params["min_coverage"]]
            self.df = self.df[self.df.Frequency >= params["min_frequency"]]

        self.df = self.df[["Chromosome","Region","Type","Reference","Allele", "Frequency", "Coverage", "Zygosity"]] #added Frequency, Coverage, Zygosity

        final_len=len(self.df)
        logging.info(f'Removed {initial_len-final_len} DataFrame rows through quality filtering. Length of {self.name} is {final_len}')

    def removeHeterozygousCalls(self, min_frequency: int|float = 70) -> None:
        '''
        Removes rows with heterozygous calls and frequency below the specified minimum frequency.

        Args:
            min_frequency (int|float): An integer or float representing the minimum frequency threshold.
                Default value is 70.

        Returns:
            None
        '''
        df_len_initial = len(self.df)
        self.df_SNPs_heterozygous = self.df[(self.df.Zygosity == "Heterozygous") | (self.df.Frequency < min_frequency)].copy()
        self.df = self.df[self.df.Zygosity == "Homozygous"]
        self.df = self.df[self.df.Frequency >= min_frequency]
        df_len_final = len(self.df)
        logging.info(f'Removed {df_len_initial-df_len_final} heterozygous rows. DF {self.name} has now {df_len_final} rows')

    def regionAsInt(self) -> None:
        '''Converts the "Region" column of the pandas DataFrame of a TableSNP object 
        to type int.

        Returns:
            None
        '''
        self.df = self.df.astype({'Region': int})

