# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from .TableSNP import TableSNP
#not sure if this class is usefull for anything...

class TableAllSNP(TableSNP):
    """
    A subclass that represents a table of all SNPs (Single Nucleotide Polymorphisms).

    Inherits from:
    ---------------
    TableSNP

    Methods:
    --------
    minParentalSNP_fraction(self, percent: float) -> None:
        Filters the table to contain only SNPs with a percent value greater than the given threshold.

    """
    def minParentalSNP_fraction(self, percent: int|float) -> None:
        '''Filters the pandas DataFrame of a TableAllSNP object to contain only SNPs 
        with a "Percent" value greater than a given percentage.

        Args:
            percent (int|float): An integer or float representing the minimum percentage 
                of SNPs to keep.

        Returns:
            None
        '''
        self.df = self.df[self.df["Percent"] > percent]