# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from .TableSNP import TableSNP
import logging
import pandas as pd

logger = logging.getLogger('__main__.' + __name__)

class SampleTable(TableSNP):
    """
    A class for representing a table of SNP data for a sample.

    Attributes:
        name (str): The name of the sample.
        df (pandas.DataFrame): The SNP data for the sample.
        output_df (pandas.DataFrame): The merged SNP data with parental origin information.
        percent (float): The percentage of known SNPs in the merged data.
        switches_df (pandas.DataFrame): The positions of switches in the merged data.
        output_str (str): A string representation of the output data.

    Methods:
        findSwitches(df, min_SNPs_switch_support=1):
            Finds the positions of switches in the given SNP data.
        innerMergeWithParentalOrigin(parental_origin_db_df, min_SNPs_switch_support=1):
            Merges the given SNP data with parental origin information and finds the positions of switches.
    """
    def __init__(self, df: pd.DataFrame, name: str) -> None:
        """
        Initializes a new instance of the SampleTable class.

        Parameters:
        df (pd.DataFrame): A DataFrame containing SNP data.
        name (str): The name of the sample.

        Returns:
        None
        """
        super().__init__(df, name)
        self.output_df = None
        self.percent = None
        self.switches_df = None
        self.output_str = None
        self.heterozygous = None
        self.aneuploid = None
        self.aneuploid_chromosomes = None
        self.heterozygous_chromosomes = None
    
    def set_heterozygous(self, heterozygous: bool) -> None:
        """
        Sets the heterozygous flag to the given value.

        Parameters:
        heterozygous (bool): The value to set the heterozygous flag to.

        Returns:
        None
        """
        self.heterozygous = heterozygous

    def set_aneuploid(self, aneuploid: bool) -> None:
        """
        Sets the aneuploid flag to the given value.

        Parameters:
        aneuploid (bool): The value to set the aneuploid flag to.

        Returns:
        None
        """
        self.aneuploid = aneuploid

    def findSwitches(self, df: pd.DataFrame, min_SNPs_switch_support: int = 1) -> pd.DataFrame:
        """
        Finds the positions of switches in a DataFrame of SNP data.

        Parameters:
        df (pd.DataFrame): A DataFrame containing SNP data with columns 'Chromosome', 'Parent', and 'Region'.
        min_SNPs_switch_support (int): The minimum number of SNPs required to support a switch. Default is 1.

        Returns:
        pd.DataFrame: A DataFrame containing the positions of switches with columns 'Chromosome', 'Position_1', 'Position_2', and 'Switch_Center'.

        This method iterates through the rows of the input DataFrame and identifies switches between different parents on the same chromosome.
        A switch is defined as a change in parent on the same chromosome.
        The method records the positions of the switches and calculates the switch center as the midpoint between the two positions.
        If the `min_SNPs_switch_support` parameter is greater than 1, the method checks if the switch is supported by enough SNPs before recording it.
        The method returns a DataFrame with the positions of the switches and their switch centers.
        """
        previous_chr="null"
        previous_parent="null"
        previous_pos=0
        output_list=[]

        for i, j in df.iterrows():

            current_chr=j.Chromosome
            current_parent=j.Parent
            current_pos=j.Region

            #check if chromosome has changed
            if previous_chr!=current_chr:
                logging.debug(f"chrmosome {previous_chr} has changed to chromosome {current_chr}")
                previous_parent=current_parent
                previous_chr=current_chr
                previous_pos=0

            #check if parent has changed
            if current_parent != previous_parent:

                if min_SNPs_switch_support > 1:
                    
                    next_pos_check=True
                    #check if the switch is supported by enough SNPs
                    for idx_offset in range(min_SNPs_switch_support-1):

                        if isinstance(i, int):
                            next_row_idx = i + idx_offset + 1

                            #check if the next row is in the dataframe
                            if next_row_idx >= len(df):
                                logging.debug(f"not enough SNPs to support a switch between position {previous_pos} and {current_pos} on chromosome {current_chr} between {previous_parent} and {current_parent}")
                                next_pos_check = False
                                break

                            next_row = df.iloc[i + idx_offset + 1]
                        else:
                            logging.error("i is not an Integer... this should not happen")
                            break

                        next_chr = next_row.Chromosome
                        next_parent = next_row.Parent
                        if (next_parent != current_parent) or (next_chr != current_chr):
                            next_pos_check=False
                            logging.debug(f"switch between position {previous_pos} and {current_pos} on chromosome {current_chr} between {previous_parent} and {current_parent} is not supported by enough SNPs")
                            break
                
                    if next_pos_check==False:
                        continue

                logging.debug(f"found a switch between position {previous_pos} and {current_pos} on chromosome {current_chr} between {current_parent} and {previous_parent}.")
                if previous_pos !=0:
                    output_list.append([current_chr, previous_pos, current_pos])
                    previous_parent=current_parent

            previous_pos=current_pos

        df_out = pd.DataFrame(output_list, columns =['Chromosome', 'Position_1', 'Position_2'])
        df_out["Switch_Center"] = df_out["Position_1"] + ((df_out["Position_2"] - df_out["Position_1"])/2.0)
        return(df_out)

    def innerMergeWithParentalOrigin(
            self, 
            parental_origin_db_df: pd.DataFrame, 
            min_SNPs_switch_support: int = 1
            ):
        """
        Merge the current table with a parental origin database DataFrame, filtering for matching Chromosome, Region, and Allele columns.

        Args:
            parental_origin_db_df (pd.DataFrame): The parental origin database DataFrame to merge with.
            min_SNPs_switch_support (int, optional): The minimum number of SNPs required to support a switch. Defaults to 1.

        Returns:
            pd.DataFrame: The merged DataFrame, sorted by Chromosome and Region.
        """
        import os

        self.output_df = pd.merge(self.df, parental_origin_db_df, how="inner", on=['Chromosome','Region','Allele'], validate="1:1")
        self.output_df = self.output_df.sort_values(by=["Chromosome","Region"],ascending=(True,True),ignore_index=True)
        self.percent=len(self.output_df)/(len(parental_origin_db_df)/2)*100

        self.switches_df = self.findSwitches(self.output_df, min_SNPs_switch_support=min_SNPs_switch_support)

        #check if outputs directory exists, if not, create it
        if not os.path.exists("outputs"):
            os.mkdir("outputs")
        if not os.path.exists("outputs/switches"):
            os.mkdir("outputs/switches")

        self.switches_df.to_csv(f"outputs/switches/{self.name}_switches.csv", index=False)

        logging.info(f'length of filtered {self.name} is: {len(self.df)}, parents merged: {len(self.output_df)}, ({round(self.percent,3)}% of known SNPs), {len(self.switches_df)} switches')
        self.output_str = f'{self.name},{len(self.df)},{len(self.output_df)},{round(self.percent,3)},{len(self.switches_df)}\n'

        return(self.output_df)
    
    def findHeterozygousChromosomes(self) -> list:
        """
        Finds the chromosomes that are heterozygous by checking the frequency of heterozygous SNPs in each chromosome.

        Returns:
        -------
        list
            A list of the chromosomes that are heterozygous.

        Examples:
        --------
        >>> sample_table = SampleTable()
        >>> sample_table.findHeterozygousChromosomes()
        ['chr1', 'chr3', 'chr5']
        """
        self.df["zygosity_from_freq"] = self.df["Frequency"].apply(lambda x: "Homozygous" if x >= 70 else ("Heterozygous" if ((x >= 35) and (x < 70)) else "Unknown"))
        chromosome_heterozygosity_dict = {}
        for chromosome in self.df["Chromosome"].unique():
            chromosome_df = self.df[self.df["Chromosome"] == chromosome]
            #first check if "Heterozygous" is in the value counts, if not, then the chromosome is not heterozygous
            logging.debug(f"{self.name} chromosome {chromosome} value counts: {chromosome_df['zygosity_from_freq'].value_counts()}")
            if "Heterozygous" not in chromosome_df["zygosity_from_freq"].value_counts().index:
                continue
            
            chromosome_heterozygosity_dict[chromosome] = chromosome_df["zygosity_from_freq"].value_counts(normalize=True)["Heterozygous"]
       
        logging.debug(f"{self.name} chromosome heterozygosity dictionary: {chromosome_heterozygosity_dict}")
        heterozygous_chromosomes = [chromosome for chromosome in chromosome_heterozygosity_dict if chromosome_heterozygosity_dict[chromosome] > 0.2]
        logging.info(f"{self.name} heterozygous chromosomes: {heterozygous_chromosomes}")
        
        #if heterozygous_chromosomes is not empty, set self.heterozygous_chromosomes to heterozygous_chromosomes
        if len(heterozygous_chromosomes) > 0:
            self.set_heterozygous(True)
            self.heterozygous_chromosomes = heterozygous_chromosomes
        else:
            self.set_heterozygous(False)
            self.heterozygous_chromosomes = []

        return heterozygous_chromosomes

    def findAnuploidChromosomesBySNPCoverage(
            self,
            show_details: bool = False,
            plot_coverage: bool = False,
            save_plot: bool = False,
            save_high_coverage_SNP_positions: bool = False,
            high_coverage_threshold: int|float = 2,
            low_coverage_threshold: float = 0.5,
            remove_aneuploid_from_high_coverage_SNP_positions: bool = True
            ) -> list:
        """
        Finds the chromosomes that are aneuploid by comparing their SNP coverage to the mean coverage.

        Args:
        ----------
        show_details : bool, optional
            If True, prints the chromosome coverage dataframe (default is False).
        plot_coverage : bool, optional
            If True, plots the genomic SNP coverage (default is False).
        save_plot : bool, optional
            If True, saves the plot to a file (default is False).
        save_high_coverage_SNP_positions : bool, optional
            If True, saves the high/low coverage SNP positions to a file (default is False).
        high_coverage_threshold : int|float, optional
            The threshold for high coverage (default is 2).
        low_coverage_threshold : float, optional
            The threshold for low coverage (default is 0.5).
        remove_aneuploid_from_high_coverage_SNP_positions : bool, optional
            If True, removes the aneuploid chromosomes from the high/low coverage SNP positions (default is True).

        Returns:
        -------
        list
            A list of the chromosomes that are aneuploid.

        Examples:
        --------
        >>> sample_table = TableSNP("sample", pd.DataFrame())
        >>> sample_table.findAnuploidChromosomesBySNPCoverage(show_details=True, plot_coverage=True, save_plot=True)
        ['chr2', 'chr4', 'chr6']
        """
        #find the average coverage for each chromosome ("Coverage" column)
        #subset only chromosome and coverage columns
        #find the overall mean coverage

        #if self.df is empty, then return an empty list
        import os
        from mayo import plot_genomic_SNP_coverage

        if len(self.df) == 0:
            logging.warning(f"{self.name} is empty, cannot find aneuploid chromosomes by SNP coverage")
            self.aneuploid_chromosomes = []
            self.heterozygous_chromosomes = []
            return []

        try:
            from settings import config as cfg
            mito_chromosomes = cfg["mitochondrial_chr_name"]
        except:
            mito_chromosomes = ["chrM", "chrMT", "M", "MT", "ref|NC_001224|"]
            logging.warning("mitochondrial_chr_name not found in config file, using default mitochondrial chromosome names")

        #find the average coverage for each chromosome ("Coverage" column, excluding mitochondrial chromosomes)
        df_SNV_coverage = self.df[["Chromosome", "Coverage"]].copy()
        df_SNV_coverage = df_SNV_coverage[~df_SNV_coverage["Chromosome"].isin(mito_chromosomes)]
        mean_coverage = df_SNV_coverage["Coverage"].mean()

        #group by chromosome and find the mean coverage
        df_chromosome_coverage = df_SNV_coverage.groupby("Chromosome").mean()

        #normalize the coverage by dividing by the mean coverage
        df_chromosome_coverage["Coverage_norm"] = df_chromosome_coverage["Coverage"].apply(lambda x: x/mean_coverage)
        
        #find the chromosomes that have coverage that is .5x higher or lower than the mean coverage
        aneuploid_chromosomes = df_chromosome_coverage[(df_chromosome_coverage["Coverage_norm"] >= 1.5) | (df_chromosome_coverage["Coverage_norm"] <= 0.5)].index.tolist()
        mean_coverage_recalc = mean_coverage
        if len(aneuploid_chromosomes) > 0:
            aneuploid_chromosomes_len = len(aneuploid_chromosomes)
            reculculation_count = 1
            while True:
                #recalculate the mean coverage, excluding the aneuploid chromosomes
                df_chromosome_coverage_recalc = df_SNV_coverage[~df_SNV_coverage["Chromosome"].isin(aneuploid_chromosomes)]
                mean_coverage_recalc = df_chromosome_coverage_recalc["Coverage"].mean()
                logging.debug(f"{self.name} mean coverage: {mean_coverage}, recalculated ({reculculation_count}) mean coverage: {mean_coverage_recalc}")

                #normalize the coverage again by dividing by the mean coverage
                df_chromosome_coverage["Coverage_norm"] = df_chromosome_coverage["Coverage"].apply(lambda x: x/mean_coverage_recalc)

                #again, find the chromosomes that have coverage that is .5x higher or lower than the mean coverage
                # this time, the aneuploid chromosomes are calculated using the recalculated mean coverage, just in case the aneuploid chromosomes are affecting the mean coverage and causing false negatives
                aneuploid_chromosomes = df_chromosome_coverage[(df_chromosome_coverage["Coverage_norm"] > 1.5) | (df_chromosome_coverage["Coverage_norm"] < 0.5)].index.tolist()

                if len(aneuploid_chromosomes) == aneuploid_chromosomes_len:
                    break

                #if the reculculation_count is greater than the number of chromosomes, then break
                elif reculculation_count > len(df_chromosome_coverage):
                    logging.warning(f"{self.name} reculculation_count ({reculculation_count}) is greater than the number of chromosomes ({len(df_chromosome_coverage)}), breaking")
                    break

                else:
                    aneuploid_chromosomes_len = len(aneuploid_chromosomes)
                    reculculation_count += 1

        #if aneuploid_chromosomes is not empty, set self.aneuploid_chromosomes to aneuploid_chromosomes
        if len(aneuploid_chromosomes) > 0:
            self.set_aneuploid(True)
            self.aneuploid_chromosomes = aneuploid_chromosomes
        else:
            self.set_aneuploid(False)
            self.aneuploid_chromosomes = []

        logging.info(f"{self.name} chr with coverage 2x different than mean (all:{mean_coverage}, anu-excluded:{mean_coverage_recalc}): {aneuploid_chromosomes}")
        if show_details:
            logging.info(f"{self.name} chromosome coverage dataframe: {df_chromosome_coverage}")

        if plot_coverage or save_high_coverage_SNP_positions:
            df_coverage = self.df[["Chromosome", "Region", "Coverage","zygosity_from_freq"]].copy()
            df_coverage = df_coverage[~df_coverage["Chromosome"].isin(mito_chromosomes)]
            df_coverage["Coverage_norm"] = df_coverage["Coverage"].apply(lambda x: x/mean_coverage_recalc)
            
            if plot_coverage:
                plot_genomic_SNP_coverage(df_coverage, f"{self.name}_(APSC-{round(mean_coverage_recalc, 2)})", save=save_plot, show=False) #APSC = average per SNP coverage

            if save_high_coverage_SNP_positions:
                df_coverage_high = df_coverage[df_coverage["Coverage_norm"] > high_coverage_threshold]
                df_coverage_low = df_coverage[df_coverage["Coverage_norm"] < low_coverage_threshold]

                #remove chromosomes that are aneuploid from the high coverage dataframe, if flag is set
                if remove_aneuploid_from_high_coverage_SNP_positions:
                    df_coverage_high = df_coverage_high[~df_coverage_high["Chromosome"].isin(self.aneuploid_chromosomes)]
                    df_coverage_low = df_coverage_low[~df_coverage_low["Chromosome"].isin(self.aneuploid_chromosomes)]

                df_coverage_high["sample_name"] = self.name
                #check if the directory exists, if not, create it
                if not os.path.exists("outputs/high_coverage_SNP_positions"):
                    os.makedirs("outputs/high_coverage_SNP_positions")
                df_coverage_high.to_csv(f"outputs/high_coverage_SNP_positions/{self.name}_high_coverage_SNP_positions.csv", index=False)

                if not os.path.exists("outputs/low_coverage_SNP_positions"):
                    os.makedirs("outputs/low_coverage_SNP_positions")
                df_coverage_low.to_csv(f"outputs/low_coverage_SNP_positions/{self.name}_low_coverage_SNP_positions.csv", index=False)
            

        #return a list of the chromosomes that are aneuploid
        return aneuploid_chromosomes
    
