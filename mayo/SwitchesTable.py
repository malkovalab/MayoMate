# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from .TableSNP import TableSNP
import random
import pandas as pd
import logging

logger = logging.getLogger('__main__.' + __name__)

class SwitchesTable(TableSNP):
    """
    A class representing a table of switches (recombinations) in a sample.

    Attributes:
    -----------
    df : pandas.DataFrame
        A DataFrame containing the switches data.
    name : str
        The name of the sample.
    percent : float
        The percentage of switches in the sample.

    Methods:
    --------
    __init__(self, df, name, percent):
        Initializes a SwitchesTable object with the given DataFrame, sample name, and percentage of switches.
    __repr__(self):
        Returns a string representation of the SwitchesTable object.
    simulateRandomSwitches(self, chromosomes):
        Simulates random switches for the sample based on the given chromosomes.
    """

    def __init__(self,
        df: pd.DataFrame,
        name: str, 
        percent: int|float,
        aneuploid: bool | None = None,
        heterozygous: bool | None = None,
        aneuploid_chromosomes: list | None = None,
        heterozygous_chromosomes: list | None = None,
        usable: bool | None = None,
        ) -> None:
        '''Initializes a SwitchesTable object with a pandas DataFrame, a name, and a 
        percentage value.

        Args:
            df (pandas.DataFrame): A pandas DataFrame containing SNP information, 
                including the chromosome, region, reference allele, and alternate allele.
            name (str): A string representing the name of the sample.
            percent (int|float): An integer or float representing the percentage of 
                SNPs to keep.

        Returns:
            None
        '''
        TableSNP.__init__(self, df, name)
        self.percent = percent
        self.aneuploid = aneuploid
        self.heterozygous = heterozygous
        self.aneuploid_chromosomes = aneuploid_chromosomes
        self.heterozygous_chromosomes = heterozygous_chromosomes
        self.usable = usable

    def __repr__(self) -> str:
        '''Returns a string representation of a SwitchesTable object.

        Returns:
            str: A string representation of the SwitchesTable object, including the name 
                attribute, the length of the pandas DataFrame, the percentage of parental 
                SNPs to keep, and the usable attribute.
        '''
        return f"SwitchesTable({self.name}), df len: {len(self.df)}, {self.percent}% of Parental SNPs, usable: {self.usable}"

    def simulateRandomSwitches(self, chromosomes: list[list]) -> None:
        """
        Simulates random switches for each chromosome in the sample.

        Parameters:
        chromosomes (list[list]): A list of lists containing chromosome names and lengths.

        Returns:
        None

        This method generates a new DataFrame with random switch positions for each chromosome in the sample.
        For each chromosome, the method finds a random position and records it as a switch center.
        The number of switches per chromosome is determined by the number of real switches in the sample.
        If a chromosome in the sample is not found in the provided list of chromosomes, an error message is logged
        and the chromosome is omitted from the simulation.
        """
        chromosome_dict = { k[0]: k[1] for k in chromosomes }
        columns=['Chromosome', 'Switch_Center']

        random_switches_list=[]

        #iterate through the sample's real switches chromosomes
        switches_series = self.df['Chromosome'].value_counts()
        switches_dict = dict(switches_series.to_dict()) #convert to dict, in old version switches_dict = dict(switches_series)
        for chr in switches_dict:
            #record the number of switches per chromosome
            chromosome_switches = switches_dict[chr]
            chromosome_switches = int(chromosome_switches)

            if chr not in chromosome_dict:
                #if it's a mitochondrial chromosome, skip the warning, but don't add it to the random switches either
                if (chr == 'chrM') | (chr == 'chrMT') | (chr == 'ref|NC_001224|'):
                    logging.debug(f"Chromosome {chr} is mitochondrial and not found in sample {self.name}'. Omitting...")
                else:
                    logging.warning(f"Chromosome {chr} not found in sample '{self.name}'. Omitting...")
                continue

            chromosome_length = chromosome_dict[chr]

            #for each real switch, find a random position for that chromosome
            for switch in range(chromosome_switches):
                random_switch_position = random.randint(0, chromosome_length)
                random_switches_list.append([chr, random_switch_position])

        #create a new dataframe with random switch positions
        self.random_switch_df = pd.DataFrame(random_switches_list, columns = columns)

    def readjustSwitchesToSPO11Signal(
            self, 
            spo11_oligo_map_df: pd.DataFrame, 
            adjust_window = 2000,
            mode = 'nearest_peak'
            
            ) -> pd.DataFrame:
        """
        Adjusts the center of each switch in the switches dataframe to the nearest SPO11-oligo peak within a specified window.

        Parameters:
        spo11_oligo_map_df (pd.DataFrame): A DataFrame containing the SPO11-oligo map data. 
                                        It should have a "Chromosome" column and a "Position" column.
        adjust_window (int, optional): The window within which to search for the nearest SPO11-oligo peak. 
                                    Defaults to 2000.

        Returns:
        pd.DataFrame: The switches dataframe with an additional "Adjusted_Center" column. 
                    The "Adjusted_Center" column contains the position of the nearest SPO11-oligo peak within the window for each switch.
                    If no peak is found within the window for a switch, the "Adjusted_Center" is set to the original switch center.

        """
        #by default, the adjusted center is the same as the switch center
        self.df["Adjusted_Center"] = self.df["Switch_Center"]

        for chromosome in self.df["Chromosome"].unique():
            print(f"adjusting switches on chromosome {chromosome} for sample {self.name}")

            chromosome_switches_df = self.df[self.df["Chromosome"] == chromosome]
            chromosome_spo11_oligo_df = spo11_oligo_map_df[spo11_oligo_map_df["Chromosome"] == chromosome]
            print(f"There are {len(chromosome_spo11_oligo_df)} SPO11-oligo peaks on chromosome {chromosome}")

            for i, switch in chromosome_switches_df.iterrows():
                switch_center = switch["Switch_Center"]
                switch_pos1 = switch["Position_1"]
                switch_pos2 = switch["Position_2"]
                
                #find the SPO11-oligo peaks within the switch positions
                switches_peaks_df = chromosome_spo11_oligo_df[(chromosome_spo11_oligo_df["Position"] > switch_pos1) & (chromosome_spo11_oligo_df["Position"] < switch_pos2)].copy()
                if len(switches_peaks_df) == 0:
                    logging.debug(f"no SPO11-oligo peak found within the switch positions for switch center {switch_center} on chromosome {chromosome}")
                    continue

                # find the SPO11-oligo peaks within the window
                window_peak_df = switches_peaks_df[(switches_peaks_df["Position"] >= switch_center - adjust_window) & (switches_peaks_df["Position"] <= switch_center + adjust_window)].copy()

                if len(window_peak_df) == 0:
                    logging.debug(f"no SPO11-oligo peak found within the {adjust_window}bp window for switch center {switch_center} on chromosome {chromosome}")
                    continue

                # find row with highest value (strongest peak) within the window
                df_max = window_peak_df[window_peak_df["Value"] == window_peak_df["Value"].max()].copy()

                #check if there are multiple peaks with the same value
                if len(df_max) > 1:
                    logging.debug(f"multiple SPO11-oligo peaks with the same value within the {adjust_window}bp window for switch center {switch_center} on chromosome {chromosome}")
                    print(f"multiple SPO11-oligo peaks with the same value within the {adjust_window}bp window for switch center {switch_center} on chromosome {chromosome}")
                    #if there are multiple peaks with the same value, choose the closest one to the original switch center
                    df_max["Distance"] = abs(df_max["Position"] - switch_center)
                    df_max = df_max[df_max["Distance"] == df_max["Distance"].min()].copy()

                strongest_peak_pos = df_max["Position"].values[0]
                strongest_peak_val = df_max["Value"].values[0]

                #find the nearest peak within the window and the switch positions
                nearest_peak_pos = window_peak_df.iloc[(window_peak_df["Position"]-switch_center).abs().argsort()[:1]]["Position"].values[0]
                nearest_peak_val = window_peak_df.iloc[(window_peak_df["Position"]-switch_center).abs().argsort()[:1]]["Value"].values[0]


                if mode == 'nearest_peak':
                    output_peak = nearest_peak_pos
                    output_val = nearest_peak_val

                elif mode == 'strongest_peak':
                    output_peak = strongest_peak_pos
                    output_val = strongest_peak_val

                else:
                    logging.error(f"Invalid mode '{mode}' for adjusting switches to SPO11-oligo peaks. Using default mode 'nearest_peak'.")
                    print(f"Invalid mode '{mode}' for adjusting switches to SPO11-oligo peaks. Using default mode 'nearest_peak'.")
                    output_peak = nearest_peak_pos
                    output_val = nearest_peak_val

                #set the adjusted center to the nearest peak position if the peak has more than the minimum number of reads
                self.df.loc[i, "Adjusted_Center"] = output_peak
                logging.info(f"adjusted switch center pos. {switch_center} (between SNPs {switch_pos1} and {switch_pos2}) by {abs(switch_center - output_peak)}bp to {output_peak} value {output_val} on chromosome {chromosome}")
                print(f"adjusted switch center pos. {switch_center} (between SNPs {switch_pos1} and {switch_pos2}) by {abs(switch_center - output_peak)}bp to {output_peak} value {output_val} on chromosome {chromosome}")
