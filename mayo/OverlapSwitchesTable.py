# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from .TableSNP import TableSNP
from .SwitchesTable import SwitchesTable
from typing import Union, Literal, overload, TypeVar, List
from numpy import NaN
import pandas as pd
import logging

logger = logging.getLogger('__main__.' + __name__)

class OverlapSwitchesTable(TableSNP):
    # Instances of this class are used for Overlapping SwitchesTable (switches) with
    # with a SNP tables and cluster location tables

    #inherit init from TableSNP but add self.clusters_dict_JT = {} and self.clusters_dict_PMACD = {}
    def __init__(self, name: str, df: pd.DataFrame) -> None:
        """
        Initializes an OverlapSwitchesTable object.

        Args:
            name (str): The name of the OverlapSwitchesTable object.
            df (pd.DataFrame): A pandas DataFrame containing SNP data.
        """
        super().__init__(name, df)
        self.cluster_dict_JT = {}
        self.cluster_dict_PMACD = {}
        self.totalssDNA_JT_dict = {}
        self.totalssDNA_PMACD_dict = {}

        self.GC_dict = {}
        self.CO_dict = {}

    def __repr__(self) -> str:
        """
        Returns a string representation of the OverlapSwitchesTable object.

        Returns:
            str: A string representation of the OverlapSwitchesTable object.
        """
        #check what self properties are available
        print_str=f"OverlapSwitchesTable({self.name}), df len: {len(self.df)}"
        if hasattr(self, 'percent'):
            print_str+=f", {self.percent}% of Parental SNPs"
        if hasattr(self, 'df_SNPs'):
            print_str+=f", df_SNPs len: {len(self.df_SNPs)}"
        if hasattr(self, 'switchesTable'):
            print_str+=f", switchesTable len: {len(self.switchesTable)}"
        if hasattr(self, 'randomSwitchesTable'):
            print_str+=f", randomSwitchesTable len: {len(self.randomSwitchesTable)}"


        print_str+=f", usable: {self.usable}"

        return print_str
    
    def addSwitchesTable(self, SwitchesTable: SwitchesTable, low_percent_warning: int|float = 70) -> None:
        """
        Adds a SwitchesTable object to the OverlapSwitchesTable object.

        Args:
            SwitchesTable (SwitchesTable): A SwitchesTable object to add to the OverlapSwitchesTable object.
            low_percent_warning (int|float, optional): A percentage threshold for a warning message. Defaults to 70.

        Raises:
            ValueError: If the name of the SwitchesTable object does not match the name of the OverlapSwitchesTable object.
        """
        if self.name == SwitchesTable.name:
            self.switchesTable = SwitchesTable.df
            self.percent = SwitchesTable.percent
            self.aneuploid = SwitchesTable.aneuploid
            self.heterozygous = SwitchesTable.heterozygous
            self.aneuploid_chromosomes = SwitchesTable.aneuploid_chromosomes
            self.heterozygous_chromosomes = SwitchesTable.heterozygous_chromosomes
            self.usable = SwitchesTable.usable
            
            logging.info(f'The length of {self.name} is: {len(self.df)}, {len(self.switchesTable)} switches')
            
            if self.percent < low_percent_warning:
                logging.warning(f'Warning: {self.name} has a low percentage of detected parental SNPs ({self.percent}). Switches may not be accurate. Setting {self.name} as unusable.')
                self.setUnusable()

        else:
            logging.error(f'Error: cannot add switchesTable property. {self.name} and {SwitchesTable.name} have different name properties.')
            raise ValueError(f'Error: cannot add switchesTable property. {self.name} and {SwitchesTable.name} have different name properties.')

    def addRandomSwitchesTable(self, SwitchesTable: SwitchesTable) -> None:
        """
        Adds a random switches table to the OverlapSwitchesTable object.

        Args:
            SwitchesTable (SwitchesTable): A SwitchesTable object containing random switches to add to the OverlapSwitchesTable object.
        """
        if self.name == SwitchesTable.name:
            self.randomSwitchesTable = SwitchesTable.random_switch_df
            self.percent = SwitchesTable.percent
            
            logging.debug(f'Added {len(self.randomSwitchesTable)} random switches to {self.name}')
            
            if self.percent < 0.85:
                logging.warning(f'Warning: {self.name} has a low percentage of detected parental SNPs ({self.percent}). Switches may not be accurate. Setting {self.name} as unusable.')
                self.setUnusable()
        else:
            logging.error(f'Error: cannot add switchesTable property. {self.name} and {SwitchesTable.name} have different name properties.')
            raise ValueError(f'Error: cannot add switchesTable property. {self.name} and {SwitchesTable.name} have different name properties.')

    def _findClusteredSNPs_PMACD(self, pval_key="p001_def"):
        
        self.df_SNPs = self.df.copy()

        chr_dict = {'chr1':'ref|NC_001133|',
                    'chr2':'ref|NC_001134|',
                    'chr3':'ref|NC_001135|',
                    'chr4':'ref|NC_001136|',
                    'chr5':'ref|NC_001137|',
                    'chr6':'ref|NC_001138|',
                    'chr7':'ref|NC_001139|',
                    'chr8':'ref|NC_001140|',
                    'chr9':'ref|NC_001141|',
                    'chr10':'ref|NC_001142|',
                    'chr11':'ref|NC_001143|',
                    'chr12':'ref|NC_001144|',
                    'chr13':'ref|NC_001145|',
                    'chr14':'ref|NC_001146|',
                    'chr15':'ref|NC_001147|',
                    'chr16':'ref|NC_001148|'}
        
        df_clusters_PMACD = self.cluster_dict_PMACD[pval_key]
        df_clusters_PMACD = df_clusters_PMACD.replace({"Chromosome":chr_dict})

        self.df_SNPs["Is_Clustered"] = False
        self.df_SNPs["Cluster_ID"] = NaN
        self.df_SNPs["Cluster_Length"] = NaN
        self.df_SNPs["Cluster_Start"] = NaN
        self.df_SNPs["Cluster_End"] = NaN

        for snp_index, snp_row in self.df_SNPs.iterrows():
            position_snp=snp_row["Region"]
            chromosome_snp=snp_row["Chromosome"]

            for index, row in df_clusters_PMACD.iterrows():
                chromosome_cluster = row["Chromosome"]
                #cluster_start = row["cluster_start"]
                cluster_start = row["Start_position"]
                #cluster_end = row["cluster_end"]
                cluster_end = row["End_position"]
                cluster_id = row["Dataset_Cluster_ID"]
                cluster_length = row["Cluster_Length"]

                if chromosome_snp == chromosome_cluster:
                    if (position_snp >= cluster_start) & (position_snp <= cluster_end):
                        self.df_SNPs.at[snp_index, "Is_Clustered"] = True
                        self.df_SNPs.at[snp_index, "Cluster_ID"] = cluster_id
                        self.df_SNPs.at[snp_index, "Cluster_Length"] = cluster_length
                        self.df_SNPs.at[snp_index, "Cluster_Start"] = cluster_start
                        self.df_SNPs.at[snp_index, "Cluster_End"] = cluster_end
        

        list_of_cluster_ids = self.df_SNPs["Cluster_ID"].unique()
        list_of_cluster_ids = [x for x in list_of_cluster_ids if str(x) != 'nan'] #remove NaN from list
        logging.info(f'{self.name} has {len(list_of_cluster_ids)} PMACD clusters; {len(self.df_SNPs[self.df_SNPs["Is_Clustered"] == True])} clustered SNVs')

    def _findClusteredSNPs_JT(self, pval_key=0.001):

        self.df_SNPs = self.df.copy()

        df_clusters_JT = self.cluster_dict_JT[pval_key] 
        df_clusters_JT["Sample_Cluster_ID"] = df_clusters_JT.index + 1

        self.df_SNPs["Is_Clustered"] = False
        self.df_SNPs["Cluster_ID"] = NaN
        self.df_SNPs["Cluster_Length"] = NaN
        self.df_SNPs["Cluster_Start"] = NaN
        self.df_SNPs["Cluster_End"] = NaN

        for snp_index, snp_row in self.df_SNPs.iterrows():
            position_snp=snp_row["Region"]
            chromosome_snp=snp_row["Chromosome"]

            for index, row in df_clusters_JT.iterrows():
                chromosome_cluster = row["Chromosome"]
                cluster_start = row["Start"]
                cluster_end = row["End"]
                cluster_id = row["Sample_Cluster_ID"] #"Dataset_Cluster_ID"
                cluster_length = row["Length"]

                if chromosome_snp == chromosome_cluster:
                    if (position_snp >= cluster_start) & (position_snp <= cluster_end):
                        self.df_SNPs.at[snp_index, "Is_Clustered"] = True
                        self.df_SNPs.at[snp_index, "Cluster_ID"] = cluster_id
                        self.df_SNPs.at[snp_index, "Cluster_Length"] = cluster_length
                        self.df_SNPs.at[snp_index, "Cluster_Start"] = cluster_start
                        self.df_SNPs.at[snp_index, "Cluster_End"] = cluster_end


        list_of_cluster_ids = self.df_SNPs["Cluster_ID"].unique()
        list_of_cluster_ids = [x for x in list_of_cluster_ids if str(x) != 'nan'] #remove NaN from list
        logging.info(f'{self.name} has {len(list_of_cluster_ids)} JT clusters; {len(self.df_SNPs[self.df_SNPs["Is_Clustered"] == True])} clustered SNVs')

    def findClusteredSNPs(
            self,
            df_clusters_type: Literal["PMACD", "JT"] = "JT",
            cluster_significance: Union[float, str] = 0.001
            ) -> None:
        """
        Finds clustered SNPs in the OverlapSwitchesTable object.

        Args:
            df_clusters_PMACD (pd.DataFrame): A pandas DataFrame containing information about SNP clusters.

        Returns:
            None
        """
        if df_clusters_type == "PMACD":
            self._findClusteredSNPs_PMACD(pval_key=cluster_significance)

        if df_clusters_type == "JT":
            self._findClusteredSNPs_JT(pval_key=cluster_significance)

    def findClusters(self, window_size: int) -> None:

        #get all clustered SNPs (i.e. SNPs that are in a cluster).
        self.df_SNPs = self.df_SNPs[self.df_SNPs["Is_Clustered"] == True]

        #get the Chromosome and Switch_Center columns from the switchesTable, convert to a dictionary
        #where the keys are the chromosome names and the value is a list containing all switch centers on that chromosome
        switchesDict = self.switchesTable[["Chromosome","Switch_Center"]].groupby("Chromosome").agg(list).to_dict()["Switch_Center"]
        
        clusters_num = self.df_SNPs["Cluster_ID"].unique().shape[0]
        clustered_snvs_num = len(self.df_SNPs)
        chromosomes_with_switches = len(switchesDict)
        logging.debug(f'{self.name}: clusters={clusters_num}, clustered SNVs={clustered_snvs_num}, chromosomes with switches={chromosomes_with_switches}... check debug for details')
        logging.debug(f'switches_dict: {switchesDict}')
        logging.debug(f'Finding clusters within {window_size} bp of switches...')

        #group df_SNPs by Cluster_ID, aggregate the min and max of the Region column, and the Chromosome column
        #this will give us the min and max positions of each cluster, and the chromosome it's on
        df_SNPs_clusters = self.df_SNPs.copy().groupby("Cluster_ID").agg({"Region":["min","max"],"Chromosome":"first"}).reset_index()
        df_SNPs_clusters.columns = ["Cluster_ID","Cluster_Start","Cluster_End","Chromosome"]

        #cluster_start and cluster_end are the min and max positions of the cluster and should be integers
        df_SNPs_clusters["Cluster_Start"] = df_SNPs_clusters["Cluster_Start"].astype(int)
        df_SNPs_clusters["Cluster_End"] = df_SNPs_clusters["Cluster_End"].astype(int)

        df_SNPs_clusters["Switch_Proximal"] = False
        self.df_SNPs["Switch_Proximal"] = False
        self.df_SNPs["Clustered_Proximal"] = False
        self.df_SNPs["Window_Size"] = window_size

        logging.debug(df_SNPs_clusters.head())

        #if the chromosome of the cluster is in the switchesDict, then we want to check if the switch center is within the window size of the cluster's min and max positions
        #if it is, then we want to set the Switch_Proximal column to True for all SNPs in that cluster
        #first, we need to create a new column in df_SNPs_clusters that will hold the switch centers for each cluster
        #we do this by using the switchesDict to get the switch center for the chromosome of the cluster
        #then we use the apply function to apply the lambda function to each row of the df_SNPs_clusters dataframe
        #the lambda function takes the switch center and the cluster start and end positions, and checks if the switch center is within the window size of the cluster start and end positions
        #if it is, then it returns True, otherwise it returns False
        #do it only if df_SNPs_clusters has rows in it

        if df_SNPs_clusters.shape[0] == 0:
            #it is possible that the chromosome of the cluster is not in the switchesDict
            #if it is not, then we need to handle that
            logging.warning(f"No clustered SNPs found in {self.name}")
        
        else:
            #df_SNPs_clusters["Candidate_Switches"] = df_SNPs_clusters.apply(lambda x: switchesDict[x["Chromosome"]] if x["Chromosome"] in switchesDict.keys() else [], axis=1)
            df_SNPs_clusters["Candidate_Switches"] = df_SNPs_clusters["Chromosome"].apply(lambda x: switchesDict.get(x, []))
            logging.debug(df_SNPs_clusters.head())

            #if the cluster is within the window size of a switch (provided in col ["Candidate_Switches"]), then we want to set the Switch_Proximal column to True for all SNPs in that cluster
            df_SNPs_clusters["Switch_Proximal"] = df_SNPs_clusters.apply(lambda x: any(((x["Cluster_Start"]  <= switch + window_size) & (switch - window_size <= x["Cluster_End"])) for switch in x["Candidate_Switches"]), axis=1)
            logging.debug(df_SNPs_clusters.head())

            #now we want to set the Switch_Proximal column in the df_SNPs dataframe to True for all SNPs in clusters that are proximal to a switch
            #we do this by using the isin function to get all the Cluster_IDs that are proximal to a switch
            #then we use the isin function to set the Switch_Proximal column to True for all SNPs in those clusters
            self.df_SNPs.loc[self.df_SNPs["Cluster_ID"].isin(df_SNPs_clusters[df_SNPs_clusters["Switch_Proximal"] == True]["Cluster_ID"]), "Switch_Proximal"] = True

            #now we want to set the Clustered_Proximal column to True for all SNPs in clusters that are proximal to a switch
            self.df_SNPs["Clustered_Proximal"] = (self.df_SNPs["Switch_Proximal"] == True) & (self.df_SNPs["Is_Clustered"] == True)

    def find_GC_Clusters(self, window_size: int) -> None:

        #get all clustered SNPs (i.e. SNPs that are in a cluster). That's what we're interested in
        self.df_SNPs = self.df_SNPs[self.df_SNPs["Is_Clustered"] == True]
        logging.info(f'{self.name} has {len(self.df_SNPs)} clustered SNVs')
        logging.info(f'{self.name} has {len(self.GC_dict)} chromosomes with GCs.')

        #group df_SNPs by Cluster_ID, aggregate the min and max of the Region column, and the Chromosome column
        #this will give us the min and max positions of each cluster, and the chromosome it's on

        df_SNPs_clusters = self.df_SNPs.copy()

        # df_SNPs_clusters["GC_Proximal"] = False
        # df_SNPs_clusters["Clustered_GC_Proximal"] = False
        # df_SNPs_clusters["Window_Size_GC"] = window_size
        
        df_SNPs_clusters = df_SNPs_clusters.groupby("Cluster_ID").agg({"Region":["min","max"],"Chromosome":"first"}).reset_index()
        df_SNPs_clusters.columns = ["Cluster_ID","Cluster_Start","Cluster_End","Chromosome"]

        #cluster_start and cluster_end are the min and max positions of the cluster and should be integers
        df_SNPs_clusters["Cluster_Start"] = df_SNPs_clusters["Cluster_Start"].astype(int)
        df_SNPs_clusters["Cluster_End"] = df_SNPs_clusters["Cluster_End"].astype(int)
        df_SNPs_clusters["GC_Proximal"] = False


        logging.debug(df_SNPs_clusters.head())

        if df_SNPs_clusters.shape[0] == 0:
            #it is possible that the chromosome of the cluster is not in the switchesDict
            #if it is not, then we need to handle that
            logging.warning(f"No clustered SNPs found in {self.name}")
            self.df_SNPs["GC_Proximal"] = False
            
        else:
            #df_SNPs_clusters["Candidate_Switches"] = df_SNPs_clusters.apply(lambda x: switchesDict[x["Chromosome"]] if x["Chromosome"] in switchesDict.keys() else [], axis=1)
            df_SNPs_clusters["Candidate_GC"] = df_SNPs_clusters["Chromosome"].apply(lambda x: self.GC_dict.get(x, []))
            logging.debug(df_SNPs_clusters.head())

            #if the cluster is within the window size of a switch (provided in col ["Candidate_Switches"]), then we want to set the Switch_Proximal column to True for all SNPs in that cluster
            #df_SNPs_clusters["Switch_Proximal"] = df_SNPs_clusters.apply(lambda x: any([True if ((x["Cluster_Start"]  <= switch + window_size) & (switch - window_size <= x["Cluster_End"]) ) else False for switch in x["Candidate_Switches"]]), axis=1)
            df_SNPs_clusters["GC_Proximal"] = df_SNPs_clusters.apply(lambda x: any(((x["Cluster_Start"]  <= GC + window_size) & (GC - window_size <= x["Cluster_End"])) for GC in x["Candidate_GC"]), axis=1)
            logging.debug(df_SNPs_clusters.head())

            #now we want to set the Switch_Proximal column in the df_SNPs dataframe to True for all SNPs in clusters that are proximal to a switch
            #we do this by using the isin function to get all the Cluster_IDs that are proximal to a switch
            #then we use the isin function to set the Switch_Proximal column to True for all SNPs in those clusters
            current_df_SNPs = self.df_SNPs.copy()

            current_df_SNPs.loc[current_df_SNPs["Cluster_ID"].isin(df_SNPs_clusters[df_SNPs_clusters["GC_Proximal"] == True]["Cluster_ID"]), "GC_Proximal"] = True

            #now we want to set the Clustered_Proximal column to True for all SNPs in clusters that are proximal to a switch
            current_df_SNPs["Clustered_GC_Proximal"] = (current_df_SNPs["GC_Proximal"] == True) & (current_df_SNPs["Is_Clustered"] == True)

            #merge (outer) current_df_SNPs with self.df_SNPs, this will add the GC_Proximal and Clustered_GC_Proximal columns to self.df_SNPs
            self.df_SNPs = self.df_SNPs.merge(current_df_SNPs, left_index=True, right_index=True, how="outer", suffixes=("", "_y"))
            self.df_SNPs.drop(self.df_SNPs.filter(regex='_y$').columns, axis=1, inplace=True)
            self.df_SNPs["Is_Clustered"] = True
        
    def find_CO_Clusters(self, window_size: int) -> None:

        #get all clustered SNPs (i.e. SNPs that are in a cluster). That's what we're interested in
        self.df_SNPs = self.df_SNPs[self.df_SNPs["Is_Clustered"] == True]
        logging.info(f'{self.name} has {len(self.df_SNPs)} clustered SNVs')
        logging.info(f'{self.name} has {len(self.CO_dict)} chromosomes with COs.')

        #group df_SNPs by Cluster_ID, aggregate the min and max of the Region column, and the Chromosome column
        #this will give us the min and max positions of each cluster, and the chromosome it's on
        df_SNPs_clusters = self.df_SNPs.copy()

        # df_SNPs_clusters["CO_Proximal"] = False
        # df_SNPs_clusters["Clustered_CO_Proximal"] = False
        # df_SNPs_clusters["Window_Size_CO"] = window_size

        df_SNPs_clusters = df_SNPs_clusters.groupby("Cluster_ID").agg({"Region":["min","max"],"Chromosome":"first"}).reset_index()
        df_SNPs_clusters.columns = ["Cluster_ID","Cluster_Start","Cluster_End","Chromosome"]

        #cluster_start and cluster_end are the min and max positions of the cluster and should be integers
        df_SNPs_clusters["Cluster_Start"] = df_SNPs_clusters["Cluster_Start"].astype(int)
        df_SNPs_clusters["Cluster_End"] = df_SNPs_clusters["Cluster_End"].astype(int)
        df_SNPs_clusters["CO_Proximal"] = False
        
        logging.debug(df_SNPs_clusters.head())

        if df_SNPs_clusters.shape[0] == 0:
            #it is possible that the chromosome of the cluster is not in the switchesDict
            #if it is not, then we need to handle that
            logging.warning(f"No clustered SNPs found in {self.name}")
            self.df_SNPs["CO_Proximal"] = False
        
        else:
            #df_SNPs_clusters["Candidate_Switches"] = df_SNPs_clusters.apply(lambda x: switchesDict[x["Chromosome"]] if x["Chromosome"] in switchesDict.keys() else [], axis=1)
            df_SNPs_clusters["Candidate_CO"] = df_SNPs_clusters["Chromosome"].apply(lambda x: self.CO_dict.get(x, []))
            logging.debug(df_SNPs_clusters.head())

            #if the cluster is within the window size of a switch (provided in col ["Candidate_Switches"]), then we want to set the Switch_Proximal column to True for all SNPs in that cluster
            #df_SNPs_clusters["Switch_Proximal"] = df_SNPs_clusters.apply(lambda x: any([True if ((x["Cluster_Start"]  <= switch + window_size) & (switch - window_size <= x["Cluster_End"]) ) else False for switch in x["Candidate_Switches"]]), axis=1)
            df_SNPs_clusters["CO_Proximal"] = df_SNPs_clusters.apply(lambda x: any(((x["Cluster_Start"]  <= CO + window_size) & (CO - window_size <= x["Cluster_End"])) for CO in x["Candidate_CO"]), axis=1)
            logging.debug(df_SNPs_clusters.head())

            #now we want to set the Switch_Proximal column in the df_SNPs dataframe to True for all SNPs in clusters that are proximal to a switch
            #we do this by using the isin function to get all the Cluster_IDs that are proximal to a switch
            #then we use the isin function to set the Switch_Proximal column to True for all SNPs in those clusters
            current_df_SNPs = self.df_SNPs.copy()

            current_df_SNPs.loc[current_df_SNPs["Cluster_ID"].isin(df_SNPs_clusters[df_SNPs_clusters["CO_Proximal"] == True]["Cluster_ID"]), "CO_Proximal"] = True

            #now we want to set the Clustered_Proximal column to True for all SNPs in clusters that are proximal to a switch
            current_df_SNPs["Clustered_CO_Proximal"] = (current_df_SNPs["CO_Proximal"] == True) & (current_df_SNPs["Is_Clustered"] == True)

            #merge current_df_SNPs with self.df_SNPs, this will add the CO_Proximal and Clustered_CO_Proximal columns to self.df_SNPs
            self.df_SNPs = self.df_SNPs.merge(current_df_SNPs, left_index=True, right_index=True, how="outer", suffixes=("", "_y"))
            self.df_SNPs.drop(self.df_SNPs.filter(regex='_y$').columns, axis=1, inplace=True)
            self.df_SNPs["Is_Clustered"] = True

    def markLowSNPDensityClusters(self, lowDensityMap: dict, window_size: int = 1000) -> None:

        #if cluster edges + window_size in self.df_SNPs overlaps with a low density region, then mark it low_density
        #for each cluster, we need to get upper and lower bounds
        #then we need to check if those bounds +/- window_size overlap with a low density region

        self.df_SNPs["Low_Density"] = False
        self.df_SNPs["Low_Density_Window_Size"] = window_size

        #get all Cluster_IDs that are clustered
        clustered_SNPs = self.df_SNPs[self.df_SNPs["Is_Clustered"] == True]
        unique_cluster_ids = clustered_SNPs["Cluster_ID"].unique()

        #iterate through each unique cluster id
        for cluster_id in unique_cluster_ids:

            cluster = clustered_SNPs[clustered_SNPs["Cluster_ID"] == cluster_id]
            cluster_chromosome = cluster["Chromosome"].iloc[0]
            cluster_start = cluster["Cluster_Start"].iloc[0]
            cluster_end = cluster["Cluster_End"].iloc[0]

            upper_bound = cluster_end + window_size
            lower_bound = cluster_start - window_size

            #check if the upper and lower bounds overlap with a low density region
            lowDensityMap_chr = lowDensityMap[cluster_chromosome]

            for low_density_region in lowDensityMap_chr:
                if low_density_region[0] <= upper_bound and low_density_region[1] >= lower_bound:
                    self.df_SNPs.loc[(self.df_SNPs.Cluster_ID == cluster_id), 'Low_Density'] = True
                    break
        
    def findRandomClustersLegacy(self, window_size: int) -> None:
        #use if something goes wrong with the new way of doing it...

        self.df_SNPs["Random_Switch_Proximal"] = False
        self.df_SNPs["RandomSwitches_Window_Size"] = window_size

        for index, row in self.randomSwitchesTable.iterrows():
            chromosome = row["Chromosome"]
            switch = row["Switch_Center"]
            self.df_SNPs["Random_Switch_Proximal"] = self.df_SNPs["Random_Switch_Proximal"] | (self.df_SNPs["Region"].between(switch-window_size,switch+window_size) & (self.df_SNPs["Chromosome"] == chromosome))

        for snp_index, snp_row in self.df_SNPs.iterrows():
            if snp_row["Is_Clustered"] == True:
                if snp_row["Random_Switch_Proximal"] == True:
                    cluster_id = snp_row["Cluster_ID"]
                    self.df_SNPs.loc[(self.df_SNPs.Cluster_ID == cluster_id), 'Random_Switch_Proximal'] = True
        
        self.df_SNPs["Random_Clustered_Proximal"] = False
        self.df_SNPs["Random_Clustered_Proximal"] = (self.df_SNPs["Random_Switch_Proximal"] == True) & (self.df_SNPs["Is_Clustered"] == True)

    def findRandomClusters(self, window_size: int) -> None:

        #get all clustered SNPs (i.e. SNPs that are in a cluster). That's what we're interested in
        self.df_SNPs = self.df_SNPs[self.df_SNPs["Is_Clustered"] == True]

        #get the Chromosome and Switch_Center columns from the randomSwitchesTable, convert to a dictionary
        #where the keys are the chromosome names and the values are the switch centers
        randomSwitchesDict = self.randomSwitchesTable.groupby("Chromosome")["Switch_Center"].apply(list).to_dict()

        #group df_SNPs by Cluster_ID, aggregate the min and max of the Region column, and the Chromosome column
        #this will give us the min and max positions of each cluster, and the chromosome it's on
        df_SNPs_clusters = self.df_SNPs.groupby("Cluster_ID").agg(
            Cluster_Start=("Region", "min"),
            Cluster_End=("Region", "max"),
            Chromosome=("Chromosome", "first")
        ).reset_index()

        #by default, all SNPs are not proximal to a switch/recombination event
        df_SNPs_clusters["Random_Switch_Proximal"] = False
        self.df_SNPs["Random_Switch_Proximal"] = False
        self.df_SNPs["Random_Clustered_Proximal"] = False
        self.df_SNPs["RandomSwitches_Window_Size"] = window_size

        #cluster_start and cluster_end are the min and max positions of the cluster and should be integers
        df_SNPs_clusters["Cluster_Start"] = df_SNPs_clusters["Cluster_Start"].astype(int)
        df_SNPs_clusters["Cluster_End"] = df_SNPs_clusters["Cluster_End"].astype(int)

        #make sure that df_SNPs_clusters["Chromosome"] is in the randomSwitchesDict
        #if it's not, then we don't need to include it (otherwise we'll get an error when we try to get the switch center for that chromosome)
        df_SNPs_clusters = df_SNPs_clusters[df_SNPs_clusters["Chromosome"].isin(randomSwitchesDict.keys())]
        
        #if the chromosome of the cluster is in the randomSwitchesDict, then we want to check if the switch center is within the window size of the cluster's min and max positions
        #if it is, then we want to set the Random_Switch_Proximal column to True for all SNPs in that cluster
        #first, we need to create a new column in df_SNPs_clusters that will hold the switch centers for each cluster
        #we do this by using the randomSwitchesDict to get the switch center for the chromosome of the cluster
        #then we use the apply function to apply the lambda function to each row of the df_SNPs_clusters dataframe
        #the lambda function takes the switch center and the cluster start and end positions, and checks if the switch center is within the window size of the cluster start and end positions
        #if it is, then it returns True, otherwise it returns False
        #do it only if df_SNPs_clusters has rows in it
        
        if df_SNPs_clusters.empty:
            logging.warning(f"No clustered SNPs found in {self.name}")
        
        else:
            #it is possible that the chromosome of the cluster is not in the switchesDict
            #if it is not, then we need to handle that
            df_SNPs_clusters["Candidate_Switches"] = df_SNPs_clusters.apply(lambda x: randomSwitchesDict[x["Chromosome"]] if x["Chromosome"] in randomSwitchesDict.keys() else [], axis=1)

            #if the cluster is within the window size of a switch (provided in col ["Candidate_Switches"]), then we want to set the Switch_Proximal column to True for all SNPs in that cluster
            
            # old method
            #df_SNPs_clusters["Random_Switch_Proximal"] = df_SNPs_clusters.apply(lambda x: any([True if ((x["Cluster_Start"]  <= switch + window_size) & (switch - window_size <= x["Cluster_End"]) ) else False for switch in x["Candidate_Switches"]]), axis=1)
            # new method
            def is_proximal(cluster_start, cluster_end, switches):
                return any((cluster_start <= switch + window_size) & (switch - window_size <= cluster_end) for switch in switches)

            df_SNPs_clusters["Random_Switch_Proximal"] = df_SNPs_clusters.apply(
                lambda x: is_proximal(x["Cluster_Start"], x["Cluster_End"], x["Candidate_Switches"]), axis=1
            )

            #now we want to set the Random_Switch_Proximal column in the df_SNPs dataframe to True for all SNPs in clusters that are proximal to a switch
            #we do this by using the isin function to get all the Cluster_IDs that are proximal to a switch
            #then we use the isin function to set the Random_Switch_Proximal column to True for all SNPs in those clusters
            
            #old and new method
            #self.df_SNPs.loc[self.df_SNPs["Cluster_ID"].isin(df_SNPs_clusters[df_SNPs_clusters["Random_Switch_Proximal"] == True]["Cluster_ID"]), "Random_Switch_Proximal"] = True
            proximal_clusters = df_SNPs_clusters[df_SNPs_clusters["Random_Switch_Proximal"]]["Cluster_ID"]
            self.df_SNPs.loc[self.df_SNPs["Cluster_ID"].isin(proximal_clusters), "Random_Switch_Proximal"] = True

            #now we want to set the Random_Clustered_Proximal column to True for all SNPs in clusters that are proximal to a switch
            self.df_SNPs["Random_Clustered_Proximal"] = (self.df_SNPs["Random_Switch_Proximal"] == True) & (self.df_SNPs["Is_Clustered"] == True)

    def output_df_table(self) -> pd.DataFrame:
        output_df = self.df
        output_df["Percent"] = self.percent
        return output_df
    
    def classify_GC_CO(self, threshold_GC_max: int = 5000):
        """
        Classifies genetic crossover (CO) and gene conversion (GC) events based on the given threshold.

        This method uses the get_GC_CO function from the mayo module to classify the SNPs in the df_SNPs DataFrame 
        into gene conversion (GC) and crossover (CO) events. The classification is based on a given threshold for 
        the maximum distance for GC events. The results are stored in the GC_dict and CO_dict attributes of the 
        instance.

        Parameters:
        threshold_GC_max (int, optional): The maximum distance for classifying a SNP as a GC event. Defaults to 5000.

        Returns:
        None
        """
        from mayo import get_GC_CO, get_haplotype_lengths

        self.haplotype_lengths = get_haplotype_lengths(self.switchesTable)
        self.GC_dict, self.CO_dict = get_GC_CO(self.haplotype_lengths, threshold_GC_max=threshold_GC_max)
    
    def calculateTotalssDNA(self, cluster_kind: Literal["JT", "PMACD"] = "JT"):
        """
        Prints the total length of ssDNA in the OverlapSwitchesTable object for each p-value.

        Returns:
            A list of lists containing the name of the OverlapSwitchesTable object, the p-value, and the total length of ssDNA.
        """
        output_list = []

        if cluster_kind == "JT":
            for pval_key in self.cluster_dict_JT.keys():
                df = self.cluster_dict_JT[pval_key]
                total_ssDNA = df["Length"].sum()
                self.totalssDNA_JT_dict[pval_key] = total_ssDNA

                logging.info(f"Total ssDNA in {self.name}, pval={pval_key}): {total_ssDNA} bp")
                output_list.append([self.name, pval_key, total_ssDNA])

        elif cluster_kind == "PMACD":
            for pval_key in self.cluster_dict_PMACD.keys():
                df = self.cluster_dict_PMACD[pval_key]
                total_ssDNA = df["Cluster_Length"].sum()
                self.totalssDNA_PMACD_dict[pval_key] = total_ssDNA

                logging.info(f"Total ssDNA in {self.name}, pval={pval_key}): {total_ssDNA} bp")
                output_list.append([self.name, pval_key, total_ssDNA])

        
        return output_list
            
