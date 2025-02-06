from mayo import *
import pandas as pd
import logging
import numpy as np
import time
from multiprocessing import Process
from mayo.settings import chromosomes as chromosomes
from mayo.settings import haploid_samples as haploid_samples
from mayo.settings import haploid_samples_ung1d as haploid_samples_ung1d

# Set Reference Chromosomal Intervals for figure drawing and create a DataFrame
from mayo.settings import chromosomes as chromosomes
df_chr = pd.DataFrame(chromosomes, columns=["chromosome", "end_position"])

def get_chromosomes_cumulative_sum(chromosomes: list) -> tuple[int, dict]:
    """
    This function calculates the genome size and the cumulative sum of the chromosomes.

    Parameters:
    chromosomes (list): A list of tuples where each tuple represents a chromosome. The first element of the tuple is the chromosome name and the second element is the chromosome size.

    Returns:
    genome_size (int): The total size of the genome.
    chromosomes_cumulative_sum (dict): A dictionary where the keys are the chromosome names and the values are tuples representing the range of the chromosome in the genome.
    """
    genome_size = 0
    chromosomes_cumulative_sum = {}

    for chr in chromosomes:
        current_chromosome = chr[0]
        genome_size += chr[1]
        chromosome_range = (genome_size - chr[1], genome_size)
        chromosomes_cumulative_sum[current_chromosome] = chromosome_range

    return genome_size, chromosomes_cumulative_sum

def get_chromosome_from_random_position(random_position, chromosomes_cumulative_sum):
    """
    This function returns the chromosome that contains the given random position.

    Parameters:
    random_position (int): A random position in the genome.
    chromosomes_cumulative_sum (dict): A dictionary where the keys are the chromosome names and the values are tuples representing the range of the chromosome in the genome.

    Returns:
    random_chromosome (str): The name of the chromosome that contains the random position.
    """
    for chromosome, chromosome_range in chromosomes_cumulative_sum.items():
        if chromosome_range[0] <= random_position <= chromosome_range[1]:
            random_chromosome = chromosome
            break

    return random_chromosome

def generate_random_SNVs_positions(samples_SNVs, number_of_simulations = 1000, simulation_index_add = 0) -> None:
    """
    This function generates random SNVs positions for each sample.

    Parameters:
    samples_SNVs (list): A list of sample objects. Each sample object should have a 'df' attribute which is a DataFrame containing the SNVs for that sample.
    number_of_simulations (int, optional): The number of simulations to run. Default is 1000.

    Returns:
    None. The function modifies the sample objects in place, adding a 'df_random' attribute to each one which is a DataFrame containing the random SNVs positions.
    """
    import random
    genome_size, chromosomes_cumulative_sum = get_chromosomes_cumulative_sum(chromosomes)

    #creseed the random state
    #np.random.seed() #this doesn't work for some reason :(
    random_float = random.random()
    new_random_state = int(random_float * 1000000)
    np.random.seed(new_random_state)
    print(f"new subprocess numpy random state: {np.random.get_state()[1][0]}")

    for sample in samples_SNVs:
        #print(f"Generating random SNVs positions for sample {sample.name}...")

        total_mutations = len(sample.df)
        if total_mutations == 0:
            #print(f"Sample {sample.name} has no mutations. Skipping...")
            sample.df_random = pd.DataFrame(columns=["simulation_index", "Chromosome", "Region"])
            continue

        #this will be used for multiprocessing, let's make sure that the random state is different for each process
        #print the current numpy random state

        simulation_indices = np.repeat(np.arange(number_of_simulations), total_mutations)
        random_positions = np.random.randint(0, genome_size, total_mutations * number_of_simulations)
        random_chromosomes = np.vectorize(get_chromosome_from_random_position)(random_positions, chromosomes_cumulative_sum)

        #add simulation_index_add to the simulation_indices
        simulation_indices = simulation_indices + simulation_index_add

        sample.df_random = pd.DataFrame({
            "simulation_index": simulation_indices,
            "Chromosome": random_chromosomes,
            "Region": random_positions
        })

def update_random_cluster_dict_JT(sample, pval, add_dict):
    """
    This function updates the random cluster dictionary for a given sample and p-value.

    Parameters:
    sample (object): The sample object to update.
    pval (float): The p-value to use as a key in the dictionary.
    add_dict (dict): The dictionary to add to the sample's random cluster dictionary.

    Returns:
    None. The function updates the sample's random cluster dictionary in place.
    """
    #retrieve the random clusters for the given sample and pval
    random_clusters = sample.random_cluster_dict_JT[pval]

    #add in the new clusters while preserving the existing clusters
    random_clusters.update(add_dict)
    
    #add the sample's random_cluster_dict_JT
    sample.random_cluster_dict_JT[pval] = random_clusters

def find_clusters_in_sample_simulated_data(
        sample,
        pval,
        max_distance_between_mutations,
        min_cluster_mutations,
        simulation_index):
    """
    This function finds mutation clusters in simulated data for a given sample and p-value.

    Parameters:
    sample (object): The sample object to update.
    pval (float): The p-value to use for finding clusters.
    max_distance_between_mutations (int): The maximum distance between mutations to consider a cluster.
    min_cluster_mutations (int): The minimum number of mutations to consider a cluster.
    simulation_index (int): The index of the simulation run.

    Returns:
    random_clusters (DataFrame): A DataFrame containing the random mutation clusters found in the sample's simulated data.
    """
    from mayo.nbinom import analyzeSampleForMutationClusters

    df_random_mutations = sample.df_random[sample.df_random["simulation_index"] == simulation_index].copy()
    
    random_clusters = analyzeSampleForMutationClusters(
        df_random_mutations,
        min_p_value = pval,
        max_distance_between_mutations = max_distance_between_mutations,
        min_cluster_mutations = min_cluster_mutations)
    
    #update a nested dictionary with pval as primary key and simulation index as secondary key, and add the random clusters as a value
    update_random_cluster_dict_JT(
        sample,
        pval,
        {simulation_index: random_clusters})
    
    random_clusters["sample"] = sample.name
    random_clusters["all_mutations"] = len(sample.df)
    random_clusters["simulation_index"] = simulation_index

    return random_clusters

def compute_and_add_JT_clusters_random(
        samples_SNVs, 
        pvals_to_test, 
        max_distance_between_mutations = 10000, 
        min_cluster_mutations = 3):
    """
    This function computes and adds JT clusters for a list of samples and a list of p-values to test.

    Parameters:
    samples_SNVs (list): A list of sample objects. Each sample object should have a 'df_random' attribute which is a DataFrame containing the random SNVs for that sample.
    pvals_to_test (list): A list of p-values to test.
    max_distance_between_mutations (int, optional): The maximum distance between mutations to consider. Default is 10000.
    min_cluster_mutations (int, optional): The minimum number of mutations to consider a cluster. Default is 3.

    Returns:
    None. The function updates the sample objects in place, adding a 'random_cluster_dict_JT' attribute to each one which is a dictionary containing the JT clusters.
    """
    for sample in samples_SNVs:
        sample.random_cluster_dict_JT = {}

    for pval in pvals_to_test:
        #print(f"Testing pval = {pval}")
        pval_master_df_random = pd.DataFrame()

        for sample in samples_SNVs:
            df_random_mutations = sample.df_random.copy()
            #create a dictionary with simulation index as key and random clusters as value
            sample.random_cluster_dict_JT.update({pval: {}})

            for simulation_index in df_random_mutations["simulation_index"].unique():
            
                random_clusters = find_clusters_in_sample_simulated_data(sample, pval, max_distance_between_mutations, min_cluster_mutations, simulation_index)
                pval_master_df_random = pd.concat([pval_master_df_random, random_clusters])
                pval_master_df_random.reset_index(drop=True, inplace=True) #reset the index so that it is unique
                logging.debug(f"{sample.name} has {len(sample.df):3} mutations and {len(random_clusters):3} random JT clusters")

            #ID each row (cluster) by index into "Dataset_Cluster_ID" column
            pval_master_df_random["Dataset_Cluster_ID"] = pval_master_df_random.index
            #since the first index is 0, add 1 to each index to avoid confusion with 0
            pval_master_df_random["Dataset_Cluster_ID"] = pval_master_df_random["Dataset_Cluster_ID"] + 1

        pval_master_df_random.to_csv(f"data/new_clusters/mutation_clusters/random_mutation_clusters_pval_{str(pval)}.tsv", sep="\t", index=False)

def run_simulation_process(
        process_index, 
        simulations_per_process, 
        pvals_to_test, 
        payload_Samples_SNVs, 
        max_distance_between_mutations=10000, 
        min_cluster_mutations=3):

    samples_SNVs = load_payload_from_pickle(payload_Samples_SNVs)
    
    generate_random_SNVs_positions(
        samples_SNVs, 
        number_of_simulations=simulations_per_process, 
        simulation_index_add = process_index*simulations_per_process)

    compute_and_add_JT_clusters_random(
        samples_SNVs,
        pvals_to_test,
        max_distance_between_mutations=max_distance_between_mutations,
        min_cluster_mutations=min_cluster_mutations)

    slice_save_name = f"null_distribution_of_SNVs_slice{process_index}.pkl"
    #save the lusters to an output csv file
    slice_output_dictionary = {}
    for sample in samples_SNVs:
        random_JT_cluster_dict_slice = sample.random_cluster_dict_JT
        slice_output_dictionary[sample.name] = random_JT_cluster_dict_slice
    
    save_payload_to_pickle(slice_output_dictionary, slice_save_name)

    print(f"Process {process_index} finished successfully..., saved to {slice_save_name}")
    
def update_combined_dictionary(all_slices_dict, slice_dict) -> dict:
    """
    This function updates the combined dictionary with the data from a slice dictionary.

    Parameters:
    all_slices_dict (dict): The combined dictionary to update.
    slice_dict (dict): The slice dictionary to add to the combined dictionary.

    Returns:
    None. The function updates the combined dictionary in place.
    """
    for sample, random_clusters in slice_dict.items():
        if sample in all_slices_dict:

            current_sample_random_clusters = all_slices_dict[sample]
            #current_sample_random_clusters contains a dictionary with pval as key and a dictionary with simulation index as key and random clusters as value

            for pval, simulation_index_dict in random_clusters.items():
                if pval in current_sample_random_clusters:
                    current_sample_pval_dict = current_sample_random_clusters[pval]
                    #add the new simulation indices to the existing simulation indices

                    for simulation_index, random_clusters in simulation_index_dict.items():
                        current_sample_pval_dict.update({simulation_index: random_clusters})

                    #update the combined dictionary with the new simulation indices
                    current_sample_random_clusters.update({pval: current_sample_pval_dict})

                else:
                    #add the new pval to the existing pvals
                    current_sample_random_clusters.update({pval: simulation_index_dict})
            
            all_slices_dict[sample] = current_sample_random_clusters

        else:
            all_slices_dict[sample] = random_clusters

    return all_slices_dict

def collect_pickle_slices(
        pkl_files_key="null_distribution_of_SNVs_slice*.pkl", 
        output_pkl_file_name="null_distribution_of_SNVs.pkl"):
    """
    This function joins the pickle files created by the simulation process.

    Parameters:
    None.

    Returns:
    None. The function saves the combined dictionary to a pickle file.
    """
    import glob
    import pickle

    print("Joining pickle files...")

    all_files = glob.glob(pkl_files_key)
    all_files.sort()

    print(f"Found {len(all_files)} pickle files: {all_files}")

    all_slices = []
    for file in all_files:
        all_slices.append(load_payload_from_pickle(file))

    save_payload_to_pickle(all_slices, output_pkl_file_name)

def clean_up_pickle_slices(pkl_files_key="null_distribution_of_SNVs_slice*.pkl") -> None:
    """
    This function deletes the pickle files created by the simulation process.

    Parameters:
    None.

    Returns:
    None.
    """
    import glob
    import os

    print("Cleaning up pickle files...")

    all_files = glob.glob(pkl_files_key)
    all_files.sort()

    for file in all_files:
        os.remove(file)

if __name__ == "__main__":

    #SETTINGS
    genome_size, chromosomes_cumulative_sum = get_chromosomes_cumulative_sum(chromosomes)
    simulations_per_process = 2
    processes = 6
    payload_Samples_SNVs = "outputs/partial_pickles/samples_SNVs.pkl"
    output_pickle_file_name = "outputs/partial_pickles/combined_simulations_test_new.pkl"

    pvals_to_test = [0.01, 0.001, 0.0001, 0.00001, 0.000001]
    max_distance_between_mutations=10000, 
    min_cluster_mutations=3
    #END SETTINGS


    time_start = time.time()
    processes_list = []

    for i in range(processes):
        p = Process(target=run_simulation_process, args=(
            i, 
            simulations_per_process, 
            pvals_to_test, 
            payload_Samples_SNVs,
            max_distance_between_mutations,
            min_cluster_mutations))
        p.start()
        processes_list.append(p)

    print(f"Scheduled {len(processes_list)} processes... Waiting for them to finish...")

    for p in processes_list:
        p.join()
    
    print("All sub-processes finished successfully...")

    collect_pickle_slices(output_pkl_file_name="outputs/partial_pickles/null_distribution_of_SNVs.pkl")
    clean_up_pickle_slices(pkl_files_key="null_distribution_of_SNVs_slice*.pkl")

    print(f"Done! Total time: {time.time() - time_start} seconds for {processes} processes with {simulations_per_process} simulations each.")


    all_slices = load_payload_from_pickle("outputs/partial_pickles/null_distribution_of_SNVs.pkl")

    sample_names = list(all_slices[0].keys())

    combined_simulations = {}
    for sample_name in sample_names:
        sample_dict_combined = {}
        for slice in all_slices:
            slice_item = slice[sample_name]
            
            sample_dict_combined = update_combined_dictionary(sample_dict_combined, slice_item)

        combined_simulations[sample_name] = sample_dict_combined

    save_payload_to_pickle(combined_simulations, output_pickle_file_name)

    print(f"Done! Saved combined simulations to {output_pickle_file_name}")