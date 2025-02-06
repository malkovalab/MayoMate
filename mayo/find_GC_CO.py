# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd

def calculate_distance(switch_centers: list) -> list:
    """
    Calculates the distance between switch centers.

    Args:
        switch_centers (list): A list of integers representing the positions of switch centers.

    Returns:
        list: A list of tuples, where each tuple contains the distance between two switch centers and the positions of those switch centers.
    """
    #make sure switch_centers are sorted in ascending order
    switch_centers.sort()

    list_of_distances = []

    if len(switch_centers) == 1:
    #length 0 means it is a last block, which will be a CO event
        list_of_distances.append((0, (switch_centers[0], switch_centers[0])))

    for i in range(0, len(switch_centers)-1):

        position_start = switch_centers[i]
        position_end = switch_centers[i+1]

        distance = position_end - position_start
        
        list_of_distances.append((distance, (position_start, position_end)))
    
    return list_of_distances

def get_haplotype_lengths(sample_df) -> dict:
    """
    Calculates the distance between switch centers.

    Args:
        switch_centers (list): A list of integers representing the positions of switch centers.

    Returns:
        list: A list of tuples, where each tuple contains the distance between two switch centers and the positions of those switch centers.
    """
    haplotype_blocks = {}

    #first, iterate through each chromosome in the dataframe and sort by Switch_Center
    for chromosome in sample_df.Chromosome.unique():
        df_chrom = sample_df.query('Chromosome == @chromosome')

        df_chrom_sorted = df_chrom.sort_values(by=['Switch_Center'])

        #store df_chrom_sorted Switch_Center column as a list
        switch_centers = df_chrom_sorted['Switch_Center'].tolist()

        #second, iterate through each Switch_Center in the sorted df and calculate the distance between each Switch_Center
        list_of_distances = calculate_distance(switch_centers)

        add_to_dict = {chromosome: list_of_distances}
        haplotype_blocks.update(add_to_dict)

    return haplotype_blocks

def update_event_dictionary(event_dictionary: dict, chromosome: str, event_pos: tuple) -> dict:
    """
    Updates the event dictionary with a new event position for a given chromosome.

    Args:
        event_dictionary (dict): A dictionary containing the event positions for each chromosome.
        chromosome (str): The name of the chromosome to update.
        event_pos (tuple): A tuple containing the start and end positions of the event.

    Returns:
        dict: The updated event dictionary.
    """
    if chromosome in event_dictionary.keys():
        event_dictionary[chromosome].append(event_pos)
    else:
        event_dictionary.update({chromosome: [event_pos]})

    return event_dictionary

def recalculate_distances(chromosome_data, shortest_block) -> list:
    """
    Recalculates the distances between haplotype blocks after removing the shortest block.

    Args:
        chromosome_data (list): A list of tuples, where each tuple contains the length of a haplotype block and the positions of the switch centers.
        shortest_block (tuple): A tuple containing the length of the shortest haplotype block and the positions of the switch centers.

    Returns:
        list: A list of tuples, where each tuple contains the distance between two switch centers and the positions of those switch centers.
    """
    #first need to extract the Switch_Center positions from the chromosome_data
    switch_centers = []

    for block in chromosome_data:
        switch_centers.append(block[1][0])
        switch_centers.append(block[1][1])

    #remove duplicates
    switch_centers = list(dict.fromkeys(switch_centers))

    #remove the shortest block bounds from the switch_centers
    switch_centers.remove(shortest_block[1][0])
    switch_centers.remove(shortest_block[1][1])

    #sort the switch_centers
    switch_centers.sort()

    #recompute the distances between the remaining haplotype blocks
    chromosome_data = calculate_distance(switch_centers)

    return chromosome_data

def get_GC_CO(haplotype_blocks: dict, threshold_GC_max = 1000) -> tuple:
    """
    Calculates the positions of CO and GC events for each chromosome in the input haplotype blocks.

    Args:
        haplotype_blocks (dict): A dictionary where the keys are chromosome names and the values are lists of tuples, where each tuple contains the distance between two switch centers and the positions of those switch centers.
        threshold_GC_max (int): The maximum length of a haplotype block for it to be considered a GC event. Defaults to 1000.

    Returns:
        tuple: A tuple containing two dictionaries. The first dictionary contains the positions of GC events for each chromosome, and the second dictionary contains the positions of CO events for each chromosome.
    """
    GC_dictionary = {}
    CO_dictionary = {}

    for chr_key in haplotype_blocks.keys():

        chromosome_data = haplotype_blocks[chr_key]

        while len(chromosome_data) > 0:

            #if there is only one block left and it is of length 0, it is 1 residual CO event
            if len(chromosome_data) == 1 and chromosome_data[0][0] == 0:
                assert chromosome_data[0][1][0] == chromosome_data[0][1][1], f"chromosome_data[0][1][0]: {chromosome_data[0][1][0]}, chromosome_data[0][1][1]: {chromosome_data[0][1][1]}"
                CO_dictionary = update_event_dictionary(CO_dictionary, chr_key, chromosome_data[0][1][0])
                break

            #find the shortest haplotype block
            shortest_block = min(chromosome_data, key=lambda x: x[0])

            if shortest_block[0] > threshold_GC_max:
                #it is a CO event
                CO_dictionary = update_event_dictionary(CO_dictionary, chr_key, shortest_block[1][0])
                CO_dictionary = update_event_dictionary(CO_dictionary, chr_key, shortest_block[1][1])
            else:
                #it is a GC event
                GC_center = (shortest_block[1][0] + shortest_block[1][1]) / 2
                GC_dictionary = update_event_dictionary(GC_dictionary, chr_key, GC_center)

            #recompute the distances with the shortest block removed
            chromosome_data = recalculate_distances(chromosome_data, shortest_block)

    return GC_dictionary, CO_dictionary
