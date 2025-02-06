from mayo import *
from mayo.settings import config as cfg
import logging
import multiprocessing

logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s:.%(funcName)s: %(message)s')

def prepare_data_for_simulation(saved_state_path):
    
    samples_SNVs, switches_list, df_clusters, chromosomes = load_payload_from_pickle(saved_state_path)

    #samples for simulation:
    haploid_samples_ung1d = [
        'Mal_JK11',
        'Mal_JK14',
        'Mal_JK16', 
        'Mal_JK17', 
        'Mal_JK19', 
        'Mal_JK21', 
        'Mal_JK23', 
        'Mal_JK24', 
        'Mal_JK26', 
        'Mal_JK27',
        'Mal_JK29', 
        'Mal_JK30', 
        'Mal_JK31', 
        #'Mal_JK35', No clusters detected, so no need to simulate
        'Mal_JK5-', 
        'Mal_JK58', 
        'Mal_JK59', 
        'Mal_JK6-', 
        'Mal_JK62', 
        'Mal_JK63', 
        'Mal_JK64', 
        'Mal_JK66', 
        'Mal_JK68',
        'Mal_JK7-', 
        'Mal_JK8-',
        'Mal_JK9-', 
        'Mal_JK96', 
        'JT31_CKD',
        'JT36_CKD', 
        'JT37_CKD', 
        'JK2_10B_', 
        'JK2_11B_', 
        'JK2_12B_', 
        'JK2_13_G', 
        'JK2_14_A', 
        'JK2_15_C', 
        'JK2_19B_', 
        'JK2_1_CT',
        'JK2_2_TA',
        'JK2_56_C', 
        'JK2_57_A', 
        'JK2_58_A', 
        'JK2_59_C', 
        'JK2_5_GA', 
        'JK2_60_A', 
        'JK2_61_C', 
        'JK2_62_C', 
        'JK2_63_A', 
        'JK2_64_G', 
        'JK2_65_C', 
        'JK2_68_A',
        'JK2_71_C',
        'JK2_72_G', 
        'JK2_73_G', 
        'JK2_74_G', 
        #'JK2_75_A', no clusters detected, so no need to simulate
        'JK2_77_C', 
        'JK2_78_T', 
        'JK2_7_GT', 
        'JK2_9B_C',

        'JMT49_CK',
        'JMT50_CK',
        'JMT51_CK',
        'JMT52_CK',
        'JMT53_CK',
        'JMT54_CK',
        'JMT55_CK',
        'JMT56_CK',
        
        'JT31_CKD',
        #'JT34_CKD', No clusters detected, so no need to simulate
        'JT36_CKD',
        'JT37_CKD',
        
        ]

    samples_SNVs[:] = [i for i in samples_SNVs if i.name in (haploid_samples_ung1d)]
    logging.info("Created a subset of samples for simulation (samples_SNVs)")

    switches_list[:] = [i for i in switches_list if i.name in (haploid_samples_ung1d)]
    logging.info("Created a subset of samples for simulation (switches_list)")

    for sample in samples_SNVs:
        if sample.aneuploid == True:
            logging.warning(f"{sample.name} is aneuploid, but it was included in the simulation dataset. Please check the data.")

    total_clusters = 0
    for i in samples_SNVs:
        i_clusters = i.df_SNPs["Cluster_ID"].nunique()
        total_clusters += i_clusters
        logging.info(f"{i.name} has {i_clusters} clusters")
    logging.info(f"######## Total number of clusters in the simulation dataset: {total_clusters} ########")

    return samples_SNVs, switches_list, df_clusters, chromosomes

def run_simulation(switches_list, chromosomes, samples_SNVs, simulation_generations, save_file_name):
    conductNullHypothesisSimulation(
        switches_list=switches_list, 
        chromosomes=chromosomes, 
        frames_list_samples_normal_ref=samples_SNVs, 
        simulation_generations=simulation_generations, 
        window_params=(0, 50001, 500),
        save_file_name=save_file_name)
    logging.info(f'Finished writing everything to the file: {save_file_name}')
    logging.info('Completed null hypothesis simulation')

if __name__ == "__main__":

    # before starting, save the state of samples_SNVs, switches_list, and df_clusters, and chromosomes
    # save_payload_to_pickle([samples_SNVs, switches_list, df_clusters, chromosomes], "saved_state_2snps.pkl")

    SNP_NUM = 1
    subprocesses=8
    simulation_generations_per_process = 1
    saved_state_path = f"outputs/partial_pickles/saved_state_for_sim_JT_10000_0.0001_{SNP_NUM}snp.pkl" 

    samples_SNVs, switches_list, df_clusters, chromosomes = prepare_data_for_simulation(saved_state_path)

    processes = []

    for i in range(subprocesses):
        p = multiprocessing.Process(target=run_simulation, args=(switches_list, chromosomes, samples_SNVs, simulation_generations_per_process, f"proximal_clusters_{simulation_generations_per_process}_{SNP_NUM}snp_sub_{i}.csv"))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    logging.info("Finished all subprocesses")

    #open all the files and combine them into one
    for i in range(subprocesses):
        if i == 0:
            df = pd.read_csv(f"proximal_clusters_{simulation_generations_per_process}_{SNP_NUM}snp_sub_{i}.csv")
        else:
            df2 = pd.read_csv(f"proximal_clusters_{simulation_generations_per_process}_{SNP_NUM}snp_sub_{i}.csv")
            #remove the Generation==0 rows
            #df2 = df2[df2.Generation != 0]
            #add the (subprocess_index * simulation_generations_per_process) to the Generation column
            df2['Generation'] = df2['Generation'] + (i * simulation_generations_per_process)

            df = pd.concat([df, df2])

    #transfer generation0 values in Associated_Clusters and Other_Clusters columns to all generations based on the Window_Size column
    gen0 = df[df.Generation == 0].copy()
    gen0 = gen0[['Window_Size', 'Associated_Clusters', 'Other_Clusters']]
    df = df.merge(gen0, on='Window_Size', how='left', suffixes=('', '_gen0'))
    df['Associated_Clusters'] = df['Associated_Clusters_gen0']
    df['Other_Clusters'] = df['Other_Clusters_gen0']
    df.drop(['Associated_Clusters_gen0', 'Other_Clusters_gen0'], axis=1, inplace=True)
    df.to_csv(f"proximal_clusters_full_{simulation_generations_per_process}_{SNP_NUM}snp_combined.csv", index=False)