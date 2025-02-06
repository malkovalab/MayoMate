# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
import pandas as pd
from numpy import NaN

chromosomes: list[list] = [
    ['ref|NC_001133|',230218], #1.91%
    ['ref|NC_001134|',813184], #6.74%
    ['ref|NC_001135|',319343], #2.64%
    ['ref|NC_001136|',1531933],#12.69%
    ['ref|NC_001137|',576070], #4.77%
    ['ref|NC_001138|',270161], #2.24%
    ['ref|NC_001139|',1090940],#9.03%
    ['ref|NC_001140|',562643], #4.66%
    ['ref|NC_001141|',439888], #3.64%
    ['ref|NC_001142|',745751], #6.17%
    ['ref|NC_001143|',666816], #5.52%
    ['ref|NC_001144|',1078177],#8.93%
    ['ref|NC_001145|',924431], #7.65%
    ['ref|NC_001146|',784333], #6.49%
    ['ref|NC_001147|',1091291],#9.03%
    ['ref|NC_001148|',948066]] #7.85%

centromeres = pd.DataFrame([
    ['chr1', 151582, "cen", NaN, NaN, NaN, "cen"],
    ['chr2',238323,  "cen", NaN, NaN, NaN, "cen"],
    ['chr3',114501,  "cen", NaN, NaN, NaN, "cen"],
    ['chr4',449821,  "cen", NaN, NaN, NaN, "cen"],
    ['chr5',152104,  "cen", NaN, NaN, NaN, "cen"],
    ['chr6',148627,  "cen", NaN, NaN, NaN, "cen"],
    ['chr7',497038,  "cen", NaN, NaN, NaN, "cen"],
    ['chr8',105703,  "cen", NaN, NaN, NaN, "cen"],
    ['chr9',355745,  "cen", NaN, NaN, NaN, "cen"],
    ['chr10',436425, "cen", NaN, NaN, NaN, "cen"],
    ['chr11',440246, "cen", NaN, NaN, NaN, "cen"],
    ['chr12',150947, "cen", NaN, NaN, NaN, "cen"],
    ['chr13',268149, "cen", NaN, NaN, NaN, "cen"],
    ['chr14',628875, "cen", NaN, NaN, NaN, "cen"],
    ['chr15',326702, "cen", NaN, NaN, NaN, "cen"],
    ['chr16',556073, "cen", NaN, NaN, NaN, "cen"]], 
    columns=["Chromosome", "Region", "Reference", "Allele", "Frequency", "Zygosity", "SPECTRA_STRANDWISE"])

telomeres = [
    #chromosome name, length of left telomere, length of right telomere
    ['ref|NC_001133|', 800,   808],
    ['ref|NC_001134|', 6608,  806],
    ['ref|NC_001135|', 1098,  838],
    ['ref|NC_001136|', 904,   7309],
    ['ref|NC_001137|', 6473,  7276],
    ['ref|NC_001138|', 5530,  431],
    ['ref|NC_001139|', 781,   7306],
    ['ref|NC_001140|', 5505,  6539],
    ['ref|NC_001141|', 7784,  821],
    ['ref|NC_001142|', 7767,  850],
    ['ref|NC_001143|', 807,   913],
    ['ref|NC_001144|', 12085, 13897],
    ['ref|NC_001145|', 6344,  891],
    ['ref|NC_001146|', 7428,  1056],
    ['ref|NC_001147|', 847,   7370],
    ['ref|NC_001148|', 7223,  5615]]

numeric_to_roman: dict[str,str] = {
    "chr1": "I",
    "chr2": "II",
    "chr3": "III",
    "chr4": "IV",
    "chr5": "V",
    "chr6": "VI",
    "chr7": "VII",
    "chr8": "VIII",
    "chr9": "IX",
    "chr10": "X",
    "chr11": "XI",
    "chr12": "XII",
    "chr13": "XIII",
    "chr14": "XIV",
    "chr15": "XV",
    "chr16": "XVI"}

S288C_to_roman : dict[str,str] = {
    'ref|NC_001133|':"I",
    'ref|NC_001134|':"II",
    'ref|NC_001135|':"III",
    'ref|NC_001136|':"IV",
    'ref|NC_001137|':"V",
    'ref|NC_001138|':"VI",
    'ref|NC_001139|':"VII",
    'ref|NC_001140|':"VIII",
    'ref|NC_001141|':"IX",
    'ref|NC_001142|':"X",
    'ref|NC_001143|':"XI",
    'ref|NC_001144|':"XII",
    'ref|NC_001145|':"XIII",
    'ref|NC_001146|':"XIV",
    'ref|NC_001147|':"XV",
    'ref|NC_001148|':"XVI"}

roman_to_S288C: dict[str,str] = {
    "I":    'ref|NC_001133|',
    "II":   'ref|NC_001134|',
    "III":  'ref|NC_001135|',
    "IV":   'ref|NC_001136|',
    "V":    'ref|NC_001137|',
    "VI":   'ref|NC_001138|',
    "VII":  'ref|NC_001139|',
    "VIII": 'ref|NC_001140|',
    "IX":   'ref|NC_001141|',
    "X":    'ref|NC_001142|',
    "XI":   'ref|NC_001143|',
    "XII":  'ref|NC_001144|',
    "XIII": 'ref|NC_001145|',
    "XIV":  'ref|NC_001146|',
    "XV":   'ref|NC_001147|',
    "XVI":  'ref|NC_001148|'}

S288C_to_numeric: dict[str,str] = {
    'ref|NC_001133|':"chr1",
    'ref|NC_001134|':"chr2",
    'ref|NC_001135|':"chr3",
    'ref|NC_001136|':"chr4",
    'ref|NC_001137|':"chr5",
    'ref|NC_001138|':"chr6",
    'ref|NC_001139|':"chr7",
    'ref|NC_001140|':"chr8",
    'ref|NC_001141|':"chr9",
    'ref|NC_001142|':"chr10",
    'ref|NC_001143|':"chr11",
    'ref|NC_001144|':"chr12",
    'ref|NC_001145|':"chr13",
    'ref|NC_001146|':"chr14",
    'ref|NC_001147|':"chr15",
    'ref|NC_001148|':"chr16"}

haploid_samples = [
    'Mal_JK11', 'Mal_JK14', 'Mal_JK16', 'Mal_JK17', 'Mal_JK19', 'Mal_JK21', 'Mal_JK23', 'Mal_JK24', 'Mal_JK26', 'Mal_JK27',
    'Mal_JK29', 'Mal_JK30', 'Mal_JK31', 'Mal_JK35', 'Mal_JK38', 'Mal_JK40', 'Mal_JK41', 'Mal_JK42', 'Mal_JK44', 'Mal_JK45',
    'Mal_JK5-', 'Mal_JK58', 'Mal_JK59', 'Mal_JK6-', 'Mal_JK60', 'Mal_JK62', 'Mal_JK63', 'Mal_JK64', 'Mal_JK66', 'Mal_JK68',
    'Mal_JK69', 'Mal_JK7-', 'Mal_JK70', 'Mal_JK71', 'Mal_JK72', 'Mal_JK73', 'Mal_JK74', 'Mal_JK75', 'Mal_JK79', 'Mal_JK8-',
    'Mal_JK80', 'Mal_JK85', 'Mal_JK87', 'Mal_JK88', 'Mal_JK9-', 'Mal_JK90', 'Mal_JK92', 'Mal_JK93', 'Mal_JK96', 'JT13_CKD',
    'JT14_CKD', 'JT15_CKD', 'JT16_CKD', 'JT17_CKD', 'JT21_CKD', 'JT24_CKD', 'JT25_CKD', 'JT27_CKD', 'JT28_CKD', 'JT31_CKD',
    'JT36_CKD', 'JT37_CKD', 'JT38_CKD', 'JT39_Mer', 'JT44_Mer', 'JT49_CKD', 'JT50_CKD', 'JT53_CKD', 'JT55_CKD', 'JT57_CKD',
    'JT61_CKD', 'JT62_CKD', 'JT63_CKD', 'JMTD_10_', 'JMTD_12_', 'JMTD_13_', 'JMTD_15_', 'JMTD_16_', 'JMTD_17_', 'JMTD_18_',
    'JMTD_19_', 'JMTD_1_C', 'JMTD_21_', 'JMTD_22_', 'JMTD_23_', 'JMTD_2_C', 'JMTD_3_C', 'JMTD_5_C', 'JMTD_6_C', 'JMTD_8_C',
    'JMTD_9_C', 'JK2_10B_', 'JK2_11B_', 'JK2_12B_', 'JK2_13_G', 'JK2_14_A', 'JK2_15_C', 'JK2_19B_', 'JK2_1_CT', 'JK2_21_T',
    'JK2_22_A', 'JK2_23_A', 'JK2_24_A', 'JK2_25_T', 'JK2_26_G', 'JK2_27_A', 'JK2_28_T', 'JK2_29_G', 'JK2_2_TA', 'JK2_30_A',
    'JK2_31_G', 'JK2_32_G', 'JK2_33_A', 'JK2_34_G', 'JK2_35_C', 'JK2_36_A', 'JK2_37_C', 'JK2_38_G', 'JK2_39_G', 'JK2_40_G',
    'JK2_41_T', 'JK2_42_C', 'JK2_43_A', 'JK2_44_G', 'JK2_45_T', 'JK2_47_A', 'JK2_49_C', 'JK2_4_TG', 'JK2_51_T', 'JK2_52_G',
    'JK2_53_T', 'JK2_56_C', 'JK2_57_A', 'JK2_58_A', 'JK2_59_C', 'JK2_5_GA', 'JK2_60_A', 'JK2_61_C', 'JK2_62_C', 'JK2_63_A',
    'JK2_64_G', 'JK2_65_C', 'JK2_68_A', 'JK2_70_G', 'JK2_71_C', 'JK2_72_G', 'JK2_73_G', 'JK2_74_G', 'JK2_75_A', 'JK2_77_C',
    'JK2_78_T', 'JK2_7_GT', 'JK2_91_A', 'JK2_98_A', 'JK2_9B_C']

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
    'Mal_JK35', 
    # 'Mal_JK38', #just empty snps file...
    # 'Mal_JK40', #no clusters
    # 'Mal_JK41', #no clusters
    # 'Mal_JK42', #UNG1
    # 'Mal_JK44', #UNG1
    # 'Mal_JK45', #UNG1
    'Mal_JK5-', 
    'Mal_JK58', 
    'Mal_JK59', 
    'Mal_JK6-', 
    # 'Mal_JK60', #no coverage...
    'Mal_JK62', 
    'Mal_JK63', 
    'Mal_JK64', 
    'Mal_JK66', 
    'Mal_JK68',
    'Mal_JK7-', 
    # 'Mal_JK74', #UNG1
    # 'Mal_JK75', #UNG1
    # 'Mal_JK79', #UNG1
    'Mal_JK8-',
    # 'Mal_JK80', #ung1∆ EV
    # 'Mal_JK85', #UNG1 EV
    # 'Mal_JK87', #UNG1 EV
    # 'Mal_JK88', #UNG1 EV
    'Mal_JK9-', 
    # 'Mal_JK90', #UNG1 EV badquality?
    # 'Mal_JK92', #UNG1 EV
    # 'Mal_JK93', #UNG1 EV missing?
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
    # 'JK2_4_TG', #no clusters

    # 'JK2_21_T', #tetrad 1
    # 'JK2_22_A', #tetrad 1
    # 'JK2_23_A', #tetrad 1
    # 'JK2_24_A', #tetrad 1
    # 'JK2_25_T', #tetrad 2
    # 'JK2_26_G', #tetrad 2
    # 'JK2_27_A', #tetrad 2
    # 'JK2_28_T', #tetrad 2
    # 'JK2_29_G', #tetrad 3
    # 'JK2_30_A', #tetrad 3
    # 'JK2_31_G', #tetrad 3
    # 'JK2_32_G', #tetrad 3
    # 'JK2_33_A', #tetrad 4
    # 'JK2_34_G', #tetrad 4
    # 'JK2_35_C', #tetrad 4
    # 'JK2_36_A', #tetrad 4
    # 'JK2_37_C', #tetrad 5
    # 'JK2_38_G', #tetrad 5
    # 'JK2_39_G', #tetrad 5
    # 'JK2_40_G', #tetrad 5
    # 'JK2_41_T', #tetrad non-selected
    # 'JK2_42_C', #tetrad non-selected
    # 'JK2_43_A', #tetrad non-selected
    # 'JK2_44_G', #tetrad non-selected
    # 'JK2_45_T', #tetrad non-selected
    # 'JK2_47_A', #tetrad non-selected
    # 'JK2_49_C', #tetrad non-selected
    # 'JK2_51_T', #tetrad non-selected
    # 'JK2_52_G', #tetrad non-selected
    # 'JK2_53_T', #tetrad non-selected

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
    'JK2_68_A', #? noisy SNPS, map loos ok
    # 'JK2_70_G', # very noisy SNPs
    'JK2_71_C', #? noisy SNPS, map  looks ok
    'JK2_72_G', 
    'JK2_73_G', 
    'JK2_74_G', 
    'JK2_75_A', 
    'JK2_77_C', 
    'JK2_78_T', 
    'JK2_7_GT', 
    # 'JK2_91_A', #non-selected
    # 'JK2_98_A', #non-selected
    'JK2_9B_C'] #61?

spo13_spo11_samples = [    
    "A50-3.cs",
    "A51-3.cs",
    "A52-3.cs",
    "A53-3.cs",
    "A56-3.cs",
    "A57-3.cs",
    "A58-3.cs",
    "A63-3.cs",
    "A64-3.cs",
    "JK5_10-3",
    "JK5_11-3",
    "JK5_12-3",
    "JK5_1-3.",
    "JK5_13-3",
    "JK5_14-3",
    "JK5_15-3",
    "JK5_16-3",
    "JK5_17-3",
    "JK5_18-3",
    "JK5_19-3",
    "JK5_20-3",
    "JK5_21-3",
    "JK5_22-3",
    "JK5_23-3",
    "JK5_24-A",
    "JK5_25-B",
    "JK5_26-C",
    "JK5_27-D",
    "JK5_28-E",
    "JK5_29-F",
    "JK5_31-G",
    "JK5_32-H",
    "JK5_33-A",
    "JK5_34-B",
    "JK5_35-C",
    "JK5_36-D",
    "JK5_37-E",
    "JK5_38-F",
    "JK5_39-G",
    "JK5_40-H",
    "JK5_51-A",
    "JK5_52-C",
    "JK5_53-D",
    "JK5_57-E",
    "JK5_58-B",
    "JK5_5-F4",
    "JK5_62-C",
    "JK5_66-D",
    "JK5_6-G4",
    "JK5_7-H4",
    "JK5_8-A5",
    "JMT65_CK",
    "JMT66_CK",
    "JMT67_CK",
    "JMT68_CK",
    "JMT69_CK",
    "JMT70_CK",
    "JMT71_CK",
    "JMT72_CK"]

specific_genotypes_4 = [
        [
        'spo13∆spo11∆_sectors_1-2_1-1', #there is 4 of these
        'spo13∆spo11∆_sectors_1-2_1-2', #there is 4 of these
        ],
        [
        'spo13∆_sectors_3-1_2-1', #there is 4 of these
        'spo13∆_sectors_3-1_2-2', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_sectors_1-1_2-2', #there is 4 of these
        'spo13∆spo11∆_sectors_1-1_2-1', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_sectors_1-1_1-1', #there is 4 of these
        'spo13∆spo11∆_sectors_1-1_1-2', #there is 4 of these
        ],
        [
        'spo13∆_sectors_2-1_1-1', #there is 4 of these
        'spo13∆_sectors_2-1_1-2', #there is 4 of these
        ],
        [
        'spo13∆_sectors_3-1_1-1', #there is 4 of these
        'spo13∆_sectors_3-1_1-2', #there is 4 of these
        ],
        [
        'ung1∆ tetrad 1', #there is 4 of these
        'ung1∆ tetrad 2', #there is 4 of these
        'ung1∆ tetrad 3', #there is 4 of these
        'ung1∆ tetrad 4', #there is 4 of these
        'ung1∆ tetrad 5', #there is 4 of these
        'ung1∆ tetrad 6', #there is 4 of these
        'ung1∆ tetrad 7', #there is 4 of these
        'ung1∆ tetrad 8', #there is 4 of these
        ],

        #JK7 SET
        [
        'spo13∆_2sectors_21_1', #there is 4 of these
        'spo13∆_2sectors_21_2', #there is 4 of these
        ],
        [
        'spo13∆_2sectors_20_1', #there is 4 of these
        'spo13∆_2sectors_20_2'  #there is 4 of these
        ],
        [
        'spo13∆_2sectors_70_1', #there is 4 of these
        'spo13∆_2sectors_70_2'  #there is 4 of these
        ],
        [
        'spo13∆_2sectors_12_1', #there is 4 of these
        'spo13∆_2sectors_12_2'  #there is 4 of these
        ],
        [
        'spo13∆spo11∆_2sectors_30_1', #there is 4 of these
        'spo13∆spo11∆_2sectors_30_2'  #there is 4 of these
        ],
        [
        'spo13∆spo11∆_2sectors_40_1', #there is 4 of these
        'spo13∆spo11∆_2sectors_40_2'  #there is 4 of these
        ],
        [
        'spo13∆spo11∆_2sectors_100_1', #there is 4 of these
        'spo13∆spo11∆_2sectors_100_2'  #there is 4 of these
        ],
        [
        'spo13∆spo11∆_2sectors_200_1', #there is 4 of these
        'spo13∆spo11∆_2sectors_200_2'  #there is 4 of these
        ],
        [
        'spo13∆spo11∆_2sectors_50_1', #there is 4 of these
        'spo13∆spo11∆_2sectors_50_2'  #there is 4 of these
        ],
        # [
        # 'spo13∆_2sectors_16_1', #there is 2 of these same as spo13∆_2sectors_160_2
        # ],
        [
        'spo13∆_2sectors_160_1', #there is 2 of these
        'spo13∆_2sectors_160_2'  #there is 4 of these #2 more of these from spo13∆_2sectors_16_1
        ],
        [
        'spo13∆_4sectors_A3_T', #there is 4 of these
        'spo13∆_4sectors_A3_B', #there is 4 of these
        ],
        [
        'spo13∆_4sectors_A4_T', #there is 4 of these
        'spo13∆_4sectors_A4_B', #there is 4 of these
        ],
        [
        'spo13∆_4sectors_A6_T', #there is 4 of these
        'spo13∆_4sectors_A6_B', #there is 4 of these
        ],
        [
        'spo13∆_4sectors_B2_T', #there is 4 of these
        'spo13∆_4sectors_B2_B', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_4sectors_C1_T', #there is 4 of these
        'spo13∆spo11∆_4sectors_C1_B', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_4sectors_C3_T', #there is 4 of these
        'spo13∆spo11∆_4sectors_C3_B', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_4sectors_C4_T', #there is 4 of these
        'spo13∆spo11∆_4sectors_C4_B', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_4sectors_C5_T', #there is 4 of these
        'spo13∆spo11∆_4sectors_C5_B', #there is 4 of these
        ],
        [
        'spo13∆spo11∆_4sectors_C6_T', #there is 4 of these
        'spo13∆spo11∆_4sectors_C6_B', #there is 4 of these
        ],
        [
        'spo13∆_4sectors_D3_T', #there is 4 of these
        'spo13∆_4sectors_D3_B', #there is 4 of these
        ],
        [
        'spo13∆_4sectors_D4_T', #there is 4 of these
        'spo13∆_4sectors_D4_B', #there is 4 of these
        ]
    ]

specific_genotypes_2 = [
        [
        'ung1∆_sectors_5-1_3-1', #there is 2 of these
        'ung1∆_sectors_5-1_3-2', #there is 2 of these
        'ung1∆_sectors_5-1_3-3', #there is 2 of these
        'ung1∆_sectors_5-1_3-4', #there is 2 of these
        ],
        [
        'ung1∆_sectors_5-2_1-1', #there is 2 of these
        'ung1∆_sectors_5-2_1-2', #there is 2 of these
        ],
        [
        'ung1∆_sectors_5-2_2-1', #there is 2 of these
        'ung1∆_sectors_5-2_2-2', #there is 2 of these
        ],
        [
        'ung1∆_sectors_5-1_5-1', #there is 2 of these
        'ung1∆_sectors_5-1_5-2', #there is 2 of these
        'ung1∆_sectors_5-1_5-3'  #there is 2 of these
        ]
    ]

specific_genotypes_3 = [
        [
        'ung1∆_inc_tet_sectors_1-1', #there is 4 of these
        'ung1∆_inc_tet_sectors_1-2', #there is 4 of these
        'ung1∆_inc_tet_sectors_1-3', #there is 4 of these
        ],
        [
        'ung1∆_inc_tet_sectors_2-1', #there is 4 of these
        'ung1∆_inc_tet_sectors_2-2', #there is 4 of these
        'ung1∆_inc_tet_sectors_2-3', #there is 4 of these
        ],
        [
        'ung1∆_inc_tet_sectors_3-1', #there is 4 of these
        'ung1∆_inc_tet_sectors_3-2', #there is 4 of these
        'ung1∆_inc_tet_sectors_3-3', #there is 4 of these
        ],
        [
        'ung1∆_inc_tet_sectors_4-1', #there is 4 of these
        'ung1∆_inc_tet_sectors_4-2', #there is 4 of these
        'ung1∆_inc_tet_sectors_4-3', #there is 4 of these
        ]
    ]

sectors_wt = [

    "JK2_21_T",
    "JK2_22_A",
    "JK2_23_A",
    "JK2_24_A",

    "JK2_25_T",
    "JK2_26_G",
    "JK2_27_A",
    "JK2_28_T",

    "JK2_29_G",
    "JK2_30_A",
    "JK2_31_G",
    "JK2_32_G",

    "JK2_33_A",
    "JK2_34_G",
    "JK2_35_C",
    "JK2_36_A",

    "JK2_37_C",
    "JK2_38_G",
    "JK2_39_G",
    "JK2_40_G",

    "JMT73_CK",
    "JMT74_CK",

    "JMT75_CK",
    "JMT76_CK",

    "JMT77_CK",
    "JMT78_me",

    "JMT79_me",
    "JMT80_CK",
    # ########### #inc. tetrad (2/4)
    "JMT81_me",
    "JMT82_CK",

    "JMT83_CK",
    "JMT84_CK",
    # ########### #inc. tetrad (2/4)
    "JMT85_me",
    "JMT86_CK",

    "JMT87_CK",
    "JMT88_CK",
    # ########### #inc. tetrad (3/4)
    "JMT89_CK",
    "JMT90_me",

    "JMT91_CK",
    "JMT92_CK",

    "JMT93_CK",
    "JMT94_CK",

    "Mal_JK46",     # "tetrad 6", #in-set tetrad 1 (what is genotype for all?)
    "Mal_JK47",     # "tetrad 6", #in-set tetrad 1
    "Mal_JK48",     # "tetrad 6", #in-set tetrad 1
    "Mal_JK49",     # "tetrad 6", #in-set tetrad 1

    "Mal_JK50",     # "tetrad 7", #in-set tetrad 2
    "Mal_JK51",     # "tetrad 7", #in-set tetrad 2
    "Mal_JK52",     # "tetrad 7", #in-set tetrad 2
    "Mal_JK53",     # "tetrad 7", #in-set tetrad 2

    "Mal_JK54",     # "tetrad 8", #in-set tetrad 3
    "Mal_JK55",     # "tetrad 8", #in-set tetrad 3
    "Mal_JK56",     # "tetrad 8", #in-set tetrad 3
    "Mal_JK57",     # "tetrad 8", #in-set tetrad 3

    "TET_1_1_",
    "TET_1_2_",
    "TET_1_3_",
    "TET_1_4_",
    "TET_1_5_",
    "TET_1_6_",
    "TET_1_7_",
    "TET_1_8_",
    "TET_1_9_",
    "TET_1_10",
    "TET_1_11",
    "TET_1_12",

    "TET_2_1_",
    "TET_2_2_",
    "TET_2_3_",
    "TET_2_4_",
    "TET_2_5_",
    "TET_2_6_",
    "TET_2_7_",
    "TET_2_8_",
    "TET_2_9_",
    "TET_2_10",
    "TET_2_11",
    "TET_2_12",

    "TET_3_1_",
    "TET_3_2_",
    "TET_3_3_",
    "TET_3_4_",
    "TET_3_5_",
    "TET_3_6_",
    "TET_3_7_",
    "TET_3_8_",
    "TET_3_9_",
    "TET_3_10",
    "TET_3_11",
    "TET_3_12",

    "TET_4_1_",
    "TET_4_2_",
    "TET_4_3_",
    "TET_4_4_",
    "TET_4_5_",
    "TET_4_6_",
    "TET_4_7_",
    "TET_4_8_",
    "TET_4_9_",
    "TET_4_10",
    "TET_4_11",
    "TET_4_12",
]

sectors_wt_reduced = [

    "JK2_21_T",
    # "JK2_22_A",
    # "JK2_23_A",
    # "JK2_24_A",

    "JK2_25_T",
    # "JK2_26_G",
    # "JK2_27_A",
    # "JK2_28_T",

    "JK2_29_G",
    # "JK2_30_A",
    # "JK2_31_G",
    # "JK2_32_G",

    "JK2_33_A",
    # "JK2_34_G",
    # "JK2_35_C",
    # "JK2_36_A",

    "JK2_37_C",
    # "JK2_38_G",
    # "JK2_39_G",
    # "JK2_40_G",

    "JMT73_CK",
    # "JMT74_CK",

    "JMT75_CK",
    # "JMT76_CK",

    "JMT77_CK",
    # "JMT78_me",

    "JMT79_me",
    # "JMT80_CK",
    # ########### #inc. tetrad (2/4)
    "JMT81_me",
    # "JMT82_CK",

    "JMT83_CK",
    # "JMT84_CK",
    # ########### #inc. tetrad (2/4)
    "JMT85_me",
    # "JMT86_CK",

    "JMT87_CK",
    # "JMT88_CK",
    # ########### #inc. tetrad (3/4)
    "JMT89_CK",
    # "JMT90_me",

    "JMT91_CK",
    # "JMT92_CK",

    "JMT93_CK",
    # "JMT94_CK",

    "Mal_JK46",     # "tetrad 6", #in-set tetrad 1 (what is genotype for all?)
    # "Mal_JK47",     # "tetrad 6", #in-set tetrad 1
    # "Mal_JK48",     # "tetrad 6", #in-set tetrad 1
    # "Mal_JK49",     # "tetrad 6", #in-set tetrad 1

    "Mal_JK50",     # "tetrad 7", #in-set tetrad 2
    # "Mal_JK51",     # "tetrad 7", #in-set tetrad 2
    # "Mal_JK52",     # "tetrad 7", #in-set tetrad 2
    # "Mal_JK53",     # "tetrad 7", #in-set tetrad 2

    "Mal_JK54",     # "tetrad 8", #in-set tetrad 3
    # "Mal_JK55",     # "tetrad 8", #in-set tetrad 3
    # "Mal_JK56",     # "tetrad 8", #in-set tetrad 3
    # "Mal_JK57",     # "tetrad 8", #in-set tetrad 3

    "TET_1_1_",
    # "TET_1_2_",
    # "TET_1_3_",
    # "TET_1_4_",
    "TET_1_5_",
    # "TET_1_6_",
    # "TET_1_7_",
    # "TET_1_8_",
    "TET_1_9_",
    # "TET_1_10",
    # "TET_1_11",
    # "TET_1_12",

    "TET_2_1_",
    # "TET_2_2_",
    # "TET_2_3_",
    # "TET_2_4_",
    "TET_2_5_",
    # "TET_2_6_",
    # "TET_2_7_",
    # "TET_2_8_",
    "TET_2_9_",
    # "TET_2_10",
    # "TET_2_11",
    # "TET_2_12",

    "TET_3_1_",
    # "TET_3_2_",
    # "TET_3_3_",
    # "TET_3_4_",
    "TET_3_5_",
    # "TET_3_6_",
    # "TET_3_7_",
    # "TET_3_8_",
    "TET_3_9_",
    # "TET_3_10",
    # "TET_3_11",
    # "TET_3_12",

    "TET_4_1_",
    # "TET_4_2_",
    # "TET_4_3_",
    # "TET_4_4_",
    "TET_4_5_",
    # "TET_4_6_",
    # "TET_4_7_",
    # "TET_4_8_",
    "TET_4_9_",
    # "TET_4_10",
    # "TET_4_11",
    # "TET_4_12",

]

sectors_spo13 = [
    "JMT9_CKD",
    "JMT10_CK",
    "JMT11_CK",
    "JMT12_CK",

    "JMT13_CK",
    "JMT14_CK",
    "JMT15_CK",
    "JMT16_CK",

    "JMT33_CK",
    "JMT34_CK",
    "JMT35_CK",
    "JMT36_CK",

    "JMT37_CK",
    "JMT38_CK",
    "JMT39_CK",
    "JMT40_CK",

    "JMT41_CK", #RTG!
    "JMT42_CK", #RTG!
    "JMT43_CK", #RTG!
    "JMT44_CK", #RTG!

    "JMT45_CK", #RTG!
    "JMT46_CK", #RTG!
    "JMT47_CK", #RTG!
    "JMT48_CK", #RTG!

    "JK7_1_AT",
    "JK7_2_AT",
    "JK7_3_AT",
    "JK7_4_AT",

    "JK7_5_AT",
    "JK7_6_AT",
    "JK7_7_AT",
    "JK7_8_AT",

    "JK7_11_T",
    "JK7_12_T",
    "JK7_13_T",
    "JK7_14_T",

    "JK7_15_T",
    "JK7_16_T",
    "JK7_17_C",
    "JK7_18_C",

    "JK7_19_C",
    "JK7_20_C",
    "JK7_21_C",
    "JK7_22_C",

    "JK7_23_C",
    "JK7_24_C",
    "JK7_25_G",
    "JK7_26_G",

    "JK7_27_G",
    "JK7_28_G",
    "JK7_29_G",
    "JK7_30_G",

    "JK7_31_G",
    "JK7_32_G",
    "JK7_33_A",
    "JK7_34_A",

    "SD2_CKDN",
    "A3T2_CKD",
    "A3T3_CKD",
    "A3T4_CKD",

    "A3B1_CKD",
    "A3B2_CKD",
    "A3B3_CKD",
    "A3B4_CKD",

    "SD3_CKDN", #RTG!
    "A4T2_CKD", #RTG!
    "A4T3_CKD", #RTG!
    "A4T4_CKD", #RTG!

    "A4B1_CKD", #RTG!
    "A4B2_CKD", #RTG!
    "A4B3_CKD", #RTG!
    "A4B4_CKD", #RTG!

    "SD4_CKDN",
    "A6T2_CKD",
    "A6T3_CKD",
    "A6T4_CKD",

    "A6B1_CKD",
    "A6B2_CKD",
    "A6B3_CKD",
    "A6B4_CKD",

    "SD6_CKDN",
    "B2T2_CKD",
    "B2T3_CKD",
    "B2T4_CKD",

    "B2B1_CKD",
    "B2B2_CKD",
    "B2B3_CKD",
    "B2B4_CKD",

    "SD17_CKD",
    "D3T2_CKD",
    "D3T3_CKD",
    "D3T4_CKD",

    "D3B1_CKD",
    "D3B2_CKD",
    "D3B3_CKD",
    "D3B4_CKD",

    "SD18_Mer",
    "D4T2_CKD",
    "D4T3_CKD",
    "D4T4_CKD",

    "D4B1_CKD",
    "D4B2_CKD",
    "D4B3_CKD",
    "D4B4_CKD"
]

sectors_spo13_reduced = [
    "JMT9_CKD",
    # "JMT10_CK",
    # "JMT11_CK",
    # "JMT12_CK",

    "JMT13_CK",
    # "JMT14_CK",
    # "JMT15_CK",
    # "JMT16_CK",

    "JMT33_CK",
    # "JMT34_CK",
    # "JMT35_CK",
    # "JMT36_CK",

    "JMT37_CK",
    # "JMT38_CK",
    # "JMT39_CK",
    # "JMT40_CK",

    "JMT41_CK",   #RTG!
    # "JMT42_CK", #RTG!
    # "JMT43_CK", #RTG!
    # "JMT44_CK", #RTG!

    "JMT45_CK",   #RTG!
    # "JMT46_CK", #RTG!
    # "JMT47_CK", #RTG!
    # "JMT48_CK", #RTG!

    "JK7_1_AT",
    # "JK7_2_AT",
    # "JK7_3_AT",
    # "JK7_4_AT",

    "JK7_5_AT",
    # "JK7_6_AT",
    # "JK7_7_AT",
    # "JK7_8_AT",

    "JK7_11_T",
    # "JK7_12_T",
    # "JK7_13_T",
    # "JK7_14_T",

    "JK7_15_T",
    # "JK7_16_T",
    # "JK7_17_C",
    # "JK7_18_C",

    "JK7_19_C",
    # "JK7_20_C",
    # "JK7_21_C",
    # "JK7_22_C",

    "JK7_23_C",
    # "JK7_24_C",
    # "JK7_25_G",
    # "JK7_26_G",

    "JK7_27_G",
    # "JK7_28_G",
    # "JK7_29_G",
    # "JK7_30_G",

    "JK7_31_G",
    # "JK7_32_G",
    # "JK7_33_A",
    # "JK7_34_A",

    "SD2_CKDN",
    # "A3T2_CKD",
    # "A3T3_CKD",
    # "A3T4_CKD",

    "A3B1_CKD",
    # "A3B2_CKD",
    # "A3B3_CKD",
    # "A3B4_CKD",

    "SD3_CKDN",   #RTG!
    # "A4T2_CKD", #RTG!
    # "A4T3_CKD", #RTG!
    # "A4T4_CKD", #RTG!

    "A4B1_CKD",   #RTG!
    # "A4B2_CKD", #RTG!
    # "A4B3_CKD", #RTG!
    # "A4B4_CKD", #RTG!

    "SD4_CKDN",
    # "A6T2_CKD",
    # "A6T3_CKD",
    # "A6T4_CKD",

    "A6B1_CKD",
    # "A6B2_CKD",
    # "A6B3_CKD",
    # "A6B4_CKD",

    "SD6_CKDN",
    # "B2T2_CKD",
    # "B2T3_CKD",
    # "B2T4_CKD",

    "B2B1_CKD",
    # "B2B2_CKD",
    # "B2B3_CKD",
    # "B2B4_CKD",

    "SD17_CKD",
    # "D3T2_CKD",
    # "D3T3_CKD",
    # "D3T4_CKD",

    "D3B1_CKD",
    # "D3B2_CKD",
    # "D3B3_CKD",
    # "D3B4_CKD",

    "SD18_Mer",
    # "D4T2_CKD",
    # "D4T3_CKD",
    # "D4T4_CKD",

    "D4B1_CKD",
    # "D4B2_CKD",
    # "D4B3_CKD",
    # "D4B4_CKD"
]

sectors_spo13spo11 = [

    "JMT1_CKD",
    "JMT2_CKD",
    "JMT3_CKD",
    "JMT4_CKD",

    "JMT5_CKD",
    "JMT6_CKD",
    "JMT7_CKD",
    "JMT8_CKD",

    "JMT17_CK",
    "JMT18_CK",
    "JMT19_CK",
    "JMT20_CK",

    "JMT21_CK",
    "JMT22_CK",
    "JMT23_CK",
    "JMT24_CK",

    "JMT25_CK",
    "JMT26_CK",
    "JMT27_CK",
    "JMT28_CK",

    "JMT29_CK",
    "JMT30_CK",
    "JMT31_CK",
    "JMT32_CK",

    "JK7_35_A",
    "JK7_36_A",
    "JK7_37_A",
    "JK7_38_A",

    "JK7_39_A",
    "JK7_40_A",
    "JK7_41_G",
    "JK7_42_G",
    # ###########
    "JK7_43_G",
    "JK7_44_G",
    "JK7_45_G",
    "JK7_46_G",

    "JK7_47_G",
    "JK7_48_G",
    "JK7_49_C",
    "JK7_50_C",
    # ###########
    "JK7_51_C",
    "JK7_52_C",
    "JK7_53_C",
    "JK7_54_C",

    "JK7_55_C",
    "JK7_56_C",
    "JK7_57_T",
    "JK7_58_T",
    # ###########
    "JK7_59_T",
    "JK7_60_T",
    "JK7_61_T",
    "JK7_62_T",

    "JK7_63_T",
    "JK7_64_T",
    "JK7_65_C",
    "JK7_66_C",
    # ###########
    "JK7_67_C",
    "JK7_68_C",
    "JK7_69_C",
    "JK7_70_C",

    "JK7_71_C",
    "JK7_72_C",
    "JK7_73_T",
    "JK7_74_T",
    # ###########
    "JK7_9_TC",
    "JK7_10_T",
    # ###########
    "JK7_75_T",
    "JK7_76_T",

    "JK7_77_T",
    "JK7_78_T",
    "JK7_79_T",
    "JK7_80_T",

    "SD10_Mer",
    "C1T2_CKD",
    "C1T3_CKD",
    "C1T4_CKD",

    "C1B1_mer",
    "C1B2_CKD",
    "C1B3_CKD",
    "C1B4_CKD",
    # ###########
    "SD12_CKD",
    "C3T2_CKD",
    "C3T3_CKD",
    "C3T4_CKD",

    "C3B1_CKD",
    "C3B2_CKD",
    "C3B3_CKD",
    "C3B4_CKD",
    # ###########
    "SD13_CKD",
    "C4T2_CKD",
    "C4T3_CKD",
    "C4T4_CKD",

    "C4B1_CKD",
    "C4B2_CKD",
    "C4B3_CKD",
    "C4B4_CKD",
    # ###########
    "SD14_Mer",
    "C5T2_CKD",
    "C5T3_CKD",
    "C5T4_CKD",

    "C5B1_CKD",
    "C5B2_CKD",
    "C5B3_CKD",
    "C5B4_CKD",
    # ###########
    "SD15_Mer",
    "C6T2_CKD",
    "C6T3_CKD",
    "C6T4_CKD",

    "C6B1_CKD",
    "C6B2_CKD",
    "C6B3_CKD",
    "C6B4_CKD"
]

sectors_spo13spo11_reduced = [

    "JMT1_CKD",
    # "JMT2_CKD",
    # "JMT3_CKD",
    # "JMT4_CKD",

    "JMT5_CKD",
    # "JMT6_CKD",
    # "JMT7_CKD",
    # "JMT8_CKD",

    "JMT17_CK",
    # "JMT18_CK",
    # "JMT19_CK",
    # "JMT20_CK",

    "JMT21_CK",
    # "JMT22_CK",
    # "JMT23_CK",
    # "JMT24_CK",

    "JMT25_CK",
    # "JMT26_CK",
    # "JMT27_CK",
    # "JMT28_CK",

    "JMT29_CK",
    # "JMT30_CK",
    # "JMT31_CK",
    # "JMT32_CK",

    "JK7_35_A",
    # "JK7_36_A",
    # "JK7_37_A",
    # "JK7_38_A",

    "JK7_39_A",
    # "JK7_40_A",
    # "JK7_41_G",
    # "JK7_42_G",
    # ###########
    "JK7_43_G",
    # "JK7_44_G",
    # "JK7_45_G",
    # "JK7_46_G",

    "JK7_47_G",
    # "JK7_48_G",
    # "JK7_49_C",
    # "JK7_50_C",
    # ###########
    "JK7_51_C",
    # "JK7_52_C",
    # "JK7_53_C",
    # "JK7_54_C",

    "JK7_55_C",
    # "JK7_56_C",
    # "JK7_57_T",
    # "JK7_58_T",
    # ###########
    "JK7_59_T",
    # "JK7_60_T",
    # "JK7_61_T",
    # "JK7_62_T",

    "JK7_63_T",
    # "JK7_64_T",
    # "JK7_65_C",
    # "JK7_66_C",
    # ###########
    "JK7_67_C",
    # "JK7_68_C",
    # "JK7_69_C",
    # "JK7_70_C",

    "JK7_71_C",
    # "JK7_72_C",
    # "JK7_73_T",
    # "JK7_74_T",
    # ###########
    "JK7_9_TC",
    # "JK7_10_T",
    # ###########
    "JK7_75_T",
    # "JK7_76_T",

    "JK7_77_T",
    # "JK7_78_T",
    # "JK7_79_T",
    # "JK7_80_T",

    "SD10_Mer",
    # "C1T2_CKD",
    # "C1T3_CKD",
    # "C1T4_CKD",

    "C1B1_mer",
    # "C1B2_CKD",
    # "C1B3_CKD",
    # "C1B4_CKD",
    # ###########
    "SD12_CKD",
    # "C3T2_CKD",
    # "C3T3_CKD",
    # "C3T4_CKD",

    "C3B1_CKD",
    # "C3B2_CKD",
    # "C3B3_CKD",
    # "C3B4_CKD",
    # ###########
    "SD13_CKD",
    # "C4T2_CKD",
    # "C4T3_CKD",
    # "C4T4_CKD",

    "C4B1_CKD",
    # "C4B2_CKD",
    # "C4B3_CKD",
    # "C4B4_CKD",
    # ###########
    "SD14_Mer",
    # "C5T2_CKD",
    # "C5T3_CKD",
    # "C5T4_CKD",

    "C5B1_CKD",
    # "C5B2_CKD",
    # "C5B3_CKD",
    # "C5B4_CKD",
    # ###########
    "SD15_Mer",
    # "C6T2_CKD",
    # "C6T3_CKD",
    # "C6T4_CKD",

    "C6B1_CKD",
    # "C6B2_CKD",
    # "C6B3_CKD",
    # "C6B4_CKD"
]

rtg_list = [
    "JK_92_CA",
    "JK2_48_G",
    "JK2_55_C",
    "JK2_66_T",
    "JK2_97_C",
    "JK2_99_G",
    "JK2_100_",

    "JK3_E3_A",
    "JK3_E11_",

    "JK3_P1_C",
    "JK3_P3_C",
    "JK3_P4_C",
    "JK3_P6_G",
    "JK3_P7_G",
    "JK3_P8_G",
    "JK3_P10_",
    "JK3_P11_",
    "JK3_P13_",
    "JK3_P15_",
    "JK3_P17_",
    "JK3_P18_",
    "JK3_P19_",
    "JK3_P20_",
    "JK3_P21_",
    "JK3_P22_",
    "JK3_P23_",
    "JK3_P24_",
    "JK3_P25_",
    "JK3_P26_",
    "JK3_P27_",
    "JK3_P28_",
    "JK3_P32_",
    "JK3_P33_",
    "JK3_P34_",

    "JK5_41-3",
    "JK5_42-3",
    "JK5_43-3",
    "JK5_44-3",
    "JK5_45-3",

    "JMT41_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT42_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT43_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT44_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT45_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT46_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT47_CK", #this is a dissection of spo13 spores, real rtg?
    "JMT48_CK", #this is a dissection of spo13 spores, real rtg?

    "SD3_CKDN" ,#seems like RTG/non-sporulated
    "A4T2_CKD" ,#seems like RTG/non-sporulated
    "A4T3_CKD" ,#seems like RTG/non-sporulated
    "A4T4_CKD" ,#seems like RTG/non-sporulated

    "A4B1_CKD" ,#seems like RTG/non-sporulated
    "A4B2_CKD" ,#seems like RTG/non-sporulated
    "A4B3_CKD" ,#seems like RTG/non-sporulated
    "A4B4_CKD" ,#seems like RTG/non-sporulated
    
    "JMTD_4_C",

    "JT1_CKDN",
    "JT3_CKDN",
    "JT5_CKDN",
    "JT10_CKD",
    "JT29_CKD",

    "JT40_CKD",
    "JT45_CKD",
    "JT48_CKD",
    "JT52_CKD",
    "JT56_CKD",
    "JT58_CKD",
    "JT59_CKD",

    "Mal_JK15",
    "Mal_JK18",
    "Mal_JK20",
    "Mal_JK34",

    "Sgs-10_R",
    "Sgs-13_R",
    "Sgs-21_R",
    "Sgs-2_R1",

    "SgsEx-26",
    "SgsEx-36",
    "SgsEx-40",
    "SgsEx-5_"

    ]