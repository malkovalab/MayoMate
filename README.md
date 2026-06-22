# MayoMate
[![DOI](https://zenodo.org/badge/927490229.svg)](https://doi.org/10.5281/zenodo.20195874)

`MayoMate` is a Python package designed to facilitate analysis and interpretation of meiotic mutagenesis yeast data. `MayoMate` toolkit allows for discovery and visualization of recombination events, mutations, and mutation clusters in Whole Genome Sequencing (WGS) data of Saccharomyces cerevisiae meiotic outcomes. 

`MayoMate` processes Illumina read alignments and variant calls from CLC Genomics Workbench 20 to identify high-confidence de novo SNVs and mutation clusters. It provides tools to create a Parental Reference Genome for calling parental SNPs using a single modified reference genome. `MayoMate`'s predefined classes store sample information, including SNVs, mutation clusters, parental SNPs, recombination events, sample names, and genotypes. These properties enable event classification, genomic feature correlation, hypermutable ssDNA estimation, mutation rate normalization, statistical analysis, and the creation of customized plots, summaries, and tables.

### Features
- Identifying and visualizing recombination events based on a panel of parental SNPs.
- Identifying, visualizing, and classifying mutation clusters.
- Associating mutation clusters with identified recombination events.
- Inferring the length of hypermutable single-stranded DNA (ssDNA).
- Identifying promoter mutation rates.
- Conducting statistical comparisons between various mutants.

## Requirements

This package was developed and tested on a standard desktop computer with 16GB of RAM using Python 3.10.4. The following packages are required to run `MayoMate`:

- `betterbeeswarm==0.2.0`
- `bio-aid==0.3.1`
- `matplotlib==3.5.2`
- `natsort==8.3.1`
- `numpy==1.22.3`
- `pandas==1.4.2`
- `plotly==5.11.0`
- `scipy==1.8.1`
- `seaborn==0.13.2`

## Installation
To install `MayoMate`, you can clone this repository by running the command below. This step shouldn't take more than a few minutes depending on your internet connection.

```bash
git clone https://github.com/malkovalab/MayoMate.git
cd MayoMate
```

Then, install the required dependencies.

### Installing with conda & pip (recommended)
```bash
conda env create -n mayomate -f environment.yml
conda activate mayomate
```

### Installing with pip alone
```bash
pip install -r requirements.txt
```

## Usage
To use `MayoMate`, please use the `main_playground_clean.ipynb` Jupyter notebook. This notebook provides a step-by-step guide on how to use `MayoMate`'s functions to help analyze and interpret meiotic mutagenesis yeast data. The repository contains a small portion of the dataset used in the full study. You can use that data to test that you have properly configured your environment and that `MayoMate` is working properly. You can see example outputs from the reduced dataset and their compute times within the notebook and within the log files. The outputs created by this notebook can be used for further downstream analysis and visualizations, for example, using scripts below.

In addition to the `main_playground_clean.ipynb` a few other standalone scripts are provided to help with specific tasks and post-processing of the outputs. These scripts are:
- `cluster_analysis.py`: This script is used to conduct analysis on called mutation clusters and compare cluster statistics between different mutants.
- `cluster_coordination.py`: This script is used to identify and visualize the classess of mutation clusters (e.g. "single-switch", "non-switching", "multiple-switch", etc.) and their ssDNA contributions.
- `remake_parental_reference.py`: This script is used to create a Parental Reference Genome for calling parental SNPs with CLC Genomics Workbench.
- `find_imbalanced_snps.py`: This script is used to identify and visualize poorly resolved SNPs within the Parental Reference Genome, enabling their exclusion from the analysis.
- `freqgraph.py`: This script is used to create mutation frequency bar graphs for the ura3-29 reporter.
- `csv_to_paraBed.py`: This script is was used to convert the CLC Genomics Workbench variant call output to a format that can be used with `NIEHS/P-MACD` package. Not used in this version of `MayoMate`.
- `SNP_dens_map.py`: This script is used to create a graph detailing distances between SNPs in the Parental Reference Genome.
- `tracks.ipynb`: This Jupyter notebook is used to reconstruct the ssDNA tracks based on the identified recombination events, mutations, and mutation clusters. It also conducts statistical analysis on the ssDNA tracks between different mutants.
- `transcription_graph`: This script is used to create a graph summarizing sA3A mutation rates for tRNA and protein-coding genes.
- `run_association_simulation_parallel.py`: This script is used to run association simulations between recombination events and mutation clusters on a computing cluster for parallel processing and faster results.

### MayoMate Pipeline Overview
![Meiotic Pipeline_new_bg](https://github.com/user-attachments/assets/62bee6ce-a6d3-4793-94d5-b0d6fe9d154b)


### Data
The data used in the `MayoMate` package is not provided in this repository. The data used in the `MayoMate` package is generated from Whole Genome Sequencing (WGS) data of Saccharomyces cerevisiae meiotic outcomes. The data is first processed using CLC Genomics Workbench 20 and then used as input for the `MayoMate` package to identify recombination events, mutations, and mutation clusters.
Raw reads data can be found in the BioProject and the Sequence Read Archive (SRA) under the accession number PRJNA1225307

### Before you start
Before you start using `MayoMate`, please make sure you have updated the `config.py` file with the correct paths to your input files and directories. The `config.py` file is located in the `MayoMate` directory. 

In addition, please make sure other settings files are also properly configured. These files are within mayo/settings directory.


## Contributing
We welcome contributions to `MayoMate`! If you would like to contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch (git checkout -b feature-branch).
3. Make your changes.
4. Commit your changes (git commit -am 'Add new feature').
5. Push to the branch (git push origin feature-branch).
6. Create a new Pull Request.

## License
MayoMate is licensed under the MIT License.
