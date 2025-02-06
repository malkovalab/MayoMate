# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

# This script processes the data from the finalfreqs.csv file and plots a bar graph of the meiotic 
# ura3-29 reversion frequencies for each genotype and vector combination. The script also calculates
# the 95% confidence interval for each genotype and vector combination and performs statistical tests
# to compare the premeiotic and meiotic frequencies for each genotype and vector combination, the 
# meiotic frequencies between all genotypes, and the vectors within each genotype.

# PARAMETERS
FILE_PATH = "finalfreqs.csv"
SAVE_PATH = "finalfreqs_processed_mutants.csv"
FIG_SIZE = (10, 8)

GENOTYPES_INCLUDE = [
    "UNG1",
    "ung1∆",
    "ung1∆spo13∆",
    "ung1∆spo11∆spo13∆",
    "ung1∆exo1-nd",
    "ung1∆pol32∆",
    "ung1∆pol32∆exo1-nd",
    "sgs1∆C",
    "sgs1∆Cung1∆",
    ]

# END PARAMETERS

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as stats

def calculate_nth_ci(df: pd.DataFrame, genotype: str, vector: str, n: int) -> pd.DataFrame:

    #get the highest and lowest values for the genotype and vector combination
    high_pre = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "premeiotic*10**10"].nlargest(n).min()
    low_pre = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "premeiotic*10**10"].nsmallest(n).max()
    high_mei = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "meiotic*10**10"].nlargest(n).min()
    low_mei = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "meiotic*10**10"].nsmallest(n).max()

    #calculate the confidence interval
    ci_pre = f"{low_pre}-{high_pre}"
    ci_mei = f"{low_mei}-{high_mei}"

    #add the confidence interval to the dataframe
    df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "premeiotic_ci"] = ci_pre
    df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "meiotic_ci"] = ci_mei

    #print(f"genotype: {genotype}, vector: {vector}, ci_pre: {ci_pre}, ci_mei: {ci_mei}, n: {n}, high_pre: {high_pre}, low_pre: {low_pre}, high_mei: {high_mei}, low_mei: {low_mei}")

    return df

def process_data(df: pd.DataFrame) -> pd.DataFrame:

    #remove the "data_origin" column
    df = df.drop(columns=['data_origin'])

    # #get median of the premeiotic*10**10 and meiotic column for each genotype and vector and add it as a new column
    df['premeiotic_median'] = df.groupby(['genotype', 'vector'])['premeiotic*10**10'].transform('median')
    df['meiotic_median'] = df.groupby(['genotype', 'vector'])['meiotic*10**10'].transform('median')

    df['premeiotic_mean'] = df.groupby(['genotype', 'vector'])['premeiotic*10**10'].transform('mean')
    df['meiotic_mean'] = df.groupby(['genotype', 'vector'])['meiotic*10**10'].transform('mean')
    df['premeiotic_mean'] = df['premeiotic_mean'].round(2)
    df['meiotic_mean'] = df['meiotic_mean'].round(2)

    #get the count of each genotype and vector combination and add it as a new column
    df['count'] = df.groupby(['genotype', 'vector'])['premeiotic*10**10'].transform('count')

    #calculate the fold change between meiotic and premeiotic frequencies (medians) and add it as a new column
    df['fold_change'] = df['meiotic_median'] / df['premeiotic_median']

    df["fold_change"] = df["fold_change"].round(2)

    # calculate the confidence interval for each genotype and vector combination and add it as a new column:
    df["premeiotic_ci"] = None
    df["meiotic_ci"] = None

    for genotype in df["genotype"].unique():
        for vector in df["vector"].unique():

            #get the count for the genotype and vector combination
            count = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "count"].max()
            
            #if count is nan (no data for that genotype and vector combination), skip it
            if np.isnan(count):
                continue

            if count < 9:
                df = calculate_nth_ci(df, genotype, vector, 1)
            
            elif count < 12:
                df = calculate_nth_ci(df, genotype, vector, 2)
            
            elif count < 15:
                df = calculate_nth_ci(df, genotype, vector, 3)

            elif count < 17:
                df = calculate_nth_ci(df, genotype, vector, 4)
            
            elif count < 20:
                df = calculate_nth_ci(df, genotype, vector, 5)

            elif count < 23:
                df = calculate_nth_ci(df, genotype, vector, 6)
            
            elif count < 25:
                df = calculate_nth_ci(df, genotype, vector, 7)
            
            elif count < 28:
                df = calculate_nth_ci(df, genotype, vector, 8)

            elif count < 30:
                df = calculate_nth_ci(df, genotype, vector, 9)

            elif count < 33:
                df = calculate_nth_ci(df, genotype, vector, 10)

            else:
                raise ValueError(f"Many samples for genotype {genotype} and vector {vector}. Please modify the code to calculate the confidence interval for this case... Refer to non-parametric ranks for 95% confidence interval.")

    return df

def run_stats(df: pd.DataFrame) -> None:
    #first compare premeiotic and meiotic frequencies for each genotype and vector combination
    for genotype in df["genotype"].unique():
        for vector in df["vector"].unique():
            #get the premeiotic and meiotic frequencies for the genotype and vector combination
            premeiotic = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "premeiotic*10**10"]
            meiotic = df.loc[(df["genotype"] == genotype) & (df["vector"] == vector), "meiotic*10**10"]

            #if there are no premeiotic or meiotic frequencies for the genotype and vector combination, skip it
            if premeiotic.empty or meiotic.empty:
                continue

            #run mann whitney u test
            u, p = stats.mannwhitneyu(premeiotic, meiotic, alternative='two-sided')
            print(f"premeiotic vs meiotic for genotype {genotype} {vector}: u = {u}, p = {p}")

    #then compare meiotic frequencies between all genotypes
    for vector in df["vector"].unique():
        for genotype1 in df["genotype"].unique():
            for genotype2 in df["genotype"].unique():

                df_genotype1 = df.loc[(df["genotype"] == genotype1) & (df["vector"] == vector), "meiotic*10**10"]
                df_genotype2 = df.loc[(df["genotype"] == genotype2) & (df["vector"] == vector), "meiotic*10**10"]

                if genotype1 == genotype2:
                    continue

                #if there are no meiotic frequencies for the genotype and vector combination, skip it
                if df_genotype1.empty or df_genotype2.empty:
                    continue

                #run mann whitney u test
                u, p = stats.mannwhitneyu(df_genotype1, df_genotype2, alternative='two-sided')
                print(f"meiotic for {vector} {genotype1} vs {genotype2}: u = {u}, p = {p}")

    # lastly compare vectors within each genotype

    for genotype in df["genotype"].unique():

        for timepoint in ["meiotic*10**10", "premeiotic*10**10"]:
        
            df_vector1 = df.loc[(df["genotype"] == genotype) & (df["vector"] == "EV"), timepoint]
            df_vector2 = df.loc[(df["genotype"] == genotype) & (df["vector"] == "sA3A"), timepoint]

            if df_vector1.empty or df_vector2.empty:
                continue

            u, p = stats.mannwhitneyu(df_vector1, df_vector2, alternative='two-sided')
            print(f"{timepoint} for {genotype} EV vs sA3A: u = {u}, p = {p}")

def prep_for_plot(df: pd.DataFrame) -> pd.DataFrame:

    #remove the premeiotic*10**10 and meiotic*10**10 columns
    df = df.drop(columns=['premeiotic_frequency', 'meiotic_frequency', 'premeiotic*10**10', 'meiotic*10**10'])

    #remove duplicates
    df = df.drop_duplicates()

    #reset the index
    df = df.reset_index(drop=True)

    #melt premeiotic_median and meiotic_median into timepoint (premeiotic or meiotic) and frequency, and premeiotic_ci and meiotic_ci into timepoint (premeiotic or meiotic) and ci
    df = df.melt(id_vars=['genotype', 'vector', 'count', 'fold_change'], value_vars=['premeiotic_median', 'meiotic_median', 'premeiotic_mean', 'meiotic_mean' ,'premeiotic_ci', 'meiotic_ci'], var_name='timepoint', value_name='value')

    df_timepoint = df.loc[df['timepoint'].isin(['premeiotic_median', 'meiotic_median'])]
    df_mean = df.loc[df['timepoint'].isin(['premeiotic_mean', 'meiotic_mean'])]
    df_ci = df.loc[df['timepoint'].isin(['premeiotic_ci', 'meiotic_ci'])]

    df_timepoint = df_timepoint.rename(columns={'value': 'median_frequency'})
    df_mean = df_mean.rename(columns={'value': 'mean_frequency'})
    df_ci = df_ci.rename(columns={'value': 'ci'})

    df_timepoint['timepoint'] = df_timepoint['timepoint'].str.replace('_median', '')
    df_mean['timepoint'] = df_mean['timepoint'].str.replace('_mean', '')
    df_ci['timepoint'] = df_ci['timepoint'].str.replace('_ci', '')

    df = pd.merge(df_timepoint, df_mean, on=['genotype', 'vector', 'count', 'timepoint', 'fold_change'])
    df = pd.merge(df, df_ci, on=['genotype', 'vector', 'count', 'timepoint', 'fold_change'])
    
    #sort the dataframe, firstly by vector, genotype, and timepoint (timepoint in reverse order)
    df = df.sort_values(by=['vector', 'genotype', 'timepoint'], ascending=[True, True, False])

    #add the vector to the genotype name
    df['genotype'] =  df['vector'] + ' ' + df['genotype']
    df['genotype'] = df['genotype'] + ' ' + df['timepoint']

    df[["ci_low", "ci_high"]] = df['ci'].str.split('-', expand=True)

    df['ci_low'] = df['ci_low'].astype(float)
    df['ci_high'] = df['ci_high'].astype(float)

    return df

def plot_bar(
        df: pd.DataFrame, 
        save_name: str | None = None, 
        show: bool = False,
        figsize: tuple = (12, 12)
        ) -> None:
    """Plot a bar graph of the frequency over time data.

    Args:
        df: A pandas dataframe containing the frequency over time data.
    """

    #if the genotype includes ung1∆exo1-nd, replace it with ung1∆\nexo1-nd
    df.loc[df["genotype"].str.contains("ung1∆exo1-nd"), "genotype"] = df["genotype"].str.replace('ung1∆exo1-nd', 'ung1∆\nexo1-nd')

    #if the genotype includes ung1∆pol32∆, replace it with ung1∆\npol32∆
    df.loc[df["genotype"].str.contains("ung1∆pol32∆"), "genotype"] = df["genotype"].str.replace('ung1∆pol32∆', 'ung1∆\npol32∆')

    #if the genotype contains "ung1∆spo" add a '\n' after "ung1∆"
    df.loc[df["genotype"].str.contains("ung1∆spo"), "genotype"] = df["genotype"].str.replace('ung1∆spo', 'ung1∆\nspo')

    #if the genotype coontains pol32∆exo1-nd, replace it with pol32∆\nexo1-nd
    df.loc[df["genotype"].str.contains("pol32∆exo1-nd"), "genotype"] = df["genotype"].str.replace('pol32∆exo1-nd', 'pol32∆\nexo1-nd')

    sns.set_style("whitegrid")
    sns.set_context("talk")
    sns.font_scale = 2

    #set the figure size
    plt.figure(figsize=figsize)

    #replace spaces with newlines in the genotype column
    df['genotype'] = df['genotype'].str.replace(' ', '\n')

    #construct the palette
    #get a list of unique genotypes
    genotypes = df['genotype'].unique()
    #for the combinations of EV and premeiotic genotypes: light blue
    #for the combinations of EV and meiotic genotypes: light red
    #for the combinations of sA3A and premeiotic genotypes: dark blue
    #for the combinations of sA3A and meiotic genotypes: dark red
    palette = {}
    for genotype in genotypes:
        if genotype.startswith('EV'):
            if genotype.endswith('premeiotic'):
                palette[genotype] = '#62B0E8' #light blue
            elif genotype.endswith('meiotic'):
                palette[genotype] = '#F07577' #light red
        elif genotype.startswith('sA3A'):
            if genotype.endswith('premeiotic'):
                palette[genotype] = '#1F77B4' #dark blue
            elif genotype.endswith('meiotic'):
                palette[genotype] = '#E31A1C' #dark red

    #plot the data as a bar graph and error bars from the ci_low to ci_high

    #first, order df by the order list
    #order = ['sA3A\nung1∆\nspo13∆\npremeiotic', 'sA3A\nung1∆\nspo13∆\nmeiotic', 'sA3A\nung1∆\nspo11∆spo13∆\npremeiotic', 'sA3A\nung1∆\nspo11∆spo13∆\nmeiotic']
    #df = df.set_index('genotype').loc[order].reset_index()
    
    ax = sns.barplot(x='genotype', y='median_frequency', data=df, errorbar=None, palette=palette) #, order=order

    #plot the error bars
    #ax.errorbar(x=np.arange(len(df)), y=df['median_frequency'], yerr=[df['median_frequency'] - df['ci_low'], df['ci_high'] - df['median_frequency']], fmt='none', ecolor='black', capsize=5)

    
    # Create a dictionary mapping 'genotype' to indices
    genotype_to_index = {genotype: index for index, genotype in enumerate(df['genotype'].unique())}

    # Add a value label to each error bar
    for i, row in df.iterrows():
        # Use the 'genotype' value to get the correct index for the label
        ax.text(genotype_to_index[row['genotype']], row['median_frequency'] - (row['median_frequency'] * 0.5), f"{row['median_frequency']:.0f}", color='black', ha="center", va="center", fontsize=14, fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', pad=0.2), rotation=0)

    #add a range label above the height of the bar (at + 50% of the height of the bar)
    for i, row in df.iterrows():
        # Use the 'genotype' value to get the correct index for the label
        ax.text(genotype_to_index[row['genotype']], row['median_frequency'] + (row['median_frequency'] * 0.1), f"({row['ci'].replace('-', ' - ')})", color='black', ha="center", va="center", fontsize=8, fontweight='bold', rotation=0)

    #set axis to log scale but still show 0
    ax.set(yscale="log")
    ax.set_ylim(bottom=1)

    plt.grid(True)
    #remove vertical grid lines
    ax.xaxis.grid(False)

    plt.title('ura3-29 reversion frequencies', fontsize=24, fontweight='bold')
    plt.ylabel('ura3-29 reversion frequency ($10^{10}$)', fontsize=18, fontweight='bold')
    plt.xticks(fontsize=14, fontweight='bold', rotation=0, fontstyle='italic')
    plt.yticks(fontsize=16, fontweight='bold')
    plt.xlabel('Genotype', fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_name is not None:
        plt.savefig(save_name, dpi=300)

    if show:
        plt.show()

    plt.close()

if __name__ == '__main__':
    #load the data
    df = pd.read_csv(FILE_PATH, sep=',')

    print(df)
    print(df["genotype"].unique())

    #if the genotype contains "spo at the beginning, prepend it with "ung1∆"
    df.loc[df["genotype"].str.startswith("spo1"), "genotype"] = "ung1∆" + df["genotype"]

    df = df.loc[df["genotype"].isin(GENOTYPES_INCLUDE)]

    df = process_data(df)

    print(df["genotype"].value_counts())

    run_stats(df)

    df = prep_for_plot(df)

    df.to_csv(SAVE_PATH, index=False)

    #df = df[df["vector"] != "EV"]

    plot_bar(df, save_name='freq_system.png', show=True, figsize=FIG_SIZE)