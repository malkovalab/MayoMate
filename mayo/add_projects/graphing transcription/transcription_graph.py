import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def prep_for_plot(df: pd.DataFrame) -> pd.DataFrame:

    #take the "Rate * 10^6" and "Control Rate * 10^6" and melt them into a single column called "Rate" and a new column called "Region" that indicates if the rate is the "Rate" or the "Control Rate"
    #df = df.rename(columns={'Rate * 10^6': 'Test', 'Control Rate * 10^6': 'Control'})
    df = df.rename(columns={'Rate per site x10^6': 'Test', 'Control Rate x10^6': 'Control'})
    # if the Spo11 is "spo11" then divide the Test and Control columns by 2
    df.loc[df['Spo11'] == 'spo11', 'Test'] = df.loc[df['Spo11'] == 'spo11', 'Test'] / 2
    df.loc[df['Spo11'] == 'spo11', 'Control'] = df.loc[df['Spo11'] == 'spo11', 'Control'] / 2

    df = df.melt(id_vars=['Spo11', 'Spectra', 'P-Value', 'Feature'], value_vars=['Test', 'Control'], value_name='Rate * 10^6', var_name='Rate_Type')
    
    #convert the P-value row from "1.41×10-16" notation to float
    #df['P-Value'] = df['P-Value'].str.replace('×10', 'E').astype(float)
    df["P_symbol"] = df["P-Value"].apply(lambda x: '****' if x < 0.0001 else '***' if x < 0.001 else '**' if x < 0.01 else '*' if x < 0.05 else 'ns')

    #rename the ['C→T', 'G→A', 'C→T & G→A'] in the "Spectra" column to ["Transcribed", "Non-transcribed", "Both"]
    df['Spectra'] = df['Spectra'].replace({
        'C→T': 'Non-transcribed', 
        'G→A': 'Transcribed', 
        'C→T & G→A': 'Both'})


    df["plot_string"] = df["Spo11"] + " " + df["Spectra"] + " " + df["Feature"] + " " + df["Rate_Type"]
    df["condition"] = df["Spectra"] + " " + df["Rate_Type"]

    order_list_spectr = ['Transcribed', 'Non-transcribed', 'Both']
    order_list_Spo11 = ["SPO11", "spo11∆"]
    order_list_feature = ['TSS', 'tRNA']
    order_list_Rate_Type = ['Control','Test']

    df['Spectra'] = pd.Categorical(df['Spectra'], categories=order_list_spectr, ordered=True)
    df['Spo11'] = pd.Categorical(df['Spo11'], categories=order_list_Spo11, ordered=True)
    df['Feature'] = pd.Categorical(df['Feature'], categories=order_list_feature, ordered=True)
    df['Rate_Type'] = pd.Categorical(df['Rate_Type'], categories=order_list_Rate_Type, ordered=True)

    df = df.sort_values(['Feature','Spectra','Spo11','Rate_Type'])

    return df

def plot_bar(
        df: pd.DataFrame, 
        save_name: str | None = None, 
        show: bool = False,
        figsize: tuple = (16, 8)
        ) -> None:

    sns.set_style("whitegrid")
    sns.set_context("talk")
    sns.font_scale = 2

    #set the figure size
    plt.figure(figsize=figsize)

    #replace spaces with newlines in the genotype column
    df['plot_string'] = df['plot_string'].str.replace(' ', '\n')

    #construct the palette
    #get a list of unique plot_strings
    vars = df['plot_string'].unique()
    #for the combinations of 'Transcribed' and 'Control': light blue
    #for the combinations of 'Transcribed' and 'Test': tab:blue
    #for the combinations of 'Non-transcribed' and 'Control': light orange
    #for the combinations of 'Non-transcribed' and 'Test': tab:orange
    #for the combinations of 'Both' and 'Control': light green
    #for the combinations of 'Both' and 'Test': tab:green
    palette = {}
    for var in vars:
        if 'Both' in var and 'Control' in var:
            palette[var] = '#b2df8a'
        elif 'Both' in var and 'Test' in var:
            palette[var] = 'tab:green'
        elif 'Transcribed' in var and 'Control' in var:
            palette[var] = '#fdbf6f'
        elif 'Transcribed' in var and 'Test' in var:
            palette[var] = 'tab:orange'
        elif 'Non-transcribed' in var and 'Control' in var:
            palette[var] = '#a6cee4'
        elif 'Non-transcribed' in var and 'Test' in var:
            palette[var] = 'tab:blue'

    # for var in vars:
    #     if 'Both' in var and 'Control' in var:
    #         palette[var] = '#b2df8a'
    #     elif 'Both' in var and 'Test' in var:
    #         palette[var] = 'tab:green'
    #     elif 'Transcribed' in var and 'Control' in var:
    #         palette[var] = '#a6cee4'
    #     elif 'Transcribed' in var and 'Test' in var:
    #         palette[var] = 'tab:blue'
    #     elif 'Non-transcribed' in var and 'Control' in var:
    #         palette[var] = '#fdbf6f'
    #     elif 'Non-transcribed' in var and 'Test' in var:
    #         palette[var] = 'tab:orange'
    
    ax = sns.barplot(x='plot_string', y='Rate * 10^6', data=df, palette=palette)
    
    # # Create a dictionary to map the genotype to the index on the x-axis
    plot_string_to_index = {plot_string: index for index, plot_string in enumerate(df['plot_string'].unique())}

    # Add a value label to each error bar
    for i, row in df.iterrows():
        ax.text(plot_string_to_index[row['plot_string']], row['Rate * 10^6'] + 15, f"{row['Rate * 10^6']:.0f}", color='black', ha="center", va="center", fontsize=14, fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', pad=0.20), rotation=0)
    
    index = 0
    for i, row in df.iterrows():
        if index % 2 == 1:
            if row['P_symbol'] == 'ns':
                ax.text((index-0.5), row['Rate * 10^6'] + 40, row['P_symbol'], color='black', ha="center", va="center", fontsize=14, fontweight='bold', rotation=0)
                ax.axhline(y=row['Rate * 10^6']+ 32, color='black', linewidth=1.5, xmin=(index-0.5)/len(df), xmax=(index+0.5)/len(df))
            else:
                ax.text((index-0.5), row['Rate * 10^6'] + 32, row['P_symbol'], color='black', ha="center", va="center", fontsize=24, fontweight='bold', rotation=0)
                ax.axhline(y=row['Rate * 10^6']+ 32, color='black', linewidth=1.5, xmin=(index-0.5)/len(df), xmax=(index+0.5)/len(df))
        index += 1

    plt.grid(True)
    ax.xaxis.grid(False)
    #add one vertical line between tRNA and TSS
    ax.axvline(x=11.5, color='lightgray', linewidth=2)

    plt.title('A3A mutation frequencies around transcription units', fontsize=24, fontweight='bold')
    plt.ylabel('Normalized A3A mutation frequency × $10^{6}$', fontsize=14, fontweight='bold')
    plt.yticks(fontsize=14, fontweight='bold')
    #plt.xlabel('Condition', fontsize=14, fontweight='bold')
    # plt.xticks(fontsize=8, fontweight='bold', rotation=0, fontstyle='italic')
    # do not show x-ticks
    plt.xticks([])
    plt.xlabel('')
    plt.ylim(0, 500)
    plt.margins(x=0)

    #add a labels for the x-axis
    for i in range(0, 24, 4):
        plt.text(i+0.5, -20, 'SPO11', fontsize=12, style='italic', fontweight='bold',ha='center')
        plt.text(i+0.5+2, -20, 'spo11∆', fontsize=12, style='italic', fontweight='bold',ha='center')

    categories = ['Transcribed', 'Non-transcribed', 'Both']
    for i, category in enumerate(categories):
        plt.text(1.5 + i * 4, -40, category, fontsize=12, fontweight='bold', ha='center')
        plt.text(13.5 + i * 4, -40, category, fontsize=12, fontweight='bold', ha='center')

    plt.text(5.5, -60, 'Protein Coding', fontsize=14, fontweight='bold', ha='center')
    plt.text(17.5, -60, 'tRNA', fontsize=14, fontweight='bold', ha='center')

    lines_major = plt.vlines(x=[-0.5, 11.5, 23.5], ymin=-70, ymax=0, color='lightgray', linewidth=2)
    lines_major.set_clip_on(False)

    lines_medium = plt.vlines(x=[3.5, 7.5, 15.5, 19.5], ymin=-45, ymax=0, color='lightgray', linewidth=1.5)
    lines_medium.set_clip_on(False)

    lnes_minor = plt.vlines(x=[1.5, 5.5, 9.5, 13.5, 17.5, 21.5], ymin=-25, ymax=0, color='lightgray', linewidth=1)
    lnes_minor.set_clip_on(False)

    vars = df['condition'].unique()
    palette = {}
    for var in vars:
        if 'Both' in var and 'Control' in var:
            palette[var] = '#b2df8a'
        elif 'Both' in var and 'Test' in var:
            palette[var] = 'tab:green'
        elif 'Transcribed' in var and 'Control' in var:
            palette[var] = '#fdbf6f'
        elif 'Transcribed' in var and 'Test' in var:
            palette[var] = 'tab:orange'
        elif 'Non-transcribed' in var and 'Control' in var:
            palette[var] = '#a6cee4'
        elif 'Non-transcribed' in var and 'Test' in var:
            palette[var] = 'tab:blue'

        # if 'Both' in var and 'Control' in var:
        #     palette[var] = '#b2df8a'
        # elif 'Both' in var and 'Test' in var:
        #     palette[var] = 'tab:green'
        # elif 'Transcribed' in var and 'Control' in var:
        #     palette[var] = '#a6cee4'
        # elif 'Transcribed' in var and 'Test' in var:
        #     palette[var] = 'tab:blue'
        # elif 'Non-transcribed' in var and 'Control' in var:
        #     palette[var] = '#fdbf6f'
        # elif 'Non-transcribed' in var and 'Test' in var:
        #     palette[var] = 'tab:orange'

    handles = [plt.Rectangle((0,0),1,1, color=palette[var]) for var in vars]
    labels = vars
    plt.legend(handles, labels, loc='upper right', fontsize=14)
    plt.tight_layout()


    if save_name is not None:
        plt.savefig(save_name, dpi=300)

    if show:
        plt.show()

    plt.close()

if __name__ == '__main__':

    PATH='newmutrates.txt'
    df = pd.read_table(PATH)
    print(df)

    df = prep_for_plot(df)
    print(df)

    plot_bar(df, save_name='plot.png', show=False)
