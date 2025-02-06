import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mayo import *

#read df_plot_list 
df_plot_list = load_payload_from_pickle("outputs/df_SNPs_mutations_window_all.pkl")

for df_plot in df_plot_list:

    df_plot["SNP_count"] = df_plot["SNP_count"].astype(int)
    df_plot["Mutation_count"] = df_plot["Mutation_count"].astype(int)

    #remove windows with 0 SNPs and 0 mutations
    #df_plot = df_plot[(df_plot["SNP_count"] > 0) & (df_plot["Mutation_count"] > 0)]

    #remove windows with 0 SNPs
    df_plot = df_plot[(df_plot["SNP_count"] > 0)]

    #group by SNP_count and take average of Mutation_count and Rec_event_count
    df_plot = df_plot.groupby("SNP_count").agg({"Mutation_count": "mean", "Rec_event_count": ["mean", "count"]}).reset_index()
    df_plot.columns = ["SNP_count", "Mutation_count", "Rec_event_count", "count"]

    count_total = df_plot["count"].sum()
    df_plot["cumulative_count"] = df_plot["count"].cumsum()
    df_plot = df_plot[df_plot["cumulative_count"] <= count_total * 0.97]

    feature_to_correlate = "Mutation_count" # "Mutation_count" or "Rec_event_count"
    # # use linear regression, get the r2 and p-value

    from scipy.stats import linregress
    slope, intercept, r_value, p_value, std_err = linregress(df_plot["SNP_count"], df_plot[feature_to_correlate])
    print(f"R2={r_value**2}, R={r_value}, p={p_value}")
    #find the correlation between the number of SNPs and mutations in a window, use spearman correlation and pearson correlation
    print(df_plot.corr(method='spearman'))
    print(df_plot.corr(method='pearson'))

    plt.figure(figsize=(6, 6))
    sns.regplot(data=df_plot, x="SNP_count", y=feature_to_correlate, scatter_kws={"s": 30})
    #set origin to 0,0
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.show()

    # create a 3d scatterplot of SNPs, mutations and recombination events
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #color by count
    ax.scatter(df_plot["SNP_count"], df_plot["Mutation_count"], df_plot["Rec_event_count"], c=df_plot["count"], cmap='viridis')
    #make sure the origin is at 0,0,0
    # ax.set_xlim(left=0)
    # ax.set_ylim(bottom=0)
    # ax.set_zlim(bottom=0)

    #regress on SNP_count
    slope_a, intercept_a, r_value, p_value, std_err = linregress(df_plot["SNP_count"], df_plot["Rec_event_count"])
    slope_b, intercept_b, r_value, p_value, std_err = linregress(df_plot["SNP_count"], df_plot["Mutation_count"])
    
    import numpy as np
    x = np.linspace(0, 16, 100)
    y_a = slope_a * x + intercept_a
    y_b = slope_b * x + intercept_b
    ax.plot(x, y_b, y_a, color="red")

    #label the axes
    ax.set_xlabel("SNP count")
    ax.set_ylabel("Mutation count")
    ax.set_zlabel("Recombination event count")

    plt.show()


    #divide the data into bins based on SNP_count and plot boxplot of Mutation_count and Rec_event_count
    # df_plot["SNP_bin"] = pd.cut(df_plot["SNP_count"], bins=range(0, 17, 4), right=False)
    # df_plot["SNP_bin"] = df_plot["SNP_bin"].astype(str)

    # #plot boxplot
    # plt.figure(figsize=(10, 6))
    # #dont show outliers, show the mean
    # sns.boxplot(data=df_plot, x="SNP_bin", y="Mutation_count", showfliers=False, showmeans=True)
    # #set origin to 0,0
    # plt.ylim(bottom=0)
    # plt.xticks(rotation=45)
    # plt.show()

    # #plot boxplot
    # plt.figure(figsize=(10, 6))
    # #dont show outliers, show the mean
    # sns.boxplot(data=df_plot, x="SNP_bin", y="Rec_event_count", showfliers=False, showmeans=True)
    # #set origin to 0,0
    # plt.ylim(bottom=0)
    # plt.xticks(rotation=45)
    # plt.show()
    
