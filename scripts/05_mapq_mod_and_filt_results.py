#  __           _       _      ___  ____               _      ___  ____                       _                   _    __ _ _ _                         _ _       
# / _\ ___ _ __(_)_ __ | |_   / _ \| ___|_    /\/\    /_\    / _ \/___ \  _ __ ___   ___   __| |   __ _ _ __   __| |  / _(_) | |_   _ __ ___  ___ _   _| | |_ ___ 
# \ \ / __| '__| | '_ \| __| | | | |___ (_)  /    \  //_\\  / /_)//  / / | '_ ` _ \ / _ \ / _` |  / _` | '_ \ / _` | | |_| | | __| | '__/ _ \/ __| | | | | __/ __|
# _\ \ (__| |  | | |_) | |_  | |_| |___) |  / /\/\ \/  _  \/ ___/ \_/ /  | | | | | | (_) | (_| | | (_| | | | | (_| | |  _| | | |_  | | |  __/\__ \ |_| | | |_\__ \
# \__/\___|_|  |_| .__/ \__|  \___/|____(_) \/    \/\_/ \_/\/   \___,_\  |_| |_| |_|\___/ \__,_|  \__,_|_| |_|\__,_| |_| |_|_|\__| |_|  \___||___/\__,_|_|\__|___/
#                |_|                                                                                                                                              
                 
# Script: 05_mapq_mod_and_filt_results.sh
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: Once obtained the results of the variant calling from the modified and filtered SAM files, it is time to analyze the results. For this, each 
# dataset is loaded into a dataframe, and then they will be outer-merged with indicator set to true. By doing this, we can compare how much each dataset 
# differs from the original one, by defining true positives (common variants), false positives (new variants) and false negatives (missed variants). Finally 
# the results can be visualized as stacked bar plots where the total number and their relative percentages is shown.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read datasets.
base_bwa = pd.read_csv("./data/tsv/bwa/masterbwa.tsv", sep = "\t").iloc[:, :5]
base_gem = pd.read_csv("./data/tsv/bwa/mastergem.tsv", sep = "\t").iloc[:, :5]
bwa_20_mod = pd.read_csv("./data/modified/bwa/tsv/masterbwa_20.tsv", sep = "\t").iloc[:, :5]
bwa_40_mod = pd.read_csv("./data/modified/bwa/tsv/masterbwa_40.tsv", sep = "\t").iloc[:, :5]
bwa_60_mod = pd.read_csv("./data/modified/bwa/tsv/masterbwa_60.tsv", sep = "\t").iloc[:, :5]
bwa_20_filt = pd.read_csv("./data/filtered/bwa/tsv/masterbwa_20.tsv", sep = "\t").iloc[:, :5]
bwa_40_filt = pd.read_csv("./data/filtered/bwa/tsv/masterbwa_40.tsv", sep = "\t").iloc[:, :5]
bwa_60_filt = pd.read_csv("./data/filtered/bwa/tsv/masterbwa_60.tsv", sep = "\t").iloc[:, :5]
gem_20_mod = pd.read_csv("./data/modified/gem/tsv/mastergem_20.tsv", sep = "\t").iloc[:, :5]
gem_40_mod = pd.read_csv("./data/modified/gem/tsv/mastergem_40.tsv", sep = "\t").iloc[:, :5]
gem_60_mod = pd.read_csv("./data/modified/gem/tsv/mastergem_60.tsv", sep = "\t").iloc[:, :5]
gem_20_filt = pd.read_csv("./data/filtered/gem/tsvmastergem_20.tsv", sep = "\t").iloc[:, :5]
gem_40_filt = pd.read_csv("./data/filtered/gem/tsvmastergem_40.tsv", sep = "\t").iloc[:, :5]
gem_60_filt = pd.read_csv("./data/filtered/gem/tsvmastergem_60.tsv", sep = "\t").iloc[:, :5]

# Merge datasets.
merged_20_bwa_mod = pd.merge(base_bwa, bwa_20_mod, how = "outer", indicator = True)
merged_40_bwa_mod = pd.merge(base_bwa, bwa_40_mod, how = "outer", indicator = True)
merged_60_bwa_mod = pd.merge(base_bwa, bwa_60_mod, how = "outer", indicator = True)
merged_20_bwa_filt = pd.merge(base_bwa, bwa_20_filt, how = "outer", indicator = True)
merged_40_bwa_filt = pd.merge(base_bwa, bwa_40_filt, how = "outer", indicator = True)
merged_60_bwa_filt = pd.merge(base_bwa, bwa_60_filt, how = "outer", indicator = True)
merged_20_gem_mod = pd.merge(base_gem, gem_20_mod, how = "outer", indicator = True)
merged_40_gem_mod = pd.merge(base_gem, gem_40_mod, how = "outer", indicator = True)
merged_60_gem_mod = pd.merge(base_gem, gem_60_mod, how = "outer", indicator = True)
merged_20_gem_filt = pd.merge(base_gem, gem_20_filt, how = "outer", indicator = True)
merged_40_gem_filt = pd.merge(base_gem, gem_40_filt, how = "outer", indicator = True)
merged_60_gem_filt = pd.merge(base_gem, gem_60_filt, how = "outer", indicator = True)

# Change the indicator column names to FP, FN, TP.

def change_names(value):
    if value == "right_only":
        return "FP = novel"
    elif value == "left_only":
        return "FN = missed"
    elif value == "both":
        return "TP = common"

merged_20_bwa_mod["_merge"] = merged_20_bwa_mod["_merge"].apply(change_names)
merged_40_bwa_mod["_merge"] = merged_40_bwa_mod["_merge"].apply(change_names)
merged_60_bwa_mod["_merge"] = merged_60_bwa_mod["_merge"].apply(change_names)
merged_20_bwa_filt["_merge"] = merged_20_bwa_filt["_merge"].apply(change_names)
merged_40_bwa_filt["_merge"] = merged_40_bwa_filt["_merge"].apply(change_names)
merged_60_bwa_filt["_merge"] = merged_60_bwa_filt["_merge"].apply(change_names)
merged_20_gem_mod["_merge"] = merged_20_gem_mod["_merge"].apply(change_names)
merged_40_gem_mod["_merge"] = merged_40_gem_mod["_merge"].apply(change_names)
merged_60_gem_mod["_merge"] = merged_60_gem_mod["_merge"].apply(change_names)
merged_20_gem_filt["_merge"] = merged_20_gem_filt["_merge"].apply(change_names)
merged_40_gem_filt["_merge"] = merged_40_gem_filt["_merge"].apply(change_names)
merged_60_gem_filt["_merge"] = merged_60_gem_filt["_merge"].apply(change_names)


# Generate a list of the different dfs so the count and percentage of FN, FP and TP can be done.
list_df_bwa = [
merged_20_bwa_mod,
merged_40_bwa_mod,
merged_60_bwa_mod,
merged_20_bwa_filt,
merged_40_bwa_filt,
merged_60_bwa_filt,
]

list_df_gem = [
merged_20_gem_mod,
merged_40_gem_mod,
merged_60_gem_mod,
merged_20_gem_filt,
merged_40_gem_filt,
merged_60_gem_filt,
]

bars1_bwa = [item["_merge"].value_counts()["TP = common"] for item in list_df_bwa]
bars2_bwa = [item["_merge"].value_counts()["FN = missed"] for item in list_df_bwa]
bars3_bwa = [item["_merge"].value_counts()["FP = novel"] for item in list_df_bwa]
bars1_gem = [item["_merge"].value_counts()["TP = common"] for item in list_df_gem]
bars2_gem = [item["_merge"].value_counts()["FN = missed"] for item in list_df_gem]
bars3_gem = [item["_merge"].value_counts()["FP = novel"] for item in list_df_gem]

#Generate dataframes with the counts and percentages.
total_bars_bwa = pd.DataFrame([bars1_bwa, bars2_bwa, bars3_bwa])
other_bars_bwa = total_bars_bwa.copy()
total_bars_gem = pd.DataFrame([bars1_gem, bars2_gem, bars3_gem])
other_bars_gem = total_bars_gem.copy()
total_bars_bwa.columns = ["BWA - MOD - MAPQ20", 
                      "BWA - MOD - MAPQ40", 
                      "BWA - MOD - MAPQ60", 
                      "BWA - FILT - MAPQ20", 
                      "BWA - FILT - MAPQ40",
                      "BWA - FILT - MAPQ60",
                      ]
total_bars_gem.columns = ["GEM - MOD - MAPQ20", 
                      "GEM - MOD - MAPQ40", 
                      "GEM - MOD - MAPQ60",
                      "GEM - FILT - MAPQ20",
                      "GEM - FILT - MAPQ40",
                      "GEM - FILT - MAPQ60",
]
for column in total_bars_bwa.columns.tolist():
    total_bars_bwa[column] = total_bars_bwa[column].apply(lambda x: x / total_bars_bwa[column].sum() * 100)
for column in total_bars_gem.columns.tolist():
    total_bars_gem[column] = total_bars_gem[column].apply(lambda x: x / total_bars_gem[column].sum() * 100)

# Plot the percentage for BWA experiments.

plt.figure(figsize=(10,6))
plt.rc("font")

bars1 = total_bars_bwa.iloc[0,:]
bars2 = total_bars_bwa.iloc[1,:]
bars3 = total_bars_bwa.iloc[2,:]

names = ["BWA - MOD - MAPQ20", 
         "BWA - MOD - MAPQ40", 
         "BWA - MOD - MAPQ60", 
         "BWA - FILT - MAPQ20", 
         "BWA - FILT - MAPQ40",
         "BWA - FILT - MAPQ60",
        ]
barWidth= 1

# Heights of bars1 + bars2
bars = np.add(bars1, bars2).tolist()

r = [x for x in range(len(names))]

# Create brown bars
plt.bar(r, bars1, color='#72A2C0', edgecolor='white', width=barWidth)
# Create green bars (middle), on top of the firs ones
plt.bar(r, bars2, bottom=bars1, color='#1D65A6', edgecolor='white', width=barWidth)
# Create green bars (top)
plt.bar(r, bars3, bottom=bars, color='#192E5B', edgecolor='white', width=barWidth)

# Custom X axis
plt.xticks(r, names)
plt.legend(["TP = common", "FN = missed", "FP = novel"], loc = "upper right")
#plt.title("BWA variant profile in percentages using categorical MAPQ scores.")
plt.ylabel("Percentage of variants")
#plt.xlabel("MAPQ score")
plt.xticks(fontsize = 10, rotation = 45)
plt.tight_layout()
# Show graphic
plt.savefig("./results/perc_bwa.png")
plt.show()


# Plot the percentage for GEM experiments.
plt.figure(figsize=(10,6))
plt.rc("font", weight = "bold")

bars1 = total_bars_gem.iloc[0,:]
bars2 = total_bars_gem.iloc[1,:]
bars3 = total_bars_gem.iloc[2,:]

names = ["GEM - MOD - MAPQ20",
         "GEM - MOD - MAPQ40", 
         "GEM - MOD - MAPQ60",
         "GEM - FILT - MAPQ20",
         "GEM - FILT - MAPQ40",
         "GEM - FILT - MAPQ60",
        ]
barWidth= 1

# Heights of bars1 + bars2
bars = np.add(bars1, bars2).tolist()

r = [x for x in range(len(names))]

# Create brown bars
plt.bar(r, bars1, color='#72A2C0', edgecolor='white', width=barWidth)
# Create green bars (middle), on top of the firs ones
plt.bar(r, bars2, bottom=bars1, color='#1D65A6', edgecolor='white', width=barWidth)
# Create green bars (top)
plt.bar(r, bars3, bottom=bars, color='#192E5B', edgecolor='white', width=barWidth)

# Custom X axis
plt.xticks(r, names, fontweight='bold')
plt.legend(["TP = common", "FN = missed", "FP = novel"], loc = "upper right")
#plt.title("GEM3 variant profile in percentages using categorical MAPQ scores.")
plt.ylabel("Percentage of variants")
plt.xticks(fontsize = 10, rotation = 45)
plt.tight_layout()
#plt.xlabel("MAPQ score")

# Show graphic
plt.savefig("./results/perc_gem.png")

plt.show()


# Plot the count for BWA experiments.
plt.figure(figsize=(10,6))
plt.rc("font", weight = "bold")

bars1 = other_bars_bwa.iloc[0,:]
bars2 = other_bars_bwa.iloc[1,:]
bars3 = other_bars_bwa.iloc[2,:]

names = ["BWA - MOD - MAPQ20", 
         "BWA - MOD - MAPQ40", 
         "BWA - MOD - MAPQ60", 
         "BWA - FILT - MAPQ20", 
         "BWA - FILT - MAPQ40",
         "BWA - FILT - MAPQ60",
        ]
barWidth= 1

# Heights of bars1 + bars2
bars = np.add(bars1, bars2).tolist()

r = [x for x in range(len(names))]

# Create brown bars
plt.bar(r, bars1, color='#72A2C0', edgecolor='white', width=barWidth)
# Create green bars (middle), on top of the firs ones
plt.bar(r, bars2, bottom=bars1, color='#1D65A6', edgecolor='white', width=barWidth)
# Create green bars (top)
plt.bar(r, bars3, bottom=bars, color='#192E5B', edgecolor='white', width=barWidth)

# Custom X axis
plt.xticks(r, names, fontweight='bold')
plt.legend(["TP = common", "FN = missed", "FP = novel"], loc = "upper right")
#plt.title("BWA variant profile using categorical MAPQ scores.")
plt.ylabel("Percentage of variants")
#plt.xlabel("MAPQ score")
plt.xticks(fontsize = 10, rotation = 45)
plt.tight_layout()
# Show graphic
plt.savefig("./results/num_bwa.png")

plt.show()


# Plot the count for GEM experiments.
plt.figure(figsize=(10,6))
plt.rc("font", weight = "bold")

bars1 = other_bars_gem.iloc[0,:]
bars2 = other_bars_gem.iloc[1,:]
bars3 = other_bars_gem.iloc[2,:]

names = ["GEM - MOD - MAPQ20",
         "GEM - MOD - MAPQ40", 
         "GEM - MOD - MAPQ60",
         "GEM - FILT - MAPQ20",
         "GEM - FILT - MAPQ40",
         "GEM - FILT - MAPQ60",
        ]
barWidth= 1

# Heights of bars1 + bars2
bars = np.add(bars1, bars2).tolist()

r = [x for x in range(len(names))]

# Create brown bars
plt.bar(r, bars1, color='#72A2C0', edgecolor='white', width=barWidth)
# Create green bars (middle), on top of the firs ones
plt.bar(r, bars2, bottom=bars1, color='#1D65A6', edgecolor='white', width=barWidth)
# Create green bars (top)
plt.bar(r, bars3, bottom=bars, color='#192E5B', edgecolor='white', width=barWidth)

# Custom X axis
plt.xticks(r, names, fontweight='bold')
plt.legend(["TP = common", "FN = missed", "FP = novel"], loc = "upper right")
#plt.title("GEM3 variant profile using categorical MAPQ scores.")
plt.ylabel("Percentage of variants")
#plt.xlabel("MAPQ score")
plt.xticks(fontsize = 10, rotation = 45)
plt.tight_layout()
# Show graphic
plt.savefig("./results/num_gem.png")

plt.show()