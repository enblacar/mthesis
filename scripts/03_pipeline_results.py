#  __           _       _      ___ _____      ___ _            _ _                                  _ _       
# / _\ ___ _ __(_)_ __ | |_   / _ \___ /_    / _ (_)_ __   ___| (_)_ __   ___   _ __ ___  ___ _   _| | |_ ___ 
# \ \ / __| '__| | '_ \| __| | | | ||_ (_)  / /_)/ | '_ \ / _ \ | | '_ \ / _ \ | '__/ _ \/ __| | | | | __/ __|
# _\ \ (__| |  | | |_) | |_  | |_| |__) |  / ___/| | |_) |  __/ | | | | |  __/ | | |  __/\__ \ |_| | | |_\__ \
# \__/\___|_|  |_| .__/ \__|  \___/____(_) \/    |_| .__/ \___|_|_|_| |_|\___| |_|  \___||___/\__,_|_|\__|___/
#                |_|                               |_|                                                        

# Script: 03_pipeline_results.py
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: Generate a venn diagram out of the results of the variant calling. This is done by loading both dataframes and then merging them using
# the outer parameter so all information is retained. Furthermore, with the indicator parameter, we can see if the variant is present in one dataframe or
# in both. A variant has to be equal if CHROM POS REF ALT is the same. Then, a venn diagram is plotted using matplotlib.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

TSVPATH="./data/tsv"

base_bwa = pd.read_csv(f"{TSVPATH}/bwa/masterbwa.tsv", sep = "\t").iloc[:, :4]
base_gem = pd.read_csv(f"{TSVPATH}/bwa/mastergem.tsv", sep = "\t").iloc[:, :4]
merged_base = pd.merge(base_bwa, base_gem, how = "outer", indicator = True)

def change_values(value):
    if value == "right_only":
        return "GEM"
    elif value == "left_only":
        return "BWA"
    elif value == "both":
        return "Common"

merged_base["_merge"] = merged_base["_merge"].apply(change_values)

gem_values = merged_base["_merge"].value_counts()["GEM"]
bwa_values = merged_base["_merge"].value_counts()["BWA"]
common_values = merged_base["_merge"].value_counts()["Common"]

plt.figure(figsize=(8,8))
plt.rc("font", weight = "bold")
v = venn2(subsets = (bwa_values, gem_values, common_values), set_labels = ("", "", ""))

# Custom it
v.get_patch_by_id('01').set_alpha(0.5)
v.get_patch_by_id('01').set_color('#6699ff')
v.get_patch_by_id('01').set_edgecolor('white')
v.get_patch_by_id('10').set_alpha(0.5)
v.get_patch_by_id('10').set_color('#009999')
v.get_patch_by_id('10').set_edgecolor('white')
v.get_patch_by_id('11').set_alpha(0.1)
v.get_patch_by_id('11').set_color('purple')
v.get_patch_by_id('11').set_edgecolor('white')
v.get_label_by_id('10').set_text(bwa_values)
v.get_label_by_id('11').set_text(common_values)
v.get_label_by_id('01').set_text(gem_values)

# Add title and annotation
# plt.title("NUMBER OF COMMON AND SPECIFIC VARIANTS", fontsize = 16)
 
# Show it
plt.legend(labels = ["BWA", "GEM3", "Common"], loc = "center left")
plt.savefig("./results/venn.png")
plt.show()