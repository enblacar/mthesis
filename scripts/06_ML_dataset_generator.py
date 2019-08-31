#  __           _       _      ___   __                __       _       _                 _                                     _             
# / _\ ___ _ __(_)_ __ | |_   / _ \ / /_ _    /\/\    / /    __| | __ _| |_ __ _ ___  ___| |_    __ _  ___ _ __   ___ _ __ __ _| |_ ___  _ __ 
# \ \ / __| '__| | '_ \| __| | | | | '_ (_)  /    \  / /    / _` |/ _` | __/ _` / __|/ _ \ __|  / _` |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \| '__|
# _\ \ (__| |  | | |_) | |_  | |_| | (_) |  / /\/\ \/ /___ | (_| | (_| | || (_| \__ \  __/ |_  | (_| |  __/ | | |  __/ | | (_| | || (_) | |   
# \__/\___|_|  |_| .__/ \__|  \___/ \___(_) \/    \/\____/  \__,_|\__,_|\__\__,_|___/\___|\__|  \__, |\___|_| |_|\___|_|  \__,_|\__\___/|_|   
#                |_|                                                                            |___/                                         

# Script: 06_ML_datasets_generator.py
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: This scripts is aimed towards the generation of ML datasets. For this, the objective is to obtain a list with the position of each novel
# and missed variant. With this, the SAM filtered to MAPQ60 is parsed and all reads landing in the regions of the variant are written into a new SAM.

import pandas as pd
import sys

class BinarySearch:
    def __init__(self, mut_list, position):
        self.mut_list = mut_list
        self.read_position = position


    def binary(self):
        """This function will perform a binary search of the read into the list
        of mutations."""

        def compare_binary(read_position, midpoint):
            """Compare read position against midpoint."""

            return -1 if self.read_position < midpoint else 1 if self.read_position > midpoint else 0

        # Declare sentinels
        first = 0
        last = len(self.mut_list) - 1
        while first <= last:
            midpoint = (first + last) // 2
            result = compare_binary(self.read_position, self.mut_list[midpoint])
            if result == 0:
                return midpoint
            elif result == -1:
                last = midpoint - 1
            else:
                first = midpoint + 1
        return first - 1 # We want the lowest closest integer if the binary search fails.

class SamFile:

    def __init__(self, file, positions, mutations, type = None):
        self.sam = file
        self.mut_list = positions # POS
        self.mutations = mutations # POS, REF, ALT
        self.type = type
        self.write_sam()

    def write_sam(self):
        """Reads a sam file and given a list of mutations, recovers all reads that
        land in the range of the mutations"""

        with open(self.sam, "rt") as sam, open(f"./data/machine_learning/datasets/{self.type}.sam", "wt") as out, open (f"./data/machine_learning/datasets/{self.type}_record.txt", "wt") as report:
            for line in sam:
                if line.startswith("@"):
                    out.write(line)
                    continue
                # Define the range
                pos = int(line.split("\t")[3])
                pos_end = pos + 100
                # Get the index of the binary search
                search = BinarySearch(self.mut_list, pos)
                result = search.binary()
                result = 0 if result == -1 else result
                write = False
                while True:
                    # Get the position of the mutation.
                    if result > len(self.mut_list) - 1:
                        break
                    item = self.mut_list[result]
                    if item <= pos_end and item >= pos:
                        write = True
                        list_fields = line.split("\t")
                        record = list_fields[0] # Get read ID
                        record_pos = list_fields[3] # Get read pos
                        fields = self.mutations[self.mutations["POS"] == self.mut_list[result]] # Get the variant line
                        fields_pos = str(fields["POS"].tolist()[0])
                        fields_ref = str(fields["REF"].tolist()[0])
                        fields_alt = str(fields["ALT"].tolist()[0])
                        report.write(record + "\t" + record_pos + "\t" + str((int(record_pos) + 100)) + "\t" + fields_pos + "\t" + fields_ref + "\t" + fields_alt + "\n")

                    # Keep iterating until the postion of the mutation excees 100 pos_end.
                    if item > pos + 100:
                        if write == True:
                            out.write(line)
                        break
                    # Index of the following mutation.
                    result += 1

def Main():

    # Load files as dataframes.
    base = pd.read_csv("./data/tsv/gem/mastergem.tsv", sep = "\t")
    mapped = pd.read_csv("./data/filtered/gem/tsv/mastergem_60.tsv", sep = "\t")
    base = base.iloc[:, :4]
    mapped = mapped.iloc[:, :4]

    # See how many variants BWA has called correctly: base vs mapped
    merged = pd.merge(base, mapped, how = "outer", indicator = True)

    mask = merged["_merge"] == "both" # Common.
    mask2 = merged["_merge"] == "right_only" # Novel.
    uncalled = merged[mask]
    novel = merged[mask2]
    positions_uncalled = sorted(uncalled["POS"].apply(lambda x: int(x)).tolist())
    mutations_uncalled = uncalled.iloc[:, 1:4]
    positions_novel = sorted(novel["POS"].apply(lambda x: int(x)).tolist())
    mutations_novel = novel.iloc[:, 1:4]
    
    file = "./data/masterbam/gem/mastergem.sam"
    SamFile(file, positions_uncalled, mutations_uncalled, type = "common")
    SamFile(file, positions_novel, mutations_novel, type = "novel")

Main()



