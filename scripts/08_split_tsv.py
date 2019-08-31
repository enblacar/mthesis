#  __           _       _      ___   ___               _ _ _     _             
# / _\ ___ _ __(_)_ __ | |_   / _ \ ( _ )_   ___ _ __ | (_) |_  | |_ _____   __
# \ \ / __| '__| | '_ \| __| | | | |/ _ (_) / __| '_ \| | | __| | __/ __\ \ / /
# _\ \ (__| |  | | |_) | |_  | |_| | (_) |  \__ \ |_) | | | |_  | |_\__ \\ V / 
# \__/\___|_|  |_| .__/ \__|  \___/ \___(_) |___/ .__/|_|_|\__|  \__|___/ \_/  
#                |_|                            |_|                                       


# Script: 08_split_tsv.py
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: This first script takes as input a reference genome, create an arbitrary number of files of 1 Mb each.
# It is assumed that the input SAM is sorted. 

import sys
import os
import pandas as pd

class tsvfile():
    def __init__(self, tsv, genomesize):
        self.tsv = tsv
        self.genomesize = genomesize
        self.split_tsv(1000000)
    def split_tsv(self, chunksize):
        """Splits a tsv file into chunks according to the size of the genome."""

        os.system("mkdir ./data/machine_learning/chunks")
        times = int(self.genomesize // chunksize)
        mod = self.genomesize % chunksize
        cap = 100
        for i in range(times + (1 if mod != 0 else 0)):
            if i == 0:
                init = i
                end = i * chunksize + chunksize + cap

            elif i == times - (0 if mod != 0 else 1):
                init = i * chunksize - cap
                end = self.genomesize

            else:
                init = i * chunksize - cap
                end = i * chunksize + chunksize + cap

            os.system(f"touch ./data/machine_learning/chunks/{i * chunksize}.tsv")
            os.system(f"""echo "0" > ./data/tmp/NRstatus.txt""")
            os.system(f"""tail -n +$(cat ./data/tmp/NRstatus.txt) {self.tsv} | awk '$1 > {end} {{print NR > "./data/tmp/NRstatus.txt"; exit 1}} $1 >= {init} && $1 <= {end} {{print $0}}' > ./data/machine_learning/chunks/{i * chunksize}.tsv""")
        os.system(f"rm ./data/tmp/NRstatus.txt")

def Main():
    tsv = sys.argv[1]
    genomesize = int(sys.argv[2])
    tsvfile(tsv, genomesize)

Main()
