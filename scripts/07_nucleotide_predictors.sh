#  __           _       _      ___ _____                     _            _   _     _                           _ _      _                 
# / _\ ___ _ __(_)_ __ | |_   / _ \___  |   _ __  _   _  ___| | ___  ___ | |_(_) __| | ___   _ __  _ __ ___  __| (_) ___| |_ ___  _ __ ___ 
# \ \ / __| '__| | '_ \| __| | | | | / (_) | '_ \| | | |/ __| |/ _ \/ _ \| __| |/ _` |/ _ \ | '_ \| '__/ _ \/ _` | |/ __| __/ _ \| '__/ __|
# _\ \ (__| |  | | |_) | |_  | |_| |/ / _  | | | | |_| | (__| |  __/ (_) | |_| | (_| |  __/ | |_) | | |  __/ (_| | | (__| || (_) | |  \__ \
# \__/\___|_|  |_| .__/ \__|  \___//_/ (_) |_| |_|\__,_|\___|_|\___|\___/ \__|_|\__,_|\___| | .__/|_|  \___|\__,_|_|\___|\__\___/|_|  |___/
#                |_|                                                                        |_|                                            
#                |_|                                                                                                                   

# Script: 06_ML_datasets_generator.py
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: This script is intended to generate a fastq file from each ML dataset, so it can be fed to GEM3 so it can generate nucleotide predictors of
# each mapped read.

ML_DATASETS="./data/machine_learning/datasets"
REFPATH="./data/reference"
TMP="./data/tmp"
PREDICTORS="./data/machine_learning/predictors"

# Generate FASTQ files out of SAM files. 
for file in $ML_DATASETS/*.sam
do
	FILENAME=$(echo $file | rev | cut -d / -f 1 | rev | cut -d . -f 1).fastq
	samtools bam2fq -o $ML_DATASETS/$FILENAME $file &
done
wait 

# GENERATE THE PILEUP PREDICTORS.
for file in $ML_DATASETS/*.fastq
do
	FILENAME=$(echo $file | rev | cut -d / -f 1 | rev | cut -d . -f 1).sam 
	gem-mapper -I $REFPATH/chr1.gem -i $file -o $TMP/$FILENAME --mapq-model=dump-predictors > $PREDICTORS/$filename.nt_predictors &
done

