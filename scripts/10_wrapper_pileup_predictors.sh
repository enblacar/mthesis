#  __           _       _     _  ___                                                       _ _                                       _ _      _                 
# / _\ ___ _ __(_)_ __ | |_  / |/ _ \ _  __      ___ __ __ _ _ __  _ __   ___ _ __   _ __ (_) | ___ _   _ _ __    _ __  _ __ ___  __| (_) ___| |_ ___  _ __ ___ 
# \ \ / __| '__| | '_ \| __| | | | | (_) \ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__| | '_ \| | |/ _ \ | | | '_ \  | '_ \| '__/ _ \/ _` | |/ __| __/ _ \| '__/ __|
# _\ \ (__| |  | | |_) | |_  | | |_| |_   \ V  V /| | | (_| | |_) | |_) |  __/ |    | |_) | | |  __/ |_| | |_) | | |_) | | |  __/ (_| | | (__| || (_) | |  \__ \
# \__/\___|_|  |_| .__/ \__| |_|\___/(_)   \_/\_/ |_|  \__,_| .__/| .__/ \___|_|    | .__/|_|_|\___|\__,_| .__/  | .__/|_|  \___|\__,_|_|\___|\__\___/|_|  |___/
#                |_|                                        |_|   |_|               |_|                  |_|     |_|                                            

# Script: 10_wrapper_pileup_predictors.py
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: This is a wrapper for script 09_pileup_predictors.py. It takes the original ML datasets and divides them into 8 to perform multi-processing.
# Once it finished, it joins all the predictors into a single file.


# Split common dataset.
awk '!/^@/' ./data/machine_learning/datasets/common.sam | awk 'NR % 8 == 1 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_1.sam}}
                        NR % 8 == 2 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_2.sam}} 
                        NR % 8 == 3 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_3.sam}} 
                        NR % 8 == 4 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_4.sam}} 
                        NR % 8 == 5 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_5.sam}} 
                        NR % 8 == 6 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_6.sam}} 
                        NR % 8 == 7 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_7.sam}} 
                        NR % 8 == 0 {{print $0 > ./data/machine_learning/datasets/split_datasets/common_8.sam}}' &

# Split novel dataset.
awk '!/^@/' ./data/machine_learning/datasets/novel.sam | awk 'NR % 8 == 1 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_1.sam}}
                        NR % 8 == 2 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_2.sam}} 
                        NR % 8 == 3 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_3.sam}} 
                        NR % 8 == 4 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_4.sam}} 
                        NR % 8 == 5 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_5.sam}} 
                        NR % 8 == 6 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_6.sam}} 
                        NR % 8 == 7 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_7.sam}} 
                        NR % 8 == 0 {{print $0 > ./data/machine_learning/datasets/split_datasets/novel_8.sam}}' &

# Generate the predictors.
for file in ./data/machine_learning/datasets/split_datasets/common*
do
	python3.6 ./scripts/09_pileup_predictors.py $file &
done
wait 

# Generate the predictors.
for file in ./data/machine_learning/datasets/split_datasets/novel*
do
	python3.6 ./scripts/09_pileup_predictors.py $file &
done
wait 

# Join all splitted predictors.
for file in ./data/machine_learning/predictors/common_*
do
	touch ./data/machine_learning/predictors/common.predictors
	cat $file >> ./data/machine_learning/predictors/common.predictors
	rm $file
done

for file in ./data/machine_learning/predictors/novel_*
do
	touch ./data/machine_learning/predictors/novel.predictors
	cat $file >> ./data/machine_learning/predictors/novel.predictors
	rm $file
done