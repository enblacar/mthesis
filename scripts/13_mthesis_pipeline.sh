#  __           _       _     _ _____              _   _               _             _            _ _            
# / _\ ___ _ __(_)_ __ | |_  / |___ /_   _ __ ___ | |_| |__   ___  ___(_)___   _ __ (_)_ __   ___| (_)_ __   ___ 
# \ \ / __| '__| | '_ \| __| | | |_ (_) | '_ ` _ \| __| '_ \ / _ \/ __| / __| | '_ \| | '_ \ / _ \ | | '_ \ / _ \
# _\ \ (__| |  | | |_) | |_  | |___) |  | | | | | | |_| | | |  __/\__ \ \__ \ | |_) | | |_) |  __/ | | | | |  __/
# \__/\___|_|  |_| .__/ \__| |_|____(_) |_| |_| |_|\__|_| |_|\___||___/_|___/ | .__/|_| .__/ \___|_|_|_| |_|\___|
#                |_|                                                          |_|     |_|                        

# Script: 13_mthesis_pipeline.sh
# Author: Enrique Blanco Carmona
# Year: 2019
# Project: MSc Bioinformatics UAB final thesis. 
# Objective: This is the whole mthesis wrapper pipeline. It just loads each script, so all the work mthesis work 
# can be reproduced with this. 


nohup bash ./scripts/01_data_download.sh &
wait

nohup bash ./scripts/02_variant_calling_pipeline.sh &
wait 

nohup python3.7 ./scripts/03_pipeline_results.py &
wait 

nohup bash ./scripts/04_mapq_mod_and_filt.sh &
wait 

nohup python3.7 ./scripts/05_mapq_mod_and_filt_results.py &
wait

nohup python3.7 ./scripts/06_ML_dataset_generator.py &
wait

nohup bash ./scripts/07_nucleotide_predictors.sh &
wait

nohup python3.7 ./scripts/08_split_tsv.py &
wait

nohup bash ./scripts/09_pileup_predictors.sh &
wait

nohup bash ./scripts/10_wrapper_pileup_predictors.sh
wait

nohup python3.7 ./scripts/11_logreg_analysis.py False &
wait 

nohup python3.7 ./scripts/12_deep_learning_analysis.py True &
wait 

