



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

