################################################################################
#
# run_ml.R
#
# This script runs mikropml
# 	outputs
#
# Dependencies...
# * data/process/ml/day_0_predict_future_toxin.tsv
# * data/process/ml/day_0_moribund.tsv
# * data/process/ml/day_0_histology.tsv
#
# Output...
#
################################################################################

# load packages
library('tidyverse')
library('mikropml')

# get seed
current_seed <- as.numeric(commandArgs(TRUE)[1])
taxonomic_level <- as.character(commandArgs(TRUE)[2])
set.seed(current_seed)
fraction <- as.numeric(commandArgs(TRUE)[3])

# read in data
print('Reading in data')
day_0_predict_future_toxin <- read_tsv('data/process/ml/day_0_predict_future_toxin.tsv',
						   col_type = cols(.default = col_double(), 
						   				   toxin = col_character()))
day_0_moribund <- read_tsv('data/process/ml/day_0_moribund.tsv',
						   col_type = cols(.default = col_double(), 
						   				   early_euth = col_character())) 
day_0_histology <- read_tsv('data/process/ml/day_0_histology.tsv',
						   col_type = cols(.default = col_double(), 
						   				   hist_score = col_character())) 
taxonomy_df <- read_tsv('data/process/final.taxonomy.tidy.tsv',
                        col_type = cols(.default = col_character()))

# preprocess data
print('Preprocessing data')
setup_ml_df <- function(input_df, outcome_column){
	tax_df <- taxonomy_df %>% 
		mutate(tax_level = get(taxonomic_level)) %>% 
		select(OTU, tax_level)
	input_df %>% 
	# select taxonomic level to be used 
		mutate(row = rownames(.)) %>% 
		pivot_longer(cols = contains('Otu'), names_to = 'OTU', values_to = 'value') %>% 
		left_join(tax_df, by = c('OTU' = 'OTU')) %>% 
		group_by(across(c(-OTU, -value))) %>% 
		summarise(value = sum(value)) %>% 
		ungroup() %>% 
		pivot_wider(names_from = tax_level, values_from = value) %>% 
		arrange(as.numeric(row)) %>% 
		select(-row) %>% 
	# remove features present in fewer than one cage (3-4 mice/cage)
		remove_singleton_columns(., threshold = 4) %>% 
		.$dat %>% 
	# remove near zero variance, correlated otus, center feature values
		preprocess_data(., outcome_colname = outcome_column) %>% 
		.$dat_transformed
	}

day_0_predict_future_toxin <- setup_ml_df(day_0_predict_future_toxin, 'toxin')
day_0_moribund <- setup_ml_df(day_0_moribund, 'early_euth')
day_0_histology <- setup_ml_df(day_0_histology, 'hist_score')

# run random forest
new_hp <- list(mtry = c(1:10, 15, 20, 25, 40, 50, 100))
print('Running Random Forest on day 0 to predict toxin production')
day_0_predict_future_toxin_rf <- run_ml(day_0_predict_future_toxin,
	   'rf',
       outcome_colname = 'toxin',
       training_frac = fraction,
	   hyperparameters = new_hp,
       seed = current_seed)
print('Running Random Forest on day 0 to predict severity')
day_0_moribund_rf <- run_ml(day_0_moribund,
	   'rf',
       outcome_colname = 'early_euth',
       training_frac = fraction,
	   hyperparameters = new_hp,
       seed = current_seed)
print('Running Random Forest on day 0 to predict day 10 severity')
day_0_histology_rf <- run_ml(day_0_histology,
	   'rf',
       outcome_colname = 'hist_score',
       training_frac = fraction,
	   hyperparameters = new_hp,
       seed = current_seed)


print('Modeling complete, saving data') 

model_list <- c('day_0_predict_future_toxin_rf', 'day_0_moribund_rf', 'day_0_histology_rf')

ml_performance <- map_dfr(model_list, function(df_name){
	i <- get(df_name)
	i$performance <- i$performance %>% 
		# convert to numeric in case model results in NA/NaN which defaults as character
		mutate_at(vars("cv_metric_AUC", "logLoss", "AUC", "prAUC", "Accuracy", 
					   "Kappa", "F1", "Sensitivity", "Specificity", 
					   "Pos_Pred_Value", "Neg_Pred_Value", "Precision", "Recall", 
					   "Detection_Rate", "Balanced_Accuracy", "seed"), 
			as.numeric) %>% 
		mutate_at(vars('method'), as.character) %>% 
		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
			   taxonomic_level = taxonomic_level,
			   training_fraction = fraction)
	})
write_tsv(ml_performance, paste0('data/process/ml/temp/tune_frac_', fraction*100, '_ml_performance_', taxonomic_level, '_', current_seed, '.tsv'))	

ml_hp_performance <- map_dfr(model_list, function(df_name){
	i <- get(df_name)$trained_model$results %>% 
		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
			   seed = current_seed,
			   taxonomic_level = taxonomic_level,
			   training_fraction = fraction)
	if(any(colnames(i) %in% 'lambda')){
		i %>% 
			select(value = lambda, AUC, dataset, seed, taxonomic_level, training_fraction) %>% 
			mutate(model = 'glmnet',
				params = 'lambda')
	
		} else if(any(colnames(i) %in% 'mtry')){
		i %>% 
			select(value = mtry, AUC, dataset, seed, taxonomic_level, training_fraction) %>% 
			mutate(model = 'rf',
				params = 'mtry')
		}
	})
write_tsv(ml_hp_performance, paste0('data/process/ml/temp/tune_frac_', fraction*100, '_ml_hp_performance_', taxonomic_level, '_', current_seed, '.tsv'))