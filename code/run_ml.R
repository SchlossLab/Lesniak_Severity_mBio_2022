################################################################################
#
# run_ml.R
#
# This script runs mikropml
# 	outputs
#
# Dependencies...
# * data/process/ml/same_day_toxin.tsv
# * data/process/ml/day_0_predict_future_toxin.tsv
# * data/process/ml/day_0_moribund.tsv
# * data/process/ml/day_10_histology.tsv
#
# Output...
#
################################################################################

# load packages
library('tidyverse')
library('mikropml')

# get seed
current_seed <- as.numeric(commandArgs(TRUE))
set.seed <- current_seed

# read in data
same_day_toxin <- read_tsv('data/process/ml/same_day_toxin.tsv',
						   col_type = cols(.default = col_double(), 
						   				   toxin = col_character()))
day_0_predict_future_toxin <- read_tsv('data/process/ml/day_0_predict_future_toxin.tsv',
						   col_type = cols(.default = col_double(), 
						   				   toxin = col_character()))
day_0_moribund <- read_tsv('data/process/ml/day_0_moribund.tsv',
						   col_type = cols(.default = col_double(), 
						   				   early_euth = col_character(),
						   				   toxin_presence = col_logical()))
day_10_histology <- read_tsv('data/process/ml/day_10_histology.tsv',
						   col_type = cols(.default = col_double(),
						   				   hist_score = col_character(), 
						   				   toxin_presence = col_logical())) %>% 
	filter(hist_score != 'mid')

# preprocess data
same_day_toxin <- preprocess_data(same_day_toxin,
								  outcome_colname = 'toxin')$dat_transformed
day_0_predict_future_toxin <- preprocess_data(day_0_predict_future_toxin,
								  outcome_colname = 'toxin')$dat_transformed
day_0_moribund <- preprocess_data(day_0_moribund,
								  outcome_colname = 'early_euth')$dat_transformed
day_10_histology <- preprocess_data(day_10_histology,
								  outcome_colname = 'hist_score')$dat_transformed



# run logistic regression
new_hp <- list(alpha = 0,
			   lambda = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0, 1e1))
print('Running Logistic Regression on same day toxin presence')
same_day_toxin_lr <- run_ml(same_day_toxin,
	   'glmnet',
       outcome_colname = 'toxin',
	   #training_frac = 0.8,
       #hyperparameters = new_hp,
       find_feature_importance = TRUE,
       seed = current_seed)
print('Running Logistic Regression on day 0 to predict toxin production')
day_0_predict_future_toxin_lr <- run_ml(day_0_predict_future_toxin,
	   'glmnet',
       outcome_colname = 'toxin',
	   #hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)

new_hp <- list(alpha = 0,
			   lambda = c(1e-7, 0.000001, 0.00001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.2))
print('Running Logistic Regression on day 0 to predict severity')
day_0_moribund_lr <- run_ml(day_0_moribund,
	   'glmnet',
       outcome_colname = 'early_euth',
	   #hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)
print('Running Logistic Regression on day 10 to classify histological severity')
day_10_histology_lr <- run_ml(day_10_histology,
	   'glmnet',
       outcome_colname = 'hist_score',
	   #hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)

ml_output <- bind_rows(
	mutate(same_day_toxin_lr, ))

# run random forest

new_hp <- list(mtry = c(1:100))
print('Running Random Forest on same day toxin presence')
same_day_toxin_rf <- run_ml(same_day_toxin,
	   'rf',
       outcome_colname = 'toxin',
	   #training_frac = 0.8,
       #hyperparameters = new_hp,
       find_feature_importance = TRUE,
       seed = current_seed)
print('Running Random Forest on day 0 to predict toxin production')
day_0_predict_future_toxin_rf <- run_ml(day_0_predict_future_toxin,
	   'rf',
       outcome_colname = 'toxin',
	   #hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)
print('Running Random Forest on day 0 to predict severity')
day_0_moribund_rf <- run_ml(day_0_moribund,
	   'rf',
       outcome_colname = 'early_euth',
	   #hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)
print('Running Random Forest on day 10 to classify histological severity')
day_10_histology_rf <- run_ml(day_10_histology,
	   'rf',
       outcome_colname = 'hist_score',
	   #hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)

print('Modeling complete, saving data')
# 
ml_performance <- bind_rows(
	mutate(same_day_toxin_lr$performance, dataset = 'same_day_toxin'),
	mutate(day_0_predict_future_toxin_lr$performance, dataset = 'day_0_predict_future_toxin'),
	mutate(day_0_moribund_lr$performance, dataset = 'day_0_moribund'),
	mutate(day_10_histology_lr$performance, dataset = 'day_10_histology'),
	mutate(same_day_toxin_rf$performance, dataset = 'same_day_toxin'),
	mutate(day_0_predict_future_toxin_rf$performance, dataset = 'day_0_predict_future_toxin'),
	mutate(day_0_moribund_rf$performance, dataset = 'day_0_moribund'),
	mutate(day_10_histology_rf$performance, dataset = 'day_10_histology'))

ml_feature_imp <- bind_rows(
	mutate(same_day_toxin_lr$feature_importance, dataset = 'same_day_toxin'),
	mutate(day_0_predict_future_toxin_lr$feature_importance, dataset = 'day_0_predict_future_toxin'),
	mutate(day_0_moribund_lr$feature_importance, dataset = 'day_0_moribund'),
	mutate(day_10_histology_lr$feature_importance, dataset = 'day_10_histology'),
	mutate(same_day_toxin_rf$feature_importance, dataset = 'same_day_toxin'),
	mutate(day_0_predict_future_toxin_rf$feature_importance, dataset = 'day_0_predict_future_toxin'),
	mutate(day_0_moribund_rf$feature_importance, dataset = 'day_0_moribund'),
	mutate(day_10_histology_rf$feature_importance, dataset = 'day_10_histology'))

ml_hp_performance <- bind_rows(
	same_day_toxin_lr$trained_model$results %>% 
		select(value = lambda, AUC) %>% 
		mutate(model = 'glmnet',
			dataset = 'same_day_toxin',
			params = 'lambda'),
	day_0_predict_future_toxin_lr$trained_model$results %>% 
		select(value = lambda, AUC) %>% 
		mutate(model = 'glmnet',
			dataset = 'day_0_predict_future_toxin',
			params = 'lambda'),
	day_0_moribund_lr$trained_model$results %>% 
		select(value = lambda, AUC) %>% 
		mutate(model = 'glmnet',
			dataset = 'day_0_moribund',
			params = 'lambda'),
	day_10_histology_lr$trained_model$results %>% 
		select(value = lambda, AUC) %>% 
		mutate(model = 'glmnet',
			dataset = 'day_10_histology',
			params = 'lambda'),
	same_day_toxin_rf$trained_model$results %>% 
		select(value = mtry, AUC) %>% 
		mutate(model = 'rf',
			dataset = 'same_day_toxin',
			params = 'mtry'),
	day_0_predict_future_toxin_rf$trained_model$results %>% 
		select(value = mtry, AUC) %>% 
		mutate(model = 'rf',
			dataset = 'day_0_predict_future_toxin',
			params = 'mtry'),
	day_0_moribund_rf$trained_model$results %>% 
		select(value = mtry, AUC) %>% 
		mutate(model = 'rf',
			dataset = 'day_0_moribund',
			params = 'mtry'),
	day_10_histology$trained_model$results %>% 
		select(value = mtry, AUC) %>% 
		mutate(model = 'rf',
			dataset = 'day_10_histology',
			params = 'mtry'))

write_tsv(ml_performance, paste0('data/process/ml/temp/ml_performance_', current_seed, '.tsv')	
write_tsv(ml_feature_imp, paste0('data/process/ml/temp/ml_feature_imp_', current_seed, '.tsv')
write_tsv(ml_hp_performance, paste0('data/process/ml/temp/ml_hp_performance_', current_seed, '.tsv')