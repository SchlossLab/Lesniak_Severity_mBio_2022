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
print('Reading in data')
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
print('Preprocessing data')
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
			   lambda = c(1e-10, 1e-5, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1e0, 1e1, 1e2))
print('Running Logistic Regression on same day toxin presence')
same_day_toxin_lr <- run_ml(same_day_toxin,
	   'glmnet',
       outcome_colname = 'toxin',
	   #training_frac = 0.8,
       hyperparameters = new_hp,
       #find_feature_importance = TRUE,
       seed = current_seed)
new_hp <- list(alpha = 0,
			   lambda = c(1e-10, 1e-5, 1e-0, 1e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 1e4))
print('Running Logistic Regression on day 0 to predict toxin production')
day_0_predict_future_toxin_lr <- run_ml(day_0_predict_future_toxin,
	   'glmnet',
       outcome_colname = 'toxin',
	   hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)

new_hp <- list(alpha = 0,
			   lambda = c(1e-10, 1e-5, 1e-0, 1e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 1e4))
print('Running Logistic Regression on day 0 to predict severity')
day_0_moribund_lr <- run_ml(day_0_moribund,
	   'glmnet',
       outcome_colname = 'early_euth',
	   hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)
new_hp <- list(alpha = 0,
			   lambda = c(1e-10, 1e-5, 1e-0, 1e1, 1e2, 2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 1e4))
print('Running Logistic Regression on day 10 to classify histological severity')
day_10_histology_lr <- run_ml(day_10_histology,
	   'glmnet',
       outcome_colname = 'hist_score',
	   hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)

# run random forest
new_hp <- list(mtry = c(1:10, 20, 30, 40, 50, 75, 100, 150))
print('Running Random Forest on same day toxin presence')
same_day_toxin_rf <- run_ml(same_day_toxin,
	   'rf',
       outcome_colname = 'toxin',
	   #training_frac = 0.8,
       hyperparameters = new_hp,
       find_feature_importance = TRUE,
       seed = current_seed)
print('Running Random Forest on day 0 to predict toxin production')
day_0_predict_future_toxin_rf <- run_ml(day_0_predict_future_toxin,
	   'rf',
       outcome_colname = 'toxin',
	   hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)
print('Running Random Forest on day 0 to predict severity')
day_0_moribund_rf <- run_ml(day_0_moribund,
	   'rf',
       outcome_colname = 'early_euth',
	   hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)
print('Running Random Forest on day 10 to classify histological severity')
day_10_histology_rf <- run_ml(day_10_histology,
	   'rf',
       outcome_colname = 'hist_score',
	   hyperparameters = new_hp,
	   #find_feature_importance = TRUE,
       seed = current_seed)

print('Modeling complete, saving data')
# 

model_list <- c('same_day_toxin_lr', 'day_0_predict_future_toxin_lr', 'day_0_moribund_lr', 
		   'day_10_histology_lr', 'same_day_toxin_rf', 'day_0_predict_future_toxin_rf', 
		   'day_0_moribund_rf', 'day_10_histology_rf')

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
		mutate(dataset = gsub('(_rf|_lr)', '', df_name))
})
write_tsv(ml_performance, paste0('data/process/ml/temp/ml_performance_', current_seed, '.tsv'))	

#ml_feature_imp <- map_dfr(model_list, function(df_name){
#	i <- get(df_name)
#	i$feature_importance <- i$feature_importance %>% 
#		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
#				seed = current_seed)
#})
#write_tsv(ml_feature_imp, paste0('data/process/ml/temp/ml_feature_imp_', current_seed, '.tsv'))

ml_hp_performance <- map_dfr(model_list, function(df_name){
	i <- get(df_name)$trained_model$results %>% 
		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
			   seed = current_seed)
	if(any(colnames(i) %in% 'lambda')){
		i %>% 
			select(value = lambda, AUC, dataset, seed) %>% 
			mutate(model = 'glmnet',
				params = 'lambda')
	
		} else if(any(colnames(i) %in% 'mtry')){
		i %>% 
			select(value = mtry, AUC, dataset, seed) %>% 
			mutate(model = 'rf',
				params = 'mtry')
		}
	})
write_tsv(ml_hp_performance, paste0('data/process/ml/temp/ml_hp_performance_', current_seed, '.tsv'))