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
current_seed <- as.numeric(commandArgs(TRUE)[1])
taxonomic_level <- as.character(commandArgs(TRUE)[2])
set.seed <- current_seed
fraction <- 0.6

# read in data
print('Reading in data')
day_10_histology <- read_tsv('data/process/ml/day_10_histology.tsv',
						   col_type = cols(.default = col_double(),
						   				   hist_score = col_character(), 
						   				   toxin_presence = col_character()))
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

day_10_histology <- setup_ml_df(day_10_histology, 'hist_score')


# run logistic regression
new_hp <- list(alpha = 0,
			   lambda = c(1e-4, 1e-3, 1e-2, 0.005, 1e-1, 0.5, 1e-0, 5, 1e1, 1e2, 1e3, 1e4))
print('Running Logistic Regression on day 10 to classify histological severity')
day_10_histology_lr <- run_ml(day_10_histology,
	   'glmnet',
       outcome_colname = 'hist_score',
       training_frac = fraction,
	   hyperparameters = new_hp,
	   find_feature_importance = TRUE,
       seed = current_seed)
day_10_histology_lr_NULL <- run_ml(mutate(day_10_histology, hist_score = sample(hist_score)),
	   'glmnet',
       outcome_colname = 'hist_score',
       training_frac = fraction,
	   hyperparameters = new_hp,
	   find_feature_importance = TRUE,
       seed = current_seed)

# run random forest
new_hp <- list(mtry = c(1:10, 15, 20, 25, 40, 50, 100))
print('Running Random Forest on day 10 to classify histological severity')
day_10_histology_rf <- run_ml(day_10_histology,
	   'rf',
       outcome_colname = 'hist_score',
       training_frac = fraction,
	   hyperparameters = new_hp,
	   find_feature_importance = TRUE,
       seed = current_seed)
day_10_histology_rf_NULL <- run_ml(mutate(day_10_histology, hist_score = sample(hist_score)),
	   'rf',
       outcome_colname = 'hist_score',
       training_frac = fraction,
	   hyperparameters = new_hp,
	   find_feature_importance = TRUE,
       seed = current_seed)


print('Modeling complete, saving data')
# 

model_list <- c('day_10_histology_lr_NULL', 'day_10_histology_lr', 
	'day_10_histology_rf_NULL', 'day_10_histology_rf')

ml_performance <- map_dfr(model_list, function(df_name){
	i <- get(df_name)
	i$performance <- i$performance %>% 
		# convert to numeric in case model results in NA/NaN which defaults as character
		mutate_at(vars("cv_metric_logLoss", "logLoss", "AUC", "prAUC", "Accuracy", 
					   "Kappa", "Mean_F1", "Mean_Sensitivity", "Mean_Specificity", 
					   "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", "Mean_Precision", "Mean_Recall", 
					   "Mean_Detection_Rate", "Mean_Balanced_Accuracy", "seed"), 
			as.numeric) %>% 
		mutate_at(vars('method'), as.character) %>% 
		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
			   taxonomic_level = taxonomic_level,
		   	   null_model = grepl('NULL', df_name))
	})
write_tsv(ml_performance, paste0('data/process/ml/temp/hist_ml_performance_', taxonomic_level, '_', current_seed, '.tsv'))	

ml_feature_imp <- map_dfr(model_list, function(df_name){
	i <- get(df_name)
	if(grepl('lr', df_name)){
		feature_coef <- coef(i$trained_model$finalModel, i$trained_model$bestTune$lambda)
		feature_coef <- data.frame(names = gsub('`', '', rownames(feature_coef)),
				coefficient = feature_coef@x)
	} else {
		feature_coef <- data.frame(names = i$feature_importance$names,
				coefficient = NA)
	}
	i$feature_importance <- i$feature_importance %>% 
		left_join(feature_coef, by = c('names')) %>% 
		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
				seed = current_seed,
				taxonomic_level = taxonomic_level,
		   		null_model = grepl('NULL', df_name))
	})
write_tsv(ml_feature_imp, paste0('data/process/ml/temp/hist_ml_feature_imp_', taxonomic_level, '_', current_seed, '.tsv'))

ml_hp_performance <- map_dfr(model_list, function(df_name){
	i <- get(df_name)$trained_model$results %>% 
		mutate(dataset = gsub('(_rf|_lr)', '', df_name),
			   seed = current_seed,
			   taxonomic_level = taxonomic_level,
		   	   null_model = grepl('NULL', df_name))
	if(any(colnames(i) %in% 'lambda')){
		i %>% 
			select(value = lambda, logLoss, dataset, seed, taxonomic_level, null_model) %>% 
			mutate(model = 'glmnet',
				params = 'lambda')
	
		} else if(any(colnames(i) %in% 'mtry')){
		i %>% 
			select(value = mtry, logLoss, dataset, seed, taxonomic_level, null_model) %>% 
			mutate(model = 'rf',
				params = 'mtry')
		}
	})
write_tsv(ml_hp_performance, paste0('data/process/ml/temp/hist_ml_hp_performance_', taxonomic_level, '_', current_seed, '.tsv'))