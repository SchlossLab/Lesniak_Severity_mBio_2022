################################################################################
#
# setup_ml_data.R
#
# This script creates the files to be used to run with mikropml
# 	outputs a dataframe with the outcome and associated features
#
# Dependencies...
# * data/process/metadata_tidy.tsv
# * data/process/toxin_tidy.tsv
# * data/process/histology_tidy.tsv
# * data/mothur/sample.final.0.03.subsample.shared
#
# Output...
# * data/process/ml_data/same_day_toxin.tsv
# * data/process/ml_data/day_0_predict_future_toxin.tsv
# * data/process/ml_data/day_0_moribund.tsv
# * data/process/ml_data/day_10_histology.tsv
#
################################################################################

# load packages
library('tidyverse')

# read in data
metadata <- read_tsv('data/process/metadata_tidy.tsv',
					 col_type = 'cccccDddddDcdddddcdl') %>% 
	filter(cdiff_strain == 431) %>% 
	mutate(cdiff_cfu = ifelse(cdiff_cfu == 0, log10(60), log10(cdiff_cfu)))

toxin <- read_tsv('data/process/toxin_tidy.tsv',
				  col_type = 'cccccdc') %>% 
	mutate(toxin = ifelse(Log_repiricoal_dilution > 1, 'present', 'absent'))

histology <- read_tsv('data/process/histology_tidy.tsv',
					  col_type = 'ccccdddcdcc') %>% 
	mutate(hist_score = case_when(summary_score < 5 ~ 'low',
								  summary_score > 5 ~ 'high',
								  summary_score == 5 ~ 'mid',
								  T ~ 'NA'))

shared <- read_tsv('data/mothur/sample.final.0.03.subsample.shared',
         col_types = cols(.default = col_double(),
                          Group = col_character())) %>% 
  mutate(Group = ifelse(grepl('_D|DA', Group), Group,
      gsub('D', '_D', Group))) %>% 
  select(-label, -numOtus)



# setup data to predict toxin level (microbiome effect on toxin)
#	from day of sample
same_day_toxin_df <- toxin %>% 
	select(group, toxin) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	left_join(select(metadata, group, C_difficile_CFU = cdiff_cfu), by = c('group')) %>% 
	select(-group)
write_tsv(same_day_toxin_df, 'data/process/ml/same_day_toxin.tsv')	

#	from day 0?
day_0_predict_future_toxin_df <- toxin %>% 
	group_by(mouse_id) %>% 
	summarise(toxin = ifelse(max(Log_repiricoal_dilution) > 1, 'present', 'absent')) %>% 
	left_join(select(metadata, mouse_id, day, group), 
			  by = c('mouse_id')) %>% 
	filter(day == 0) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	select(-mouse_id, -day, -group)
write_tsv(day_0_predict_future_toxin_df, 'data/process/ml/day_0_predict_future_toxin.tsv')	

# setup data to predict moribunidty (microbiome, cdiff, toxin)
#toxin_summary <- toxin %>% 
#	filter(toxin_sample_day  %in% c(1,2)) %>% 
#	group_by(mouse_id) %>% 
#	summarise(toxin_positive_days = sum(Log_repiricoal_dilution > 1))
#cfu_summary <- metadata %>% 
#	filter(day %in% c(1,2)) %>% 
#	group_by(mouse_id, early_euth) %>% 
#	summarise(median_CFU = median(cdiff_cfu))
day_0_moribund <- metadata %>% 
	select(group, mouse_id, day, early_euth) %>% 
	filter(day == 0) %>% 
#	left_join(cfu_summary, by = c('mouse_id')) %>% 
#	left_join(toxin_summary, by = c('mouse_id')) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	mutate(early_euth = ifelse(early_euth, 'severe', 'mild')) %>% 
	select(-mouse_id, -day, -group) %>% 
	relocate(early_euth)
write_tsv(day_0_moribund, 'data/process/ml/day_0_moribund.tsv')	

# setup data to predict histology severity of colonized mice (microbiome, cdiff, toxin)
#toxin_summary <- toxin %>% 
#	group_by(mouse_id) %>% 
#	summarise(toxin_positive_days = sum(Log_repiricoal_dilution > 1))
day_0_histology_df <- metadata %>% 
	select(group, mouse_id, day, early_euth) %>% 
	filter(day == 0) %>% 
#	left_join(cfu_summary, by = c('mouse_id')) %>% 
#	left_join(toxin_summary, by = c('mouse_id')) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	#ungroup %>% 
	filter(early_euth == F) %>% 
	inner_join(select(histology, mouse_id, hist_score), by = c('mouse_id')) %>% 
	select(-mouse_id, -group, -early_euth, -day) %>% 
	filter(hist_score != 'mid') %>% 
	relocate(hist_score)
write_tsv(day_0_histology_df, 'data/process/ml/day_0_histology.tsv')