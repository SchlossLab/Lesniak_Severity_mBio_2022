###############################################################################
#
# Figure 5
#	main idea - What are important features when taken in context of each other
# 	plot - relative abundacne by outcome of top model features
#
# need files:
# 	data/process/ml/ml_feature_imp.tsv
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
library(cowplot)

source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
#metadata <- read_metadata()
#toxin <- read_toxin()
#histology <- read_histology()
taxonomy_df <- read_taxonomy()
#shared <- read_shared() %>% 
#	pivot_longer(cols = starts_with('Otu'), names_to = 'OTU', values_to = 'counts') %>% 
#	left_join(taxonomy_df, by = "OTU") %>% 
#	mutate(counts = counts/2107 * 100)
	

features_df <- read_tsv(here('data/process/ml/lr_feature_imp.tsv'),
		col_type = 'ddcccddcc')  %>% 
	mutate(sign = ifelse(coefficient < 0,'negative', 'positive'),
	    p_outcome =  1 - (exp(coefficient)/(1+exp(coefficient))), 
	    OR_outcome = p_outcome/(1-p_outcome),
	    logit_outcome = log(OR_outcome)) %>% 
	group_by(names, dataset, method, taxonomic_level) %>% 
	mutate(median_coef = median(coefficient),
		mean_diff = mean(perf_metric_diff),
		median_diff = median(perf_metric_diff)) %>% 
	ungroup

###############################################################################
#	setup data
###############################################################################
# setup data from model prediction toxin production from day 0
# 2 classes: 'absent', 'present' 
day_0_toxin_features <- features_df %>% 
		filter(dataset == 'day_0_predict_future_toxin',
			method == 'glmnet',	taxonomic_level == 'Genus') %>% 
		filter(mean_diff > 0 & median_diff > 0) %>% 
		select(names, OR_outcome, median_coef) %>% 
		mutate(taxa_label = paste0('*', names, '*'),
			taxa_label = gsub('_unclassified\\*', '**', taxa_label),
			OR_outcome = abs(OR_outcome - 2))

# setup data for model predicting moribund from day 0
# 2 classes: 'mild', 'severe' 
day_0_moribund_features <- features_df %>% 
		filter(dataset == 'day_0_moribund',
			method == 'glmnet',	taxonomic_level == 'Order') %>% 
		filter(mean_diff > 0 & median_diff > 0) %>% 
		select(names, OR_outcome, median_coef) %>% 
		mutate(taxa_label = paste0('*', names, '*'),
			taxa_label = gsub('_unclassified\\*', '**', taxa_label),
			OR_outcome = abs(OR_outcome - 2))

# setup data for model predicting histology score from day 0
# 2 classes: 'high', 'low' 
day_0_histology_features <- features_df %>% 
		filter(dataset == 'day_0_histology',
			method == 'glmnet',	taxonomic_level == 'OTU') %>% 
		filter(mean_diff > 0) %>% 
		select(names, OR_outcome, median_coef) %>% 
		left_join(select(taxonomy_df, OTU, taxa_label = tax_otu_label), 
			by = c('names' = 'OTU'))  %>%    
  	mutate(taxa_label = gsub('_sensu_stricto', ' s.s.', taxa_label),
  		taxa_label = gsub('(.*) \\(OTU (.*)', '*\\1* \\(OTU \\2', taxa_label),
    	taxa_label = gsub('_unclassified\\*', '**', taxa_label),
    	taxa_label = gsub('_', ' ', taxa_label))

###############################################################################
#	plot data
###############################################################################
# Features odds ratios
day_0_toxin_OR_plot <- day_0_toxin_features %>%
  ggplot(aes(x = reorder(taxa_label, median_coef), y = OR_outcome)) + 
  	geom_hline(yintercept = 1, linetype = 'dashed') + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    coord_flip(ylim = c(0.94,1.05)) + 
    labs(x = NULL, y = 'Toxin activity odds ratio') +
    theme_bw() + 
    theme(axis.text.y = ggtext::element_markdown())

day_0_moribund_OR_plot <- day_0_moribund_features %>% 
  ggplot(aes(x = reorder(taxa_label, median_coef), y = OR_outcome)) + 
  	geom_hline(yintercept = 1, linetype = 'dashed') + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    coord_flip() + 
    labs(x = NULL, y = 'Moribundity odds ratio') +
    theme_bw() + 
    theme(axis.text.y = ggtext::element_markdown())

day_0_hist_OR_plot <- day_0_histology_features %>% 
  ggplot(aes(x = reorder(taxa_label, -median_coef), y = OR_outcome)) + 
  	geom_hline(yintercept = 1, linetype = 'dashed') + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    coord_flip() + 
    labs(x = NULL, y = 'High histopathology odds ratio') +
    theme_bw() + 
    theme(axis.text.y = ggtext::element_markdown())

###############################################################################
#	save plot
###############################################################################
ggsave(here('submission/Figure_5.tiff'),
	plot_grid(
		plot_grid(NULL, day_0_toxin_OR_plot, nrow = 1, rel_widths = c(1, 4.17)),
		plot_grid(NULL, day_0_moribund_OR_plot, nrow = 1, rel_widths = c(1, 4.7)),
		day_0_hist_OR_plot, 
		ncol = 1, rel_heights = c(10, 17, 37), labels = c('A', 'B', 'C')),
	height = 8, width = 4, unit = 'in',
  compression = 'lzw')
###############################################################################