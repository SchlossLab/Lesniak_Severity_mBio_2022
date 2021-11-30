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
source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata()
toxin <- read_toxin()
histology <- read_histology()
taxonomy_df <- read_taxonomy()
shared <- read_shared() %>% 
	pivot_longer(cols = starts_with('Otu'), names_to = 'OTU', values_to = 'counts') %>% 
	left_join(taxonomy_df, by = "OTU") %>% 
	mutate(counts = counts/2107 * 100)
	

features_df <- read_tsv(here('data/process/ml/ml_feature_imp.tsv'),
		col_type = 'ddcccddcc') %>% 
	filter(dataset != 'day_0_histology_50') %>% 
	mutate(dataset = gsub('_80', '', dataset))

###############################################################################
#	setup data
###############################################################################
# setup data from model prediction toxin production from day 0
day_0_toxin_features <- features_df %>% 
		filter(dataset == 'day_0_predict_future_toxin',
			method == 'rf',	taxonomic_level == 'Phylum') %>% 
		group_by(names) %>% 
		mutate(mean_diff = mean(perf_metric_diff),
			median_diff = median(perf_metric_diff)) %>% 
		ungroup() %>% 
		filter(mean_diff > 0 & median_diff > 0) %>% 
		select(names, median_diff) %>% 
		unique %>% 
		arrange(median_diff) %>% 
		mutate(taxa_label = paste0('*', names, '*'))
day_0_predict_future_toxin_df <- toxin %>% 
	group_by(mouse_id) %>% 
	summarise(toxin = ifelse(max(Log_repiricoal_dilution) > 1, 
		'Toxigenic', 'Non-toxigenic')) %>% 
	left_join(select(metadata, mouse_id, day, group), 
		by = c('mouse_id')) %>% 
	filter(day == 0) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	inner_join(day_0_toxin_features, by = c('Phylum' = 'names')) %>% 
	group_by(group, toxin, taxa_label) %>% 
	summarise(counts = sum(counts)) %>% 
	mutate(taxon = factor(taxa_label, day_0_toxin_features$taxa_label),
		counts = ifelse(counts == 0, 0.045, counts))

# setup data for model predicting moribund from day 0
day_0_moribund_features <- features_df %>% 
		filter(dataset == 'day_0_moribund',
			method == 'rf',	taxonomic_level == 'Class') %>% 
		group_by(names) %>% 
		mutate(mean_diff = mean(perf_metric_diff),
			median_diff = median(perf_metric_diff)) %>% 
		ungroup() %>% 
		filter(mean_diff > 0 & median_diff > 0) %>% 
		select(names, median_diff) %>% 
		unique %>% 
		arrange(median_diff) %>% 
		mutate(taxa_label = gsub('_', ' ', names),
			taxa_label = paste0('*', taxa_label, '*'),
			taxa_label = gsub(' unclassified\\*', '**', taxa_label))
day_0_moribund <- metadata %>% 
	select(group, mouse_id, day, early_euth) %>% 
	filter(day == 0) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	mutate(early_euth = ifelse(early_euth, 'Moribund', 'Non-moribund')) %>% 
	inner_join(day_0_moribund_features, by = c('Class' = 'names')) %>% 
	group_by(group, early_euth, taxa_label) %>% 
	summarise(counts = sum(counts)) %>% 
	mutate(taxon = factor(taxa_label, day_0_moribund_features$taxa_label),
		counts = ifelse(counts == 0, 0.045, counts))

# setup data for model predicting histology score from day 0
day_0_histology_features <- features_df %>% 
		filter(dataset == 'day_0_histology',
			method == 'rf',	taxonomic_level == 'Genus') %>% 
		group_by(names) %>% 
		mutate(mean_diff = mean(perf_metric_diff)) %>% 
		ungroup() %>% 
		filter(mean_diff > 0) %>% 
		select(names, mean_diff) %>% 
		unique %>% 
		arrange(mean_diff) %>% 
		mutate(taxa_label = gsub('_', ' ', names),
			taxa_label = paste0('*', taxa_label, '*'),
			taxa_label = gsub(' unclassified\\*', '**', taxa_label))
day_0_histology <- metadata %>% 
	select(group, mouse_id, day, early_euth) %>% 
	filter(day == 0) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	ungroup %>% 
	filter(early_euth == F) %>% 
	inner_join(select(histology, mouse_id, summary_score), 
		by = c('mouse_id')) %>% 
	mutate(hist_score = case_when(summary_score > 5 ~ 'High histopathologic score',
		summary_score < 5 ~ 'Low histopathologic score',
		T ~ 'NA')) %>% 
	filter(hist_score %in% c('High histopathologic score', 
		'Low histopathologic score')) %>% 
	inner_join(day_0_histology_features, by = c('Genus' = 'names')) %>% 
	group_by(group, hist_score, taxa_label) %>% 
	summarise(counts = sum(counts)) %>% 
	mutate(taxon = factor(taxa_label, day_0_histology_features$taxa_label),
		counts = ifelse(counts == 0, 0.045, counts))

###############################################################################
#	plot data
###############################################################################
day_0_toxin_feature_abundance_plot <- day_0_predict_future_toxin_df %>% 
	ggplot(aes(x = taxon, y = counts, color = toxin)) + 
		geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
			alpha = 0.1) + 
		stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
			position = position_dodge(width = 0.7)) +
		scale_y_log10() + 
		theme_bw() + 
		coord_flip() + 
		scale_color_manual(values = c("#e4b5ff", "#6a45c2")) + 
		labs(x = NULL, y = 'Relative Abundance (%)', color	= NULL) + 
		theme(legend.position = 'top',
			axis.text.y = ggtext::element_markdown(),
			legend.margin=margin(2,0,0,0),
        	legend.box.margin=margin(2,-10,-10,-10))

day_0_moribund_feature_abundance_plot <- day_0_moribund %>% 
	ggplot(aes(x = taxon, y = counts, color = early_euth)) + 
		geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
			alpha = 0.1) + 
		stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
			position = position_dodge(width = 0.7)) +
		scale_y_log10() + 
		theme_bw() + 
		coord_flip() + 
		scale_color_manual(values = c("4EBFA6", "#0c1b37")) + 
		labs(x = NULL, y = 'Relative Abundance (%)', color	= NULL) + 
		theme(legend.position = 'top',
			axis.text.y = ggtext::element_markdown(),
			legend.margin=margin(2,0,0,0),
        	legend.box.margin=margin(2,-10,-10,-10))

day_0_hist_feature_abundance_plot <- day_0_histology %>% 
	ggplot(aes(x = taxon, y = counts, color = hist_score)) + 
		geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
			alpha = 0.1) + 
		stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
			position = position_dodge(width = 0.7)) +
		scale_y_log10() + 
		theme_bw() + 
		coord_flip() + 
		scale_color_manual(breaks = c('Low histopathologic score', 
				'High histopathologic score'),
			values = c("#62ebc9", "#399283")) + 
		labs(x = NULL, y = 'Relative Abundance (%)', color	= NULL) + 
		theme(legend.position = 'top',
			axis.text.y = ggtext::element_markdown(),
			legend.margin=margin(2,0,0,0),
        	legend.box.margin=margin(2,-10,-10,-10))

###############################################################################
#	save plot
###############################################################################
ggsave(here('results/figures/Figure_5.jpg'),
	cowplot::plot_grid(
		 cowplot::plot_grid(NULL,
			day_0_toxin_feature_abundance_plot,
			nrow = 1, rel_widths = c(1, 25)),
			day_0_moribund_feature_abundance_plot,
		 cowplot::plot_grid(NULL,		 
			day_0_hist_feature_abundance_plot,
			nrow = 1, rel_widths = c(1, 20)),
			rel_heights = c(9, 15, 14),
			labels = c('A', 'B', 'C'),
			ncol = 1), 

	height = 9, width = 6, unit = 'in')
###############################################################################