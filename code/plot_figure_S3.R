###############################################################################
#
# Figure 
#	main idea - how did we chose the model/features
# 	plot - model performance by taxonomic rank and feature mean decrease AUC
#
# need files:
# 	data/process/ml/lr_performance.tsv
#	data/process/ml/lr_feature_imp.tsv
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
library(tidyverse)
library(cowplot)

###############################################################################
# load data
###############################################################################
performance_df <- read_tsv(here::here('data/process/ml/lr_performance.tsv'),
		col_type = 'dddddddddddddddcdcc')

###############################################################################
#  setup data
###############################################################################
performance_df <- performance_df %>% 
	mutate(dataset = factor(dataset, c('day_0_predict_future_toxin',
			'day_0_moribund', 'day_0_histology')))  %>% 
	select(method, seed, dataset, `Cross\nValidation` = cv_metric_AUC, Test = AUC, taxonomic_level) %>% 
	pivot_longer(cols = c(`Cross\nValidation`, 'Test'), 
		names_to = 'metric', values_to = 'Performance') %>% 
	mutate(taxonomic_level = factor(taxonomic_level, 
		levels = c('OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum')))

###############################################################################
#  plot data
###############################################################################
day_0_toxin_perf_plot <- performance_df %>% 
  filter(dataset == 'day_0_predict_future_toxin') %>% 
  ggplot(aes(x = taxonomic_level, y = Performance, color = metric)) + 
    geom_jitter(alpha = 0.1, position = position_jitterdodge(dodge.width = 0.5)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), 
      position = position_dodge(width = .5)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.25) +
    theme_bw() +
    geom_rect(data = data.frame(ymin = 0.69, ymax = 0.91, method = 'glmnet', taxonomic_level = 2, Performance = 1),
              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
              fill = NA) + 
    labs(x = NULL, y = 'Performance (AUC)', color = NULL,
    	title = 'Toxin') + 
    theme(legend.position = 'none') + 
    coord_flip()

day_0_moribund_perf_plot <- performance_df %>% 
  filter(dataset == 'day_0_moribund') %>% 
  ggplot(aes(x = taxonomic_level, y = Performance, color = metric)) + 
    geom_jitter(alpha = 0.1, position = position_jitterdodge(dodge.width = 0.5)) +   
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), 
      position = position_dodge(width = .5)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.25) +
    theme_bw() +
    geom_rect(data = data.frame(ymin = 0.89, ymax = 1.01, method = 'glmnet', taxonomic_level = 4, Performance = 1),
              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
              fill = NA) + 
    labs(x = NULL, y = 'Performance (AUC)', color = NULL,
    	title = 'Moribundity') + 
    theme(legend.position = 'none') + 
    coord_flip()

day_0_hist_perf_plot <- performance_df %>% 
  filter(dataset == 'day_0_histology') %>% 
  ggplot(aes(x = taxonomic_level, y = Performance, color = metric)) + 
    geom_jitter(alpha = 0.1, position = position_jitterdodge(dodge.width = 0.5)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), 
      position = position_dodge(width = .5)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.25) +
    #facet_grid(dataset~method, scales = 'free_y') + 
    theme_bw() +
    geom_rect(data = data.frame(ymin = 0.95, ymax = 1.02, method = 'glmnet', taxonomic_level = 1, Performance = 1),
              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
              fill = NA) + 
    labs(x = NULL, y = 'Performance (AUC)', color = NULL,
    	title = 'Histopathologic score') + 
    theme(legend.position = 'bottom') + 
    coord_flip()
###############################################################################
#  save plot
###############################################################################
ggsave(here::here('submission/Figure_S3.tiff'),
	plot_grid(
		day_0_toxin_perf_plot,
		day_0_moribund_perf_plot,
		day_0_hist_perf_plot,
		ncol = 1, labels = c('A', 'B', 'C'), rel_heights = c(4,4,5)),  
	height = 7, width = 3, unit = 'in',
  compression = 'lzw')
###############################################################################
