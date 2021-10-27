###############################################################################
#
# Figure 
# 	plot -
#	main idea - 
#
# need files:
# 	data/process/
#	data/mothur/
#	data/process/
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
performance_df <- read_tsv(here('data/process/ml/ml_performance.tsv'),
                     col_type = 'dddddddddddddddcdcc') %>% 
  filter(dataset != 'day_0_histology_50') %>% 
  mutate(dataset = gsub('_80', '', dataset),
         dataset = factor(dataset, c('same_day_toxin', 'day_0_predict_future_toxin',
                                     'day_0_moribund', 'day_0_histology')),
         taxonomic_level = factor(taxonomic_level, 
                                  levels = c('OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum')))  

features_df <- read_tsv(here('data/process/ml/ml_feature_imp.tsv'),
                     col_type = 'ddcccddcc') %>% 
  filter(dataset != 'day_0_histology_50') %>% 
  mutate(dataset = gsub('_80', '', dataset))

###############################################################################
#  setup data
###############################################################################
day_0_toxin_production <- features_df %>% 
  filter(dataset == 'day_0_predict_future_toxin',
         (method == 'rf' & taxonomic_level == 'Phylum')) %>% 
  group_by(names) %>% 
  mutate(median_diff = median(perf_metric_diff)) %>% 
  ungroup() %>% 
  filter(median_diff > 0) 

day_0_moribund_plot_df <- features_df %>% 
  filter(dataset == 'day_0_moribund',
         method == 'rf',
         taxonomic_level == 'Class') %>% 
  group_by(names) %>% 
  mutate(median_diff = median(perf_metric_diff)) %>% 
  ungroup() %>% 
  filter(median_diff > 0) 

day_0_histology_plot_df <- features_df %>% 
  filter(dataset == 'day_0_histology',
         method == 'rf',
         taxonomic_level == 'Genus') %>% 
  group_by(names) %>% 
  mutate(median_diff = mean(perf_metric_diff)) %>% 
  ungroup() %>% 
  filter(median_diff > 0) 
###############################################################################
#  analyze data
###############################################################################

###############################################################################
#  plot data
###############################################################################
day_0_toxin_perf_plot <- performance_df %>% 
  filter(dataset == 'day_0_predict_future_toxin') %>% 
  select(method, seed, dataset, cv_auc = cv_metric_AUC, test_auc = AUC, taxonomic_level) %>% 
  pivot_longer(cols = c(cv_auc, test_auc), 
               names_to = 'metric', 
               values_to = 'Performance') %>% 
  mutate(taxonomic_level = factor(taxonomic_level, 
                                  levels = c('OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum'))) %>% 
  ggplot(aes(x = taxonomic_level, y = Performance, color = metric)) + 
    geom_boxplot() + 
    geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.25) +
    facet_grid(dataset~method, scales = 'free_y') + 
    theme_bw() +
#    geom_rect(data = data.frame(ymin = 0.35, ymax = 0.98, method = 'glmnet', taxonomic_level = 2, Performance = 1),
#              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
#              fill = NA) +
    geom_rect(data = data.frame(ymin = 0.48, ymax = 0.99, method = 'rf', taxonomic_level = 6, Performance = 1),
              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
              fill = NA)

day_0_toxin_feature_plot <- day_0_toxin_production %>%
  ggplot(aes(x = reorder(names,median_diff), y = perf_metric_diff)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    #geom_jitter(alpha = 0.1, width = 0.1) + 
    coord_flip() + 
    facet_grid(method~dataset, scales= 'free_y', space = 'free_y') + 
    labs(x = NULL, y = 'Decrease AUC') +
    theme_bw()

day_0_moribund_perf_plot <- performance_df %>% 
  filter(dataset == 'day_0_moribund') %>% 
  select(method, seed, dataset, cv_auc = cv_metric_AUC, test_auc = AUC, taxonomic_level) %>% 
  pivot_longer(cols = c(cv_auc, test_auc), 
               names_to = 'metric', 
               values_to = 'Performance') %>% 
  mutate(taxonomic_level = factor(taxonomic_level, 
                                  levels = c('OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum'))) %>% 
  ggplot(aes(x = taxonomic_level, y = Performance, color = metric)) + 
    geom_boxplot() + 
    geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.25) +
    facet_grid(dataset~method, scales = 'free_y') + 
    theme_bw() +
    geom_rect(data = data.frame(ymin = 0.65, ymax = 1.01, method = 'rf', taxonomic_level = 5, Performance = 1),
              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
              fill = NA)

day_0_moribund_feature_plot <- day_0_moribund_plot_df %>% 
  ggplot(aes(x = reorder(names,median_diff), y = perf_metric_diff)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    #geom_jitter(alpha = 0.1, width = 0.1) + 
    coord_flip() + 
    facet_grid(method~dataset) + 
    labs(x = NULL, y = 'AUC Diff with feature permuted') +
    theme_bw()

day_0_hist_perf_plot <- performance_df %>% 
  filter(dataset == 'day_0_histology') %>% 
  select(method, seed, dataset, cv_auc = cv_metric_AUC, test_auc = AUC, taxonomic_level) %>% 
  pivot_longer(cols = c(cv_auc, test_auc), 
               names_to = 'metric', 
               values_to = 'Performance') %>% 
  mutate(taxonomic_level = factor(taxonomic_level, 
                                  levels = c('OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum'))) %>% 
  ggplot(aes(x = taxonomic_level, y = Performance, color = metric)) + 
    geom_boxplot() + 
    geom_hline(yintercept = 0.5, linetype = 'dashed', size = 0.25) +
    #facet_grid(dataset~method, scales = 'free_y') + 
    theme_bw() +
    geom_rect(data = data.frame(ymin = 0.48, ymax = 1.01, method = 'rf', taxonomic_level = 2, Performance = 1),
              aes(xmin = taxonomic_level - 0.5, xmax = taxonomic_level + 0.5, ymin = ymin, ymax = ymax, color = NA), 
              fill = NA)

day_0_hist_feature_plot <- day_0_histology_plot_df %>% 
  ggplot(aes(x = reorder(names,median_diff), y = perf_metric_diff)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    geom_jitter(alpha = 0.1, width = 0.1) + 
    coord_flip() + 
    facet_grid(method~dataset) + 
    labs(x = NULL, y = 'AUC Diff with feature permuted') +
    theme_bw()

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_S2.jpg'),
  cowplot::plot_grid(
      cowplot::plot_grid(day_0_toxin_perf_plot + coord_flip(), day_0_toxin_feature_plot, nrow = 1),
      cowplot::plot_grid(day_0_moribund_perf_plot + coord_flip(), day_0_moribund_feature_plot, nrow = 1),
      cowplot::plot_grid(day_0_hist_perf_plot + coord_flip(), day_0_hist_feature_plot, nrow = 1),
      ncol = 1),  height = 9, width = 6.875, unit = 'in')
###############################################################################
