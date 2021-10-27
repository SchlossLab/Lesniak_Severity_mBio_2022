###############################################################################
#
# Figure 5
#	main idea - 
# 	plot -
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
source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata()
toxin <- read_toxin()
histology <- read_histology()
taxonomy_df <- read_taxonomy()
shared <- read_shared()


###############################################################################
#  setup data
###############################################################################
# setup data from model prediction toxin production from day 0
day_0_predict_future_toxin_df <- toxin %>% 
	group_by(mouse_id) %>% 
	summarise(toxin = ifelse(max(Log_repiricoal_dilution) > 1, 'present', 'absent')) %>% 
	left_join(select(metadata, mouse_id, day, group), 
			  by = c('mouse_id')) %>% 
	filter(day == 0) %>% 
	inner_join(shared, by = c('group' = 'Group'))
day_0_toxin_features <- day_0_toxin_production %>% 
  select(names, median_diff) %>% 
  unique %>% 
  arrange(median_diff) %>% 
  pull(names)
# setup data for model predicting moribund from day 0
toxin_summary <- toxin %>% 
	filter(toxin_sample_day  %in% c(1,2)) %>% 
	group_by(mouse_id) %>% 
	summarise(toxin_positive_days = sum(Log_repiricoal_dilution > 1))
cfu_summary <- metadata %>% 
	filter(day %in% c(1,2)) %>% 
	group_by(mouse_id, early_euth) %>% 
	summarise(median_CFU = median(cdiff_cfu))
day_0_moribund <- metadata %>% 
	select(group, mouse_id, day) %>% 
	filter(day == 0) %>% 
	left_join(cfu_summary, by = c('mouse_id')) %>% 
	left_join(toxin_summary, by = c('mouse_id')) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	mutate(early_euth = ifelse(early_euth, 'severe', 'mild'))
day_0_moribund_features <- day_0_moribund_plot_df %>% 
  select(names, median_diff) %>% 
  unique %>% 
  arrange(median_diff) %>% 
  pull(names)
# setup data for model predicting histology score from day 0
toxin_summary <- toxin %>% 
	group_by(mouse_id) %>% 
	summarise(toxin_positive_days = sum(Log_repiricoal_dilution > 1))
day_0_histology <- metadata %>% 
	select(group, mouse_id, day) %>% 
	filter(day == 0) %>% 
	left_join(cfu_summary, by = c('mouse_id')) %>% 
	left_join(toxin_summary, by = c('mouse_id')) %>% 
	inner_join(shared, by = c('group' = 'Group')) %>% 
	ungroup %>% 
	filter(early_euth == F) %>% 
	inner_join(select(histology, mouse_id, hist_score), by = c('mouse_id'))
day_0_histology_features <- day_0_histology_plot_df %>% 
  select(names, median_diff) %>% 
  separate(names, into = c('A', 'B')) %>% 
  unique %>% 
  arrange(median_diff) %>% 
  pivot_longer(cols = c('A', 'B'), values_to = 'names') %>% 
  pull(names) %>% 
  unique

###############################################################################
#  plot data
###############################################################################
day_0_toxin_feature_abundance_plot <- day_0_predict_future_toxin_df %>% 
  pivot_longer(cols = starts_with('Otu'), names_to = 'OTU', values_to = 'counts') %>% 
  left_join(taxonomy_df, by = "OTU") %>% 
  filter(Phylum %in% day_0_toxin_features) %>% 
  group_by(group, toxin, Phylum) %>% 
  summarise(counts = sum(counts)/21.07) %>% 
  mutate(taxon = factor(Phylum, day_0_toxin_features)) %>% 
  ggplot(aes(x = taxon, y = counts, color = toxin)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
                 position = position_dodge(width = 0.7)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
                alpha = 0.2) + 
    scale_y_log10() + 
    theme_bw() + 
    coord_flip() + 
    scale_color_manual(values = c('red', 'black'),
                       breaks = c('absent', 'present'), 
                       labels = c('False', 'True')) + 
    labs(x = NULL, y = 'Relative Abundance (%)', color  = 'Toxin Detected')

day_0_moribund_feature_abundance_plot <- day_0_moribund %>% 
  pivot_longer(cols = starts_with('Otu'), names_to = 'OTU', values_to = 'counts') %>% 
  left_join(taxonomy_df, by = "OTU") %>% 
  filter(Class %in% day_0_moribund_features) %>% 
  group_by(group, early_euth, Class) %>% 
  summarise(counts = sum(counts)/21.07) %>% 
  mutate(taxon = factor(Class, day_0_moribund_features)) %>% 
  ggplot(aes(x = taxon, y = counts, color = early_euth)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
                 position = position_dodge(width = 0.7)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
                alpha = 0.2) + 
    scale_y_log10() + 
    theme_bw() + 
    coord_flip() + 
#    scale_color_manual(values = c('red', 'black'),
 #                      breaks = c('absent', 'present'), 
  #                     labels = c('False', 'True')) + 
    labs(x = NULL, y = 'Relative Abundance (%)', color  = 'Disease Severity')

day_0_hist_feature_abundance_plot <- day_0_histology %>% 
  filter(hist_score %in% c('high', 'low')) %>% 
  pivot_longer(cols = starts_with('Otu'), names_to = 'OTU', values_to = 'counts') %>% 
  left_join(taxonomy_df, by = "OTU") %>% 
  filter(Genus %in% day_0_histology_features) %>% 
  group_by(group, hist_score, Genus) %>% 
  summarise(counts = sum(counts)/21.07) %>% 
  mutate(taxon = factor(Genus, day_0_histology_features)) %>% 
  ggplot(aes(x = taxon, y = counts, color = hist_score)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
                 position = position_dodge(width = 0.7)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
                alpha = 0.2) + 
    scale_y_log10() + 
    theme_bw() + 
    coord_flip() + 
#    scale_color_manual(values = c('red', 'black'),
 #                      breaks = c('absent', 'present'), 
  #                     labels = c('False', 'True')) + 
    labs(x = NULL, y = 'Relative Abundance (%)', color  = 'Disease Severity')
###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_5.jpg'),
  cowplot::plot_grid(
#    cowplot::plot_grid(
     cowplot::plot_grid(NULL,
      day_0_toxin_feature_abundance_plot + 
        theme(legend.position = 'right') + 
        labs(color = 'Toxin\nDetected'),
      nrow = 1, rel_widths = c(1, 25)),
#      same_day_toxin_feature_abundance_plot + 
#        theme(legend.position = 'none'),
#      ncol = 1),
#    cowplot::plot_grid(  
      day_0_moribund_feature_abundance_plot + 
        theme(legend.position = 'right') + 
        labs(color = 'Disease\nSeverity'),
     cowplot::plot_grid(NULL,     
      day_0_hist_feature_abundance_plot + 
        theme(legend.position = 'right') + 
        scale_color_manual(values = c('dark green', 'light green')) + 
        labs(color = 'Clinical\nScore'),
      nrow = 1, rel_widths = c(1, 20)),
      rel_heights = c(9, 15, 14),
      labels = c('A', 'B', 'C'),
      ncol = 1), 
#    nrow = 1),
  height = 9, width = 6, unit = 'in')
###############################################################################
