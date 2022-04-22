###############################################################################
#
# Figure 3
#	main idea - hWhat is the difference in toxin and histology?
# 	plot - toxin/histopatholic scores by mouse/donor
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

donor_aes <- donor_df

LOD_df <- data.frame(x = 1.15, y = 2)
DECD_df <- data.frame(x = 2, y = 2.5, day = factor('10'), 
  outcome = factor('Moribund', levels = c('Non-moribund', 'Moribund')))

###############################################################################
#  setup data
###############################################################################
toxin <- toxin %>% 
  inner_join(select(metadata, group, human_source, early_euth), by = c('group')) %>% 
  mutate(outcome = ifelse(early_euth == T, 'Moribund', 'Non-moribund'),
         outcome = factor(outcome, levels = c('Non-moribund', 'Moribund'))) %>% 
  filter(!is.na(outcome)) %>%
  mutate(day = as.numeric(toxin_sample_day)) %>% 
  left_join(donor_aes, by = 'human_source') %>% 
  filter(!day %in% c(3,4,7),
         toxin_sample_type == 'stool') %>% 
  mutate(day = ifelse(day == 1, '1*', as.character(day)),
    day = factor(day, levels = c('1*', '2', '10'))) 

histology <- histology %>% 
  left_join(distinct(select(metadata, mouse_id, early_euth)), by = 'mouse_id') %>% 
  filter(!is.na(early_euth)) %>% 
  mutate(early_euth = ifelse(early_euth, 'Moribund', 'Non-moribund'),
         early_euth = factor(early_euth, levels = c('Non-moribund', 'Moribund')),
         mouse_id = paste0(cage_id, '_', ear_tag),
         score = summary_score) %>% 
  left_join(donor_aes, by = 'human_source')

toxin_present <- read_toxin() %>% 
  inner_join(select(metadata, group, human_source, early_euth), by = c('group')) %>% 
  mutate(outcome = ifelse(early_euth == T, 'Moribund', 'Non-moribund'),
         outcome = factor(outcome, levels = c('Non-moribund', 'Moribund'))) %>% 
  filter(!is.na(outcome)) %>%
  mutate(day = as.numeric(toxin_sample_day)) %>% 
  left_join(donor_aes, by = 'human_source') %>% 
  filter(day < 4) %>% 
  group_by(mouse_id, outcome) %>% 
  summarize(toxin = max(Log_repiricoal_dilution) > 1) %>% 
  ungroup %>% 
  mutate(toxin = ifelse(toxin, 'Toxin +', 'Toxin -'))

###############################################################################
#  analyze data
###############################################################################

# difference in toxin present in mice
# chisq.test(table(toxin_present$outcome, toxin_present$toxin))
# X-squared = 7.0964, df = 1, p-value = 0.007724

# day 1 diff
#read_toxin() %>% 
#  inner_join(select(metadata, group, human_source, early_euth), by = c('group')) %>% 
#  mutate(outcome = ifelse(early_euth == T, 'Moribund', 'Non-moribund'),
#         outcome = factor(outcome, levels = c('Non-moribund', 'Moribund'))) %>% 
#  filter(!is.na(outcome)) %>%
#  mutate(day = as.numeric(toxin_sample_day)) %>% 
#  left_join(donor_aes, by = 'human_source') %>% 
#  filter(day == 1,
#         toxin_sample_type == 'stool')  %>% 
#  summarize(pvalue = 
#    wilcox.test( 
#      pull(filter(., outcome == 'Moribund'), Log_repiricoal_dilution),
#      pull(filter(., outcome == 'Non-moribund'), Log_repiricoal_dilution))$p.value)
#
# p = 0.00104
#

# endpoint histology
#histology %>% 
#  summarize(pvalue = 
#    wilcox.test( 
#      pull(filter(., early_euth == 'Moribund'), summary_score),
#      pull(filter(., early_euth == 'Non-moribund'), summary_score))$p.value)
# p = 0.000000164
#

#endpoint_summary_score <- histology %>% 
#  filter(early_euth == 'Non-moribund') %>% 
#  pull(summary_score)
#
#diptest::dip.test(endpoint_summary_score)
#  Hartigans' dip test for unimodality / multimodality
#data:  endpoint_summary_score
#D = 0.092593, p-value = 0.04857
#alternative hypothesis: non-unimodal, i.e., at least bimodal
#
#LaplacesDemon::is.unimodal(endpoint_summary_score)
# FALSE
#LaplacesDemon::is.bimodal(endpoint_summary_score)
# TRUE

###############################################################################
#  plot data
###############################################################################
toxin_plot <- toxin %>% 
  ggplot(aes(x = donor_labels, y = Log_repiricoal_dilution, 
             fill = donor_colors, color = donor_colors)) + 
    geom_hline(yintercept = 2, linetype = 'dashed', 
      size = 1, color = 'white') + 
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    labs(x = NULL, y = expression(Toxin ~ Titer ~ Log[10])) + 
    facet_grid(day~outcome, scales = 'free_x', space = 'free_x') + 
    scale_y_continuous(breaks = c(1,2, 3), 
      labels = c(' 1','2','3'),
      sec.axis = dup_axis(labels = NULL, breaks = c(1,2,3),
                                           name = "Days post-challenge")) +
    theme_bw() +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor.y = element_blank()) + 
    geom_label(data = LOD_df, aes(x = x-.25, y = y), label = "LOD", 
      fill = 'white', size = 10/.pt, color = 'white', inherit.aes = F) + 
    geom_text(data = LOD_df, aes(x = x, y = y), label = "LOD", 
      color = 'lightgray', inherit.aes = F) + 
    geom_segment(data = DECD_df, aes(x = 1, xend = 1.9, y = 1.5, yend = 1.5), 
      size = .25, color = 'black', inherit.aes = F) + 
    geom_segment(data = DECD_df, aes(x = 5.1, xend = 6, y = 1.5, yend = 1.5), 
      size = .25, color = 'black', inherit.aes = F) + 
    geom_text(data = DECD_df, aes(x = 3.5, y = 1.5), label = "Deceased", 
      color = 'black', inherit.aes = F)

histology_plot <- histology %>%  
  ggplot(aes(x = donor_labels, y = score, 
  		fill = donor_colors, color = donor_colors)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) + 
    theme_bw() + 
    labs(x = NULL, y = 'Summary Score') +
    scale_y_continuous(breaks = c(0,2,4,6,8,10), 
                       labels = c(0,2,4,6,8,10),
                       sec.axis = dup_axis(labels = NULL, 
                                           breaks = c(0,2,4,6,8,10),
                                           name = "Endpoint*")) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    facet_grid(.~early_euth, scales = 'free', space = 'free_x') + 
    theme(legend.position = 'none',
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor.y = element_blank())

toxin_diff_plot <- toxin_present %>% 
  ggplot(aes(outcome, fill = toxin)) + 
    geom_bar(position = position_dodge()) + 
    scale_fill_manual(values = c("#e4b5ff", "#6a45c2")) + 
    labs(x = NULL, y = 'Number of Mice', fill = NULL) + 
    theme_bw() + 
    theme(legend.position = 'top',
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank())

histology %>% 
  ggplot(aes(summary_score)) + 
    geom_histogram(bins = 11) + 
    facet_grid(early_euth~.) +
    theme_bw() + 
    scale_y_continuous(breaks = 0:9) + 
    scale_x_continuous(breaks = 0:10) + 
    theme(panel.grid.minor = element_blank(),
      strip.background.y = element_blank())
    
###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_3.tiff'),
  cowplot::plot_grid(
    cowplot::plot_grid(NULL, NULL, ncol = 1,
      rel_heights = c(2, 1.1), labels = c('A', 'B')),
    cowplot::plot_grid(toxin_plot,  
      cowplot::plot_grid(histology_plot, NULL, 
        nrow = 1, rel_widths = c(10, .525)), # histology plot width
      ncol = 1, rel_heights = c(8, 4)), # height of plots
    nrow = 1, rel_widths = c(1, 23)), # width of plots for labels
  width = 5, height = 4.5,, unit = 'in',
  compression = 'lzw')

ggsave(here('submission/Figure_S1.tiff'),
  toxin_diff_plot,
  width = 3, height = 5, unit = 'in',
  compression = 'lzw')
###############################################################################