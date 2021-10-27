###############################################################################
#
# Figure 3
#	main idea - 
# 	plot -
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata

toxin <- read_toxin

histology <- read_histology

donor_aes <- donor_df

###############################################################################
#  plot data
###############################################################################
histology_plot <- histology %>% 
  left_join(distinct(select(metadata, mouse_id, early_euth)), by = 'mouse_id') %>% 
  filter(!is.na(early_euth)) %>% 
  mutate(early_euth = ifelse(early_euth, 'Severe', 'Mild'),
         mouse_id = paste0(cage_id, '_', ear_tag)) %>% 
  #left_join(mouse_by_outcome, by = c('mouse_id')) %>%
  pivot_longer(cols = c(edema_tissue, inflammation_tissue, epithelial_damage, summary_score),
               names_to = 'histology', values_to = 'score') %>% 
  filter(histology == 'summary_score') %>% 
  mutate(histology = ifelse(histology == 'summary_score', 'End point', NA),
         human_source = gsub('DA0{0,2}', '', human_source)) %>% 
  ggplot(aes(x = human_source, y = score, , color = early_euth, fill = early_euth)) +
    #stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), position = position_dodge(width = 0.4)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 2) + 
    #geom_point(position = position_jitterdodge(jitter.width = 0.3), alpha = 0.5) + 
    theme_bw() + 
    labs(x = NULL, y = 'Clinical\nScore', fill = 'Outcome') +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    scale_y_continuous(breaks = c(1,3,5,7,9), 
                       labels = c(1,3,5,7,9)) + 
    facet_grid(histology~early_euth, scales = 'free', space = 'free_x') + 
    guides(color = 'none') + 
    theme(legend.position = 'bottom',
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          legend.margin=margin(t=-0.3, r=0, b=-0.2, l=0, unit="cm"),
          legend.key.height = unit(0, 'cm'))

toxin_plot <- toxin %>% 
  inner_join(select(metadata, group, human_source, early_euth), by = c('group')) %>% 
#  left_join(mouse_by_outcome, by = c('mouse_id')) %>% 
  mutate(outcome = ifelse(early_euth == T, 'Severe', 'Mild')) %>% 
  filter(!is.na(outcome)) %>%
  mutate(day = as.numeric(toxin_sample_day),
         human_source = gsub('DA0{0,2}', '', human_source)) %>% 
  filter(!day %in% c(4,7),
         toxin_sample_type == 'stool') %>% 
  ggplot(aes(x = human_source, y = Log_repiricoal_dilution, 
             fill = outcome, color = outcome)) + 
    #stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5), position = position_dodge(width = 0.4)) +
    #geom_jitter(width= 0.3, height = 0, alpha = 0.3) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 2) + 
    labs(x = NULL, y = 'Toxin activity\n ') + 
    facet_grid(day~outcome, scales = 'free_x', space = 'free_x') + 
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_blank())

# presentation plot
#presentation_histology_plot <- histology %>% 
#  left_join(distinct(select(metadata, mouse_id, early_euth)), by = 'mouse_id') %>% 
#  filter(!is.na(early_euth)) %>% 
#  mutate(early_euth = ifelse(early_euth, 'Severe', 'Non-severe'),
#         mouse_id = paste0(cage_id, '_', ear_tag)) %>% 
#  #left_join(mouse_by_outcome, by = c('mouse_id')) %>%
#  pivot_longer(cols = c(edema_tissue, inflammation_tissue, epithelial_damage, summary_score),
#               names_to = 'histology', values_to = 'score') %>% 
#  filter(histology == 'summary_score') %>% 
#  mutate(histology = ifelse(histology == 'summary_score', 'End point', NA),
#         human_source = gsub('DA0{0,2}', '', human_source),
#         human_source = factor(human_source, levels = 
#                                 c('1324', '953', '581', '1134', '10148', 
#                                   '369', '430', '10093', '10027', '10034',
#                                   '10082', '884', '578', '431', '1245'),
#                               labels = LETTERS[1:15])) %>% 
#  ggplot(aes(x = human_source, y = score, fill = human_source, color = human_source)) +
#    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) + 
#    theme_bw() + 
#    labs(x = NULL, y = 'Clinical Score') +
#    scale_y_continuous(breaks = c(1,3,5,7,9), 
#                       labels = c(1,3,5,7,9)) + 
#    scale_fill_manual(values = donor_colors) + 
#    scale_color_manual(values = donor_colors) + 
#    facet_grid(.~early_euth, scales = 'free', space = 'free_x') + 
#    theme(legend.position = 'none',
#          strip.background.x = element_blank(),
#          strip.text.x = element_text(size = 12),
#          legend.margin=margin(t=-0.3, r=0, b=-0.2, l=0, unit="cm"),
#          legend.key.height = unit(0, 'cm'))

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_3.jpg'),
  cowplot::plot_grid(toxin_plot, NULL, histology_plot, ncol = 1, 
                   rel_heights = c(8, 1, 4),
                   labels = c('A', 'B', NULL)),
  width = 6, height = 4)
#
#ggsave('~/Desktop/histology_plot.jpg',
#       presentation_histology_plot,
#       width = 8, height = 3.5)
###############################################################################