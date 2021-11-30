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
metadata <- read_metadata() %>% 
  filter(cdiff_strain == '431')

toxin <- read_toxin()

histology <- read_histology()

donor_aes <- donor_df

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
         toxin_sample_type == 'stool') 

histology <- histology %>% 
  left_join(distinct(select(metadata, mouse_id, early_euth)), by = 'mouse_id') %>% 
  filter(!is.na(early_euth)) %>% 
  mutate(early_euth = ifelse(early_euth, 'Moribund', 'Non-moribund'),
         early_euth = factor(early_euth, levels = c('Non-moribund', 'Moribund')),
         mouse_id = paste0(cage_id, '_', ear_tag),
         score = summary_score) %>% 
  left_join(donor_aes, by = 'human_source')

###############################################################################
#  plot data
###############################################################################
toxin_plot <- toxin %>% 
  ggplot(aes(x = donor_labels, y = Log_repiricoal_dilution, 
             fill = donor_colors, color = donor_colors)) + 
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
          panel.grid.minor.y = element_blank())


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
                                           name = "Endpoint")) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    facet_grid(.~early_euth, scales = 'free', space = 'free_x') + 
    theme(legend.position = 'none',
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor.y = element_blank())

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_3.jpg'),
  cowplot::plot_grid(
    cowplot::plot_grid(NULL, NULL, ncol = 1,
      rel_heights = c(2, 1.1), labels = c('A', 'B')),
    cowplot::plot_grid(toxin_plot,  
      cowplot::plot_grid(histology_plot, NULL, 
        nrow = 1, rel_widths = c(10, .525)), # histology plot width
      ncol = 1, rel_heights = c(8, 4)), # height of plots
    nrow = 1, rel_widths = c(1, 23)), # width of plots for labels
  width = 5, height = 4.5)
###############################################################################