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

LOD_df <- data.frame(x = 1.15, y = 2)
DECD_df <- data.frame(x = 2, y = 2.5, day = 10, 
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
    geom_segment(data = DECD_df, aes(x = 1, xend = 6, y = 2.5, yend = 2.5), 
      size = .25, color = 'black', inherit.aes = F) + 
    geom_label(data = DECD_df, aes(x = 3.5, y = 2.5), label = "Deceased", 
      fill = 'white', color = 'white', inherit.aes = F) + 
    geom_text(data = DECD_df, aes(x = 3.5, y = 2.5), label = "Deceased", 
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