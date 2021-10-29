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

###############################################################################
#  setup data
###############################################################################
toxin <- toxin %>% 
  inner_join(select(metadata, group, human_source, early_euth), by = c('group')) %>% 
  mutate(outcome = ifelse(early_euth == T, 'Severe', 'Non-severe')) %>% 
  filter(!is.na(outcome)) %>%
  mutate(day = as.numeric(toxin_sample_day)) %>% 
  left_join(donor_aes, by = 'human_source') %>% 
  filter(!day %in% c(4,7),
         toxin_sample_type == 'stool') 

histology <- histology %>% 
  left_join(distinct(select(metadata, mouse_id, early_euth)), by = 'mouse_id') %>% 
  filter(!is.na(early_euth)) %>% 
  mutate(early_euth = ifelse(early_euth, 'Severe', 'Non-severe'),
         mouse_id = paste0(cage_id, '_', ear_tag),
         score = summary_score) %>% 
  left_join(donor_aes, by = 'human_source')

###############################################################################
#  plot data
###############################################################################
toxin_plot <- toxin %>% 
  ggplot(aes(x = donor_labels, y = Log_repiricoal_dilution, 
             fill = donor_colors, color = donor_colors)) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 2) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    labs(x = NULL, y = 'Toxin activity') + 
    facet_grid(day~outcome, scales = 'free_x', space = 'free_x') + 
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_blank())


histology_plot <- histology %>%  
  ggplot(aes(x = donor_labels, y = score, 
  		fill = donor_colors, color = donor_colors)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) + 
    theme_bw() + 
    labs(x = NULL, y = 'Clinical Score') +
    scale_y_continuous(breaks = c(0,2,4,6,8,10), 
                       labels = c(0,2,4,6,8,10)) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    facet_grid(.~early_euth, scales = 'free', space = 'free_x') + 
    theme(legend.position = 'none',
          strip.background.x = element_blank(),
          strip.text.x = element_text(size = 12),
          legend.margin=margin(t=-0.3, r=0, b=-0.2, l=0, unit="cm"),
          legend.key.height = unit(0, 'cm'))

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_3.jpg'),
  cowplot::plot_grid(toxin_plot, NULL, histology_plot, ncol = 1, 
                   rel_heights = c(8, 1, 4),
                   labels = c('A', 'B', NULL)),
  width = 6, height = 4)
###############################################################################