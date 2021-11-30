###############################################################################
#
# Figure S1
#	main idea - moribund mice had high levels of 
#     epithelial damage, tissue edema and inflammation
# 	plot - epithelial damage, tissue edema and inflammation
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

histology <- read_histology()

donor_aes <- donor_df

###############################################################################
#  setup data
###############################################################################
histology <- histology %>% 
  left_join(distinct(select(metadata, mouse_id, early_euth)), by = 'mouse_id') %>% 
  filter(!is.na(early_euth)) %>% 
  mutate(early_euth = ifelse(early_euth, 'Moribund', 'Non-moribund'),
         early_euth = factor(early_euth, levels = c('Non-moribund', 'Moribund')),
         mouse_id = paste0(cage_id, '_', ear_tag)) %>% 
  left_join(donor_aes, by = 'human_source') %>% 
  rename('Tissue edema' = edema_tissue,
    'Tissue inflammation' = inflammation_tissue,
    'Epithelial damage' = epithelial_damage) %>% 
  pivot_longer(cols = c('Tissue edema', 'Tissue inflammation', 'Epithelial damage'),
    names_to = 'type', values_to = 'score')

###############################################################################
#  plot data
###############################################################################
histology_plot <- histology %>%  
  ggplot(aes(x = donor_labels, y = score, 
      fill = donor_colors, color = donor_colors)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1) + 
    theme_bw() + 
    labs(x = NULL, y = 'Histopathological Score') +
    scale_y_continuous(breaks = c(0,2,4,6,8,10), 
                       labels = c(0,2,4,6,8,10)) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    facet_grid(type~early_euth, scales = 'free_x', space = 'free_x') + 
    theme(legend.position = 'none',
          strip.background.x = element_blank(),
          strip.text.x = element_text(size = 12),
          legend.margin=margin(t=-0.3, r=0, b=-0.2, l=0, unit="cm"),
          legend.key.height = unit(0, 'cm'))

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_S1.jpg'),
  histology_plot,
  width = 6, height = 4.5)
###############################################################################