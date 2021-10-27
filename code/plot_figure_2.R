###############################################################################
#
# Figure 2
#	main idea - Mice were colonized without any perturbation
# 	plot - CFU of each mouse at day 0 and day 10
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

donor_aes <- donor_df

###############################################################################
#  setup data
###############################################################################

###############################################################################
#  plot data
###############################################################################
cfu_plot <- metadata %>% 
  filter(cdiff_strain == 431,
         day %in% c(1,10)) %>% 
  mutate(human_source = as.numeric(gsub('DA', '', human_source)),
         human_source = factor(human_source, levels = 
                                c('1324', '953', '581', '1134', '10148', 
                                  '369', '430', '10093', '10027', '10034',
                                  '10082', '884', '578', '431', '1245')),
         day = factor(paste('Day', day), levels = c('Day 1', 'Day 10'))) %>% 
  ggplot(aes(x = human_source, y = log10(cdiff_cfu + 60), color = human_source)) +
    geom_jitter(width = 0.2) + 
    facet_grid(day~.)+ 
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    scale_color_manual(values = donor_colors) + 
    theme_bw() + 
    labs(x = 'Donor Source', y = 'C. difficile CFU (log10)') + 
    guides(color = 'none')

# presentation plot
#presentation_cfu_plot <- metadata %>% 
#  filter(cdiff_strain == 431,
#         day %in% c(1,10)) %>% 
#  mutate(human_source = as.numeric(gsub('DA', '', human_source)),
#         human_source = factor(human_source, levels = 
#                                c('1324', '953', '581', '1134', '10148', 
#                                  '369', '430', '10093', '10027', '10034',
#                                  '10082', '884', '578', '431', '1245'),
#                               labels = LETTERS[1:15]),
#         day = factor(paste('Day', day), levels = c('Day 1', 'Day 10'))) %>% 
#  ggplot(aes(x = human_source, y = cdiff_cfu, color = human_source)) +
#    geom_jitter(width = 0.2) + 
#    facet_grid(day~.) + 
#    scale_color_manual(values = donor_colors) + 
#    scale_y_continuous(limits = c(1.5,8.5),
#                       breaks = c(2,4,6,8),
#                       labels = c('10^2', '10^4', '10^6', '10^8')) + 
#    theme_bw() + 
#    labs(x = 'Donor Source', y = '*C. difficile* CFU') + 
#    guides(color = 'none') + 
#    theme(axis.text.y = ggtext::element_markdown(),
#          axis.title.y = ggtext::element_markdown(),
#          strip.text.y = element_text(angle = 0, size = 12),
#          strip.background = element_blank(),
#          panel.spacing = unit(2, "lines")) + 
#			geom_hline(yintercept = 2, linetype = 'dashed', size = .5, color = 'white') + 
#			geom_label(x = 14, y = 2, label = "LOD", fill = 'white', color = 'white') + 
#			geom_text(x = 14, y = 2, label = "LOD", color = 'lightgray', size = 12/.pt)

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_2.jpg'),
       cfu_plot,
       width = 6, height = 4)
###############################################################################
