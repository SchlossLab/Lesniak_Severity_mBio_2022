###############################################################################
#
# Figure S1
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
metadata <- read_metadata()

###############################################################################
#  setup data
###############################################################################
wt_loss_df <- metadata %>% 
  filter(cdiff_strain == 431,
         day %in% c(2,10)) %>% 
  mutate(human_source = as.numeric(gsub('DA', '', human_source)),
         human_source = factor(human_source, levels = 
                                c('1324', '953', '581', '1134', '10148', 
                                  '369', '430', '10093', '10027', '10034',
                                  '10082', '884', '578', '431', '1245')),
         day = factor(paste('Day', day), levels = c('Day 2', 'Day 10'))) %>% 

###############################################################################
#  plot data
###############################################################################
wt_loss_plot <- wt_loss_df %>% 
  ggplot(aes(x = human_source, y = percent_weightLoss, color = human_source)) +
    geom_jitter(width = 0.2) + 
    facet_grid(day~.)+ 
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    scale_color_manual(values = donor_colors) + 
    theme_bw() + 
    labs(x = 'Donor Source', y = 'Weight (% relative to initial)') + 
    guides(color = 'none')

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_S1.jpg'),
       wt_loss_plot,
       width = 6, height = 4)
###############################################################################