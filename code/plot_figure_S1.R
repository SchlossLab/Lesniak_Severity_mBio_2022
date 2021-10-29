###############################################################################
#
# Figure S1
#	main idea - mice lost weight with challenge
# 	plot - % weight loss through experiment
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

donor_aes <- donor_df

###############################################################################
#  setup data
###############################################################################
wt_loss_df <- metadata %>% 
  filter(cdiff_strain == 431,
         day %in% c(2,10)) %>% 
  mutate(day = factor(paste('Day', day), levels = c('Day 2', 'Day 10'))) %>% 
  left_join(donor_aes, by = 'human_source')

DECD_df <- data.frame(x = 12, y = 75, day = as.factor('Day 10'))
###############################################################################
#  plot data
###############################################################################
wt_loss_plot <- wt_loss_df %>% 
  ggplot(aes(x = donor_labels, y = percent_weightLoss, color = donor_colors)) +
    geom_jitter(width = 0.2) + 
	scale_color_identity() + 
    facet_grid(day~.)+ 
    theme_bw() + 
    labs(x = 'Donor Source', y = 'Weight (% relative to initial)') + 
    guides(color = 'none') + 
	theme(strip.text.y = element_text(angle = 0, size = 12),
		strip.background = element_blank(),
		panel.spacing = unit(2, "lines")) + 
	geom_segment(data = DECD_df, aes(x = 10, xend = 15, y = 75, yend = 75), 
		size = .25, color = 'black') + 
	geom_label(data = DECD_df, aes(x = 12.5, y = 75), label = "Deceased", 
		fill = 'white', color = 'white', inherit.aes = F) + 
	geom_text(data = DECD_df, aes(x = 12.5, y = 75), label = "Deceased", 
		color = 'black', size = 12/.pt, inherit.aes = F)

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_S1.jpg'),
       wt_loss_plot,
       width = 6, height = 4)
###############################################################################