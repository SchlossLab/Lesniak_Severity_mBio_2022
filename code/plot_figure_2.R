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
metadata <- read_metadata()

donor_aes <- donor_df

###############################################################################
#  setup data
###############################################################################
LOD_df <- data.frame(x = 14.5, y = 2)
DECD_df <- data.frame(x = 12, y = 3, day = as.factor('Day 10'))

cfu_data <- metadata %>% 
	filter(cdiff_strain == 431,
		day %in% c(1,10),
		!is.na(cdiff_cfu)) %>% 
	mutate(day = factor(paste('Day', day), levels = c('Day 1', 'Day 10')),
		cdiff_cfu = ifelse(cdiff_cfu > 0, log10(cdiff_cfu), log10(60))) %>% 
	left_join(donor_aes, by = c('human_source'))

###############################################################################
#  plot data
###############################################################################
cfu_plot <- cfu_data %>% 
	ggplot(aes(x = donor_labels, y = cdiff_cfu, color = donor_colors)) +
		geom_jitter(width = 0.2) + 
		scale_color_identity() + 
		facet_grid(day~.) + 
		scale_y_continuous(limits = c(1.5,8.5),
			breaks = c(2,4,6,8),
			labels = c('10^2', '10^4', '10^6', '10^8')) + 
		theme_bw() + 
		labs(x = 'Donor Source', y = '*C. difficile* CFU') + 
		guides(color = 'none') + 
		theme(axis.text.y = ggtext::element_markdown(),
			axis.title.y = ggtext::element_markdown(),
			strip.text.y = element_text(angle = 0, size = 12),
			strip.background = element_blank(),
			panel.spacing = unit(2, "lines")) + 
		geom_hline(yintercept = 2, linetype = 'dashed', 
			size = .5, color = 'white') + 
		geom_label(data = LOD_df, aes(x = x, y = y), label = "LOD", 
			fill = 'white', color = 'white', inherit.aes = F) + 
		geom_text(data = LOD_df, aes(x = x, y = y), label = "LOD", 
			color = 'lightgray', size = 12/.pt, inherit.aes = F) + 
		geom_segment(data = DECD_df, aes(x = 10, xend = 15, y = 3, yend = 3), 
			size = .25, color = 'black') + 
		geom_label(data = DECD_df, aes(x = 12.5, y = 3), label = "Deceased", 
			fill = 'white', color = 'white', inherit.aes = F) + 
		geom_text(data = DECD_df, aes(x = 12.5, y = 3), label = "Deceased", 
			color = 'black', size = 12/.pt, inherit.aes = F)

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_2.jpg'),
       cfu_plot,
       width = 6, height = 4)
###############################################################################