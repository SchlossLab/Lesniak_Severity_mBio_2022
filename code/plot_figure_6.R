###############################################################################
#
# Figure 6
#	main idea - Is our observation conserved across isolates?
#	plot - CFU and histopathology score with different isolates
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
isolate_levels <- c('     RT         \nIsolate         ', 
	'027\n431', 'UM11\n299', '027\n395','027\n458')

metadata <- read_metadata() %>% 
	mutate(isolate_431 = ifelse(cdiff_strain == 431, '431', 'Other'),
		isolate = ifelse(cdiff_strain == 299, paste0('UM11\n', cdiff_strain),
			paste0('027\n', cdiff_strain)),
		isolate = factor(isolate, levels = isolate_levels))
	
toxin <- read_toxin()

histology <- read_histology() %>% 
	select(summary_score, mouse_id)

donor_aes <- donor_df

donor_levels <- c(' ', 'N7' ,'N6', 'M4')

###############################################################################
#	setup data
###############################################################################
non_431_donors <- metadata %>% 
	filter(cdiff_strain != 431) %>% 
	pull(human_source) %>% 
	unique

cfu_data <- metadata %>% 
	filter(human_source %in% non_431_donors,
				 day == 10) %>% 
	add_row(human_source = 'DA00578', cdiff_cfu = NA, isolate_431 = 'Other',
		isolate = factor('027\n431', levels = isolate_levels)) %>% 
	add_row(human_source = 'DA00578', cdiff_cfu = NA, isolate_431 = 'Other',
		isolate = factor('027\n395', levels = isolate_levels)) %>% 
	left_join(toxin, by = 'group') %>% 
	left_join(donor_aes, by = c('human_source')) %>% 
	add_row(donor_labels = factor(' ', levels = donor_levels),
		isolate = factor('     RT         \nIsolate         ', levels = isolate_levels)) %>% 
	mutate(cdiff_cfu = ifelse(cdiff_cfu == 0, log10(60), log10(cdiff_cfu)),
		donor_labels = factor(donor_labels, levels = donor_levels))

histology_data <- metadata %>% 
	filter(human_source %in% non_431_donors) %>% 
	left_join(histology, by = 'mouse_id') %>% 
	select(cdiff_strain, summary_score, human_source, 
		mouse_id, isolate_431, isolate) %>% 
	unique %>% 
	left_join(donor_aes, by = 'human_source') %>% 
	add_row(donor_labels = factor(' ', levels = donor_levels),
		summary_score = -1,
		isolate = factor('     RT         \nIsolate         ', levels = isolate_levels)) %>% 
	mutate(donor_labels = factor(donor_labels, levels = donor_levels))

###############################################################################
#	plot data
###############################################################################
cfu_plot <- cfu_data %>% 
	ggplot(aes(x = isolate, y = cdiff_cfu, color = isolate_431)) + 
		geom_point(position = position_jitterdodge(jitter.width = 0.1)) + 
		geom_hline(yintercept = 2, linetype = 'dashed', color = 'lightgray') + 
		geom_label(data = data.frame(x = 1, y = 2, 
			donor_labels = factor('N7', levels = donor_levels)), aes(x = x, y = y), 
			label = "LOD", fill = 'white', color = 'white', size = 8/.pt, inherit.aes = F) + 
		geom_text(data = data.frame(x = 1, y = 2, 
			donor_labels = factor('N7', levels = donor_levels)), aes(x = x, y = y), 
			label = "LOD", color = 'lightgray', size = 8/.pt, inherit.aes = F) + 
		geom_text(data = data.frame(x = 1, y = 1.65, 
			donor_labels = factor('M4', levels = donor_levels)), aes(x = x, y = y), 
			label = "DECD", color = "#3588d1", size = 8/.pt, inherit.aes = F) + 
		geom_text(data = data.frame(x = 2, y = 1.65, 
			donor_labels = factor('M4', levels = donor_levels)), aes(x = x, y = y), 
			label = "DECD", color = "#76406b", size = 8/.pt, inherit.aes = F) + 
		scale_color_manual(values = c("#3588d1", "#76406b")) + 
		scale_y_continuous(limits = c(1.5, 8.5),
			breaks = c(2,4,6,8), 
			labels = c('10^2', '10^4', '10^6', '10^8')) + 
		facet_grid(.~donor_labels, scale = 'free_x', margins = c(1, 2, 2)) + 
		theme_bw() + 
		theme(legend.position = 'none',
			panel.spacing.x=unit(c(0,1,1), 'lines'),
			strip.background = element_rect(fill = 'white'),
			axis.text.y = ggtext::element_markdown(),
			axis.title.y = ggtext::element_markdown(),
			axis.ticks.x = element_blank(),
			axis.text.x = element_text(#hjust = c(1, .5),
				size = 8),
			legend.title = ggtext::element_markdown(),
			plot.title = ggtext::element_markdown(hjust = 0.5),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()) + 
		labs(x = NULL, y = '*C. difficile* CFU', 
			title = 'Day 10 post-challenge', subtitle = 'Donor source',
			color = '*C. difficile*	
isolate') +
		ggh4x::force_panelsizes(cols = c(.001, 1, 1, 1), TRUE)

histology_score_plot <- histology_data %>% 
	ggplot(aes(x = isolate, y = summary_score, 
			color = isolate_431, fill = isolate_431)) + 
		geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5) + 
		facet_grid(.~donor_labels, scale = 'free_x') + 
		theme_bw() + 
		scale_color_manual(values = c("#3588d1", "#76406b")) + 
		scale_fill_manual(values = c("#3588d1", "#76406b")) + 
		scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
		theme(legend.position = 'none',
			panel.spacing.x=unit(c(0,1,1), 'lines'),
			strip.background = element_rect(fill = 'white'),
			plot.title = element_text(hjust = 0.5),
			axis.ticks.x = element_blank(),
			axis.text.x = element_text(#hjust = c(1, .5),
				size = 8),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()) + 
		labs(title = 'Endpoint', subtitle = 'Donor source',
			x = NULL, y = 'Histopathologic score') +
		coord_cartesian(ylim = c(0,10)) + 
		ggh4x::force_panelsizes(cols = c(0.001, 1, 1, 1), TRUE)

###############################################################################
#	plot data
###############################################################################
ggsave(here('results/figures/Figure_6.jpg'),
			 cowplot::plot_grid(cfu_plot, histology_score_plot, 
													ncol = 1#, rel_widths = c(6,4)
													),
							height = 6, width = 3.5, unit = 'in')
###############################################################################
