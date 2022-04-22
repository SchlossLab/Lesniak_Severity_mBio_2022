###############################################################################
#
# Figure 4
#	main idea - What community features differentiate outcome?
# 	plot - Relative abundacne distribution by outcome and correlation
#
#	need files (produced by code/run_lefse_analysis.R):
#		data/process/lefse/day_0_toxin_production_Genus.design
#		data/process/lefse/day_0_severity_Genus.design
#		data/process/lefse/day_0_toxin_production_Genus.0.03.lefse_summary
#		data/process/lefse/day_0_severity_Genus.0.03.lefse_summary
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

taxonomy <- read_taxonomy()

histology <- read_histology()

shared <- read_shared() %>% 
	select(-label, -numOtus) %>% 
	pivot_longer(-Group, names_to = 'OTU', values_to = 'counts') %>% 
	mutate(relative_abundance = counts/2107 * 100)

###############################################################################
#  create plot lefse data function
###############################################################################
plot_lefse <- function(file_name){
	lefse_design <- data.table::fread(here(paste0('data/process/lefse/', 
		file_name, '.design')))
		colnames(lefse_design) <- c('Group', 'class')
    
	lefse_df <- data.table::fread(here(paste0('data/process/lefse/', file_name, '.0.03.lefse_summary')),
		fill = T) %>% 
	filter(LDA > 4) %>% 
	#left_join(select(taxonomy, OTU, tax_otu_label), by = 'OTU') %>% 
	mutate(fctr_class = factor(Class),
		#order = LDA + 10 * as.numeric(fctr_class), # group by LDA effect on variable
		order = LDA,
		taxa_label = OTU)
    
	shared %>% 
		left_join(taxonomy, by = 'OTU') %>% 
		right_join(lefse_df, by = c('OTU' = 'OTU')) %>% 
		right_join(lefse_design, by = c('Group')) %>% 
		group_by(Group, tax_otu_label, order, class, LDA, pValue) %>% # removed mouse_id from interactive plot
		summarise(relative_abundance = sum(relative_abundance)) %>% 
		mutate(taxa_label = tax_otu_label,
			relative_abundance = ifelse(relative_abundance == 0, 0.045, 
				relative_abundance),
			class = case_when(class == 'severe' ~ 'Moribund',
				class == 'moderate' ~ 'Non-moribund, High histopathologic score',
				class == 'mild' ~ 'Non-moribund, Low histopathologic score',
				class == FALSE ~ 'Toxin -',
				class == TRUE ~ 'Toxin +', 
				T ~ 'NA'),
			class = factor(class, levels = c('Moribund',
				'Non-moribund, High histopathologic score', 
				'Non-moribund, Low histopathologic score',
				'Toxin +', 'Toxin -')),
			taxa_label = gsub('_', ' ', taxa_label),
			taxa_label = gsub('sensu stricto', 's.s.', taxa_label),
			taxa_label = gsub(' \\(', '\\* \\(', taxa_label),
			taxa_label = paste0('*', taxa_label),
			taxa_label = ifelse(grepl('unclassified', taxa_label), 
				gsub('unclassified\\*', '**', taxa_label),
				taxa_label),
			taxa_label = paste0(taxa_label, '  
LDA = ',
				round(LDA, 3), ', *P* = ', formatC(pValue, format = "e", 0))) %>% 
		ggplot(aes(x = reorder(taxa_label, order), y = relative_abundance, 
				color = class)) +
			stat_summary(fun.data = 'median_hilow', aes(group = class),
				fun.args = (conf.int=0.5), position = position_dodge(width = .5)) +
			geom_hline(yintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
			theme_bw() + 
			labs(x = NULL, y = 'Day 0 Relative Abundance (%)', color = NULL) + 
			scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
				labels=c("0","0.1","1","10","100")) + 
			theme(
				axis.text.y = ggtext::element_markdown(),
				panel.grid.major.y = element_line(color = 'gray95'),
				panel.grid.major.x = element_line(color = 'gray85'),
				panel.grid.minor.x = element_blank(),
				legend.position = 'top',
				legend.margin=margin(0,0,0,0),
    	    	legend.box.margin=margin(0,-10,-10,-10)) + 
			guides(color = guide_legend(reverse = TRUE)) + 
			coord_flip()
}

###############################################################################
#  analyze data
###############################################################################
day_10_hist_genus_corr <- metadata %>% 
	filter(cdiff_strain == 431,
		day == 10) %>% 
	inner_join(histology, by = 'mouse_id') %>% 
	filter(early_euth == F) %>% 
	select(Group = group, test = summary_score) %>% 
	inner_join(shared, by = c('Group')) %>% 
	group_by(OTU) %>% 
	mutate(sd = sd(counts)) %>% 
	filter(!is.na(OTU),
		sd != 0) %>% 
	nest() %>% 
	mutate(spearman_p = map(.x = data, .f = ~ cor.test(.x$counts, .x$test, 
			method = 'spearman', exact = F)$p.value),
		spearman_rho = map(.x = data, .f = ~ cor.test(.x$counts, .x$test, 
			method = 'spearman', exact = F)$estimate)) %>% # compare cleared vs colonized
	unnest(c(spearman_p, spearman_rho)) %>% 
	mutate(pvalue = p.adjust(spearman_p, method = 'BH', n = nrow(.))) %>% # correct p values
	filter(pvalue < 0.05) %>% # select only those above 0.05 after pvalue correction
	unnest(data) %>% 
	left_join(taxonomy, by = c('OTU')) %>% 
	mutate(relative_abundance = ifelse(relative_abundance == 0, 0.045, 
			relative_abundance))

# Test effect of summry score/donor on community distances w/PERMANOVA
## Create dataframe of variables
#dist_variables <- metadata %>% 
#	filter(day == 10) %>% 
#	left_join(select(histology, mouse_id, summary_score),
#		by = c('mouse_id')) %>% 
#	filter(!is.na(summary_score)) %>% 
#	select(id = group, human_source, summary_score) %>% 
#	mutate(human_source = factor(human_source),
#		summary_score = factor(summary_score)) %>% 
#	arrange(id)
#
## Read in thetayc distances of samples
#d10_dist <- read_beta() %>% 
#	filter(rows %in% samples & columns %in% samples)
#samples <- sort(unique(c(d10_dist$rows, d10_dist$columns)))
#
#dist_variables <- dist_variables %>% 
#	filter(id %in% samples)
#
#distance_matrix <- bind_rows(d10_dist, 
#	select(d10_dist, rows = columns, columns = rows, distances),
#	tibble(rows = samples, columns = samples, distances = 0)) %>% 
#	filter(rows %in% samples & columns %in% samples) %>% 
#	pivot_wider(names_from = 'columns', values_from = 'distances')  %>% 
#	arrange(factor(rows, levels = samples)) %>% 
#	select(samples) %>% 
#	as.matrix
#row.names(distance_matrix) <- samples
#adonis_matrix <- as.dist(distance_matrix)	
#
#adonis(adonis_matrix ~ summary_score, data = dist_variables, permutations = 9999)
##   p = 0.0014
##   R2 = 0.45057
#adonis(adonis_matrix ~ human_source, data = dist_variables, permutations = 9999)
##   p = 1e-4
##   R2 = 0.87304
#adonis(adonis_matrix ~ human_source/summary_score, 
#	data = dist_variables, permutations = 9999)
## human source 
##   p = 1e-4
##   R2 = 0.87304
## human source:summary_score
##   p = 3e-4
##   R2 = 0.08277

###############################################################################
#  plot data
###############################################################################
day_10_hist_genus_corr_plot <- day_10_hist_genus_corr %>% 
	mutate(point_cat = case_when(test < 5 ~ 'Low histopathologic score',
		test > 5 ~ 'High histopathologic score',
		T ~ ' '),
		point_cat = factor(point_cat, c('Low histopathologic score',
			'High histopathologic score', ' ')),
		taxa_label = gsub(' \\(', '\\*  
\\(', tax_otu_label),
		taxa_label = paste0('*', taxa_label),
		taxa_label = ifelse(grepl('unclassified', taxa_label), 
			gsub('_unclassified\\*', '* *', taxa_label),
			taxa_label)) %>% 
	ggplot(aes(x = relative_abundance, y = test, color = point_cat)) + 
		geom_point(shape = 1) +
		scale_y_continuous(breaks = c(0,2,4,6,8,10), 
			labels = c(0,2,4,6,8,10)) + 
		#geom_smooth(method = 'lm', se = F) +
		geom_vline(xintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
		scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100),
			labels=c("0","0.1","1","10","100")) + 
		facet_wrap(taxa_label ~ ., scales = 'free_x', nrow = 2) + 
		theme_bw() + 
		labs(x = 'Day 10 Relative Abundance (%)', y = 'Summary score',
			color = NULL) + 
		scale_color_manual(values = c("#62ebc9", "#399283", 'grey')) + 
		theme(strip.background = element_rect(fill = 'white'),
			panel.grid.minor.x = element_blank(),
			strip.text = ggtext::element_markdown(),
			legend.position = 'top',
				legend.margin=margin(2,0,0,0),
    	    	legend.box.margin=margin(2,-10,-10,-10)) + 
		guides(color = guide_legend(override.aes = list(alpha = 1,
					shape = 16,
					color = c("#62ebc9", "#399283", 'white'))))

day_0_toxin_production_plot <- plot_lefse('day_0_toxin_production_OTU') +
	scale_color_manual(values = c("#6a45c2", "#e4b5ff"))

day_0_severity_Genus_plot <- plot_lefse('day_0_severity_OTU') +
	scale_color_manual(values = c("#0c1b37", "#399283", "#62ebc9")) + 
	theme(legend.direction = 'vertical')

###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_4.tiff'),
	cowplot::plot_grid(day_0_severity_Genus_plot + theme(text=element_text(size = 9)),
		cowplot::plot_grid(day_0_toxin_production_plot + 
				theme(text=element_text(size = 9)),
			day_10_hist_genus_corr_plot + theme(text=element_text(size = 9)),
			ncol = 1, rel_heights = c(3,2), labels = c('B', 'C')),
		nrow = 1, rel_widths = c(6,5), labels = c('A', NULL)),
  height = 6, width = 6.5, unit = 'in',
  compression = 'lzw')

###############################################################################
