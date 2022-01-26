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
	filter(!is.na(LDA)) %>% 
	top_n(20, LDA) %>% 
	mutate(fctr_class = factor(Class),
		order = LDA + 10 * as.numeric(fctr_class),
		tax_otu_label = OTU)
    
	shared %>% 
		left_join(taxonomy, by = 'OTU') %>% 
		right_join(lefse_df, by = c('Genus' = 'OTU')) %>% 
		right_join(lefse_design, by = c('Group')) %>% 
		group_by(Group, Genus, order, class) %>% # removed mouse_id from interactive plot
		summarise(relative_abundance = sum(relative_abundance)) %>% 
		mutate(taxa_label = Genus,
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
			taxa_label = gsub('2', '', taxa_label),
			taxa_label = paste0('*', taxa_label, '*'),
			taxa_label = ifelse(grepl('unclassified', taxa_label), 
				gsub(' unclassified\\*', '* *', taxa_label),
				taxa_label)) %>% 
		ggplot(aes(x = reorder(taxa_label, desc(taxa_label)), y = relative_abundance, 
				color = class)) +
			stat_summary(fun.data = 'median_hilow', aes(group = class),
				fun.args = (conf.int=0.5), position = position_dodge(width = .5)) +
			geom_hline(yintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
			theme_bw() + 
			labs(x = NULL, y = 'Relative Abundance (%)', color = NULL) + 
			scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
				labels=c("0","0.1","1","10","100")) + 
			theme(
				axis.text.y = ggtext::element_markdown(),
				panel.grid.major.y = element_line(color = 'gray95'),
				panel.grid.major.x = element_line(color = 'gray85'),
				panel.grid.minor.x = element_blank(),
				legend.position = 'top',
				legend.margin=margin(2,0,0,0),
    	    	legend.box.margin=margin(2,-10,-10,-10)) + 
			guides(color = guide_legend(reverse = TRUE)) + 
			coord_flip()
}

###############################################################################
#  analyze data
###############################################################################
day_10_hist_genus_corr <- metadata %>% 
	filter(cdiff_strain == 431,
		day == mouse_endpoint) %>% 
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

###############################################################################
#  plot data
###############################################################################
day_10_hist_genus_corr_plot <- day_10_hist_genus_corr %>% 
	mutate(taxa_label = gsub(' \\(', '\\*  
\\(', tax_otu_label),
		taxa_label = paste0('*', taxa_label),
		taxa_label = ifelse(grepl('unclassified', taxa_label), 
			gsub('_unclassified\\*', '* *', taxa_label),
			taxa_label)) %>% 
	ggplot(aes(x = relative_abundance, y = test)) + 
		geom_point(alpha = 0.3) +
		scale_y_continuous(breaks = c(0,2,4,6,8,10), 
			labels = c(0,2,4,6,8,10)) + 
		#geom_smooth(method = 'lm', se = F) +
		geom_vline(xintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
		scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100),
			labels=c("0","0.1","1","10","100")) + 
		facet_wrap(taxa_label ~ ., scales = 'free_x') + 
		theme_bw() + 
		labs(x = 'Relative Abundance (%)', y = 'Summary score') + 
		theme(strip.background = element_rect(fill = 'white'),
			panel.grid.minor.x = element_blank(),
			strip.text = ggtext::element_markdown())

day_0_toxin_production_plot <- plot_lefse('day_0_toxin_production_Genus') +
	scale_color_manual(values = c("#6a45c2", "#e4b5ff"))

day_0_severity_Genus_plot <- plot_lefse('day_0_severity_Genus') +
	scale_color_manual(values = c("#0c1b37", "#399283", "#62ebc9")) + 
	theme(legend.direction = 'vertical')

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_4.jpg'),
	cowplot::plot_grid(day_0_severity_Genus_plot + theme(text=element_text(size = 9)),
		cowplot::plot_grid(day_0_toxin_production_plot + 
				theme(text=element_text(size = 9)),
			day_10_hist_genus_corr_plot + theme(text=element_text(size = 9)),
			ncol = 1, rel_heights = c(3,2), labels = c('B', 'C')),
		nrow = 1, labels = c('A', NULL)),
  height = 7, width = 6.875, unit = 'in')

###############################################################################
