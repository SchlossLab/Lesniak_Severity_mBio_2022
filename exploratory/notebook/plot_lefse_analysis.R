library(tidyverse)

metadata <- read_tsv('data/process/metadata_tidy.tsv',
		col_type = 'cccccDddddDcdddddcdl') 
shared <- read_tsv('data/mothur/sample.final.0.03.subsample.shared',
		col_types = cols(.default = col_double(),
			Group = col_character())) %>% 
	mutate(Group = ifelse(grepl('_D|DA', Group), Group,
		gsub('D', '_D', Group)))
taxonomy <- read_tsv('data/process/final.taxonomy.tidy.tsv',
                     col_types = cols(.default = col_character()))

 # plot lefse results
lefse_df <- data.table::fread('data/process/severe_disease_OTU.0.03.lefse_summary',
		fill = T) %>% 
	filter(!is.na(LDA)) %>% 
	top_n(20, LDA) %>% 
	left_join(select(taxonomy, OTU, tax_otu_label), by = 'OTU') %>% 
	mutate(order = LDA + 10 * as.logical(Class),
		tax_otu_label = gsub('unclassified', 'UC', tax_otu_label))

lefse_plot <- lefse_df %>% 
	ggplot(aes(y = LDA, fill = Class, 
			x = reorder(tax_otu_label, order))) +
		geom_col() + 
		theme_bw() + 
		coord_flip() + 
		guides(fill = 'none') + 
		scale_fill_manual(values = c(1,2),
			limits = c(TRUE, FALSE),
			labels = c('Severe', 'Non-severe')) + 
		labs(x = NULL, y = 'LDA Score (log 10)')

abundance_plot <- metadata %>% 
	filter(cdiff_strain == 431,
		day == 0) %>%   
	left_join(shared, by = c('group' = 'Group')) %>% 
	select(group, early_euth, one_of(lefse_df$OTU)) %>%  
	gather(OTU, abundance, matches('Otu\\d')) %>% 
	left_join(lefse_df, by = 'OTU') %>% 
	mutate(abundance = abundance / 21.07,
		abundance = ifelse(abundance == 0, 0.045, abundance)) %>% 
	ggplot(aes(x = reorder(tax_otu_label, order), 
			y = abundance, color = as.factor(early_euth))) + 
		geom_jitter(position = position_jitterdodge(dodge.width = 0.5), alpha = 0.25) + 
		stat_summary(fun.data = 'median_hilow', aes(group = early_euth),
			fun.args = (conf.int=0.5), position = position_dodge(width = 0.5)) +
		theme_bw() + labs(x = NULL, y = 'Relative Abundance', color = NULL) + 
		scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
			labels=c("0","0.1","1","10","100")) +
		scale_color_manual(values = c(1,2),
						   limits = c(TRUE, FALSE),
						   labels = c('Severe', 'Non-severe')) + 
		geom_hline(yintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
		theme(panel.grid.major.y = element_blank(), 
			panel.grid.minor = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			legend.position = 'right') + 
		coord_flip()

ggsave('exploratory/notebook/lefse_plot_severe_disease_OTU.jpg', 
	cowplot::plot_grid(lefse_plot, abundance_plot, nrow = 1), 
	width = 15, height = 10)


lefse_df <- data.table::fread(paste0('data/process/severe_disease_Genus.0.03.lefse_summary'),
		fill = T) %>% 
	filter(!is.na(LDA)) %>% 
	top_n(20, LDA) %>% 
	mutate(order = LDA + 10 * as.logical(Class),
		tax_otu_label = OTU)

lefse_plot <- lefse_df %>% 
	ggplot(aes(y = LDA, fill = Class, 
			x = reorder(tax_otu_label, order))) +
		geom_col() + 
		theme_bw() + 
		coord_flip() + 
		guides(fill = 'none') + 
		scale_fill_manual(values = c(1,2),
			limits = c(TRUE, FALSE),
			labels = c('Severe', 'Non-severe')) + 
		labs(x = NULL, y = 'LDA Score (log 10)')

abundance_plot <- metadata %>% 
	filter(cdiff_strain == 431,
		day == 0) %>% 
	left_join(shared, by = c('group' = 'Group')) %>% 
	gather(OTU, abundance, matches('Otu\\d*')) %>% 
	left_join(taxonomy, by = 'OTU') %>% 
	right_join(lefse_df, by = c('Genus' = 'OTU')) %>% 
	group_by(group, Genus, order, early_euth) %>% 
	summarise(abundance = sum(abundance)) %>% 
	mutate(tax_otu_label = Genus) %>% 
	mutate(abundance = abundance / 21.07,
		abundance = ifelse(abundance == 0, 0.045, abundance)) %>% 
	ggplot(aes(x = reorder(tax_otu_label, order), 
			y = abundance, color = as.factor(early_euth))) + 
		geom_jitter(position = position_jitterdodge(dodge.width = 0.5), alpha = 0.25) + 
		stat_summary(fun.data = 'median_hilow', aes(group = early_euth),
			fun.args = (conf.int=0.5), position = position_dodge(width = 0.5)) +
		theme_bw() + labs(x = NULL, y = 'Relative Abundance', color = NULL) + 
		scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
			labels=c("0","0.1","1","10","100")) +
		scale_color_manual(values = c(1,2),
						   limits = c(TRUE, FALSE),
						   labels = c('Severe', 'Non-severe')) + 
		geom_hline(yintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
		theme(panel.grid.major.y = element_blank(), 
			panel.grid.minor = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			legend.position = 'right') + 
		coord_flip()

ggsave('exploratory/notebook/lefse_plot_severe_disease_Genus.jpg', 
	cowplot::plot_grid(lefse_plot, abundance_plot, nrow = 1), 
	width = 15, height = 10)
