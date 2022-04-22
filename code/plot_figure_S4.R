###############################################################################
#
# Figure S4
# 	plot - relative abundance over time
#	main idea - which temporal trends are significant by high/low histopath?
#
# need files:
# 	data/process/lefse/temporal_trend_hilo.design
#	data/process/lefse/temporal_trend_hilo.0.03.lefse_summary
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
taxonomy_df <- read_taxonomy()
metadata <- read_metadata()
histology <- read_histology() %>% 
	mutate(hist_score = case_when(summary_score < 5 ~ 'Low',
		summary_score > 5 ~ 'High',
		T ~ 'NA'))

shared <- read_shared() %>% 
	select(-label, -numOtus) %>% 
	pivot_longer(-Group, names_to = 'OTU', values_to = 'counts') %>% 
	mutate(relative_abundance = counts/2107 * 100)	

lefse_design <- data.table::fread(here('data/process/lefse/temporal_trend_hilo.design'))
	colnames(lefse_design) <- c('Group', 'class')
    
lefse_df <- data.table::fread(here('data/process/lefse/temporal_trend_hilo.0.03.lefse_summary'),
		fill = T) %>% 
	filter(!is.na(LDA),
    LDA > 3) %>% 
	rename(tax_otu_label = OTU) %>% 
	mutate(fctr_class = factor(Class),
		order = LDA + 10 * as.numeric(fctr_class),
    tax_otu_label = gsub('_', ' ', tax_otu_label),
    tax_otu_label = gsub(' cla', '_cla', tax_otu_label),
		taxa_label = gsub('_unclassified', '*', tax_otu_label),
		taxa_label = gsub('_', ' ', taxa_label),
		taxa_label = paste0('*', taxa_label),
    taxa_label = gsub(' \\(', '*  
      \\(', taxa_label))
###############################################################################
#  setup data
###############################################################################
temporal_lefse_df <- taxonomy_df %>% 
  inner_join(lefse_df, by = 'tax_otu_label') %>% 
  left_join(shared, by = 'OTU') %>% 
  left_join(metadata, by = c('Group' = 'group')) %>%
  group_by(mouse_id, day, taxa_label, early_euth) %>% 
  summarise(relative_abundance = sum(relative_abundance)) %>% 
  inner_join(select(histology, mouse_id, hist_score), by = 'mouse_id') %>% 
  mutate(outcome = ifelse(early_euth, 'Severe', hist_score),
         outcome = factor(outcome, levels = c('Low', 'High', 'Severe'))) %>% 
  filter(outcome != 'NA') %>% 
  group_by(taxa_label) %>% 
  mutate(median_ra = median(relative_abundance),
    present = sum(relative_abundance > 0),
  	relative_abundance = ifelse(relative_abundance == 0, 0.045, 
  		relative_abundance)) %>% 
  filter(present > 86)

###############################################################################
#  plot data
###############################################################################
#lefse_ temporal trend hi/low day 0-10
temporal_lefse_plot <- temporal_lefse_df %>% 
  ggplot(aes(x = day, y = relative_abundance, group = interaction(day, outcome), 
  		color = outcome)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
                 position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 1/21.07, linetype = 'dashed', size = 0.25) + 
    facet_wrap(taxa_label~., scales = 'free_y', ncol = 4) + 
    scale_y_log10() + scale_x_continuous(breaks = c(0:10)) + 
    theme_bw() + 
    theme(strip.text = ggtext::element_markdown(),
    	legend.position = 'top',
    	strip.background = element_rect(fill = 'white')) + 
    scale_color_manual(breaks = c('Low', 'High', 'Severe'),
    	values = c("#62ebc9", "#399283", "#0c1b37"),
    	labels = c('Low histopatholic score',  'High histopatholic score', 'Severe')) + 
    labs(x = 'Day', y = 'Relative Abundance (%)', color = NULL)
###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_S4.tiff'),
       temporal_lefse_plot,
       height = 9, width = 6, unit = 'in',
  compression = 'lzw')
###############################################################################
