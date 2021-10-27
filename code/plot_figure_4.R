###############################################################################
#
# Figure 4
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

taxonomy <- read_taxonomy()

shared <- read_shared()

###############################################################################
#  create plot lefse data function
###############################################################################
plot_lefse <- function(file_name){
  lefse_design <- data.table::fread(here(paste0('data/process/lefse/', file_name, '.design')))
      colnames(lefse_design) <- c('Group', 'class')
    
  lefse_df <- data.table::fread(here(paste0('data/process/lefse/', file_name, '.0.03.lefse_summary')),
    	fill = T) %>% 
    filter(!is.na(LDA)) %>% 
    top_n(20, LDA) %>% 
    mutate(fctr_class = factor(Class),
      order = LDA + 10 * as.numeric(fctr_class),
    	tax_otu_label = OTU)
    
   shared %>% 
    gather(OTU, abundance, matches('Otu\\d*')) %>% 
  	left_join(taxonomy, by = 'OTU') %>% 
  	right_join(lefse_df, by = c('Genus' = 'OTU')) %>% 
    right_join(lefse_design, by = c('Group')) %>% 
  	group_by(Group, Genus, order, class) %>% # removed mouse_id from interactive plot
  	summarise(abundance = sum(abundance)) %>% 
  	mutate(tax_otu_label = Genus) %>% 
  	mutate(abundance = abundance / 21.07,
  		abundance = ifelse(abundance == 0, 0.045, abundance)) %>% 
    #highlight_key(~mouse_id) %>% 
  	ggplot(aes(x = tax_otu_label, 
  			y = abundance, color = class)) + # , shape = mouse_id remove mouse_id from interactive plot
  		geom_hline(yintercept = 0.047, linetype = 'dashed', lwd = 0.3) + 
      stat_summary(fun.data = 'median_hilow', aes(group = class),
  			fun.args = (conf.int=0.5), position = position_dodge(width = .8)) +
  		theme_bw() + 
      labs(x = NULL, y = 'Relative Abundance', color = NULL, title = file_name) + 
  		scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
  			labels=c("0","0.1","1","10","100")) +
  		#scale_color_manual(values = c('#8de4d3','#f75b5f')) + 
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_line(color = 'black'),
  			legend.position = 'top') + 
  		coord_flip() + 
    geom_jitter(position = position_jitterdodge(dodge.width = .8, jitter.width = 0.1),
                alpha = 0.1)
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
    left_join(shared, by = c('Group')) %>% 
    group_by(otu) %>%
    mutate(sd = sd(counts)) %>% 
    filter(!is.na(otu),
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
    left_join(taxonomy, by = c('otu' = 'OTU')) %>% 
    mutate(tax_otu_label = gsub(' \\(', '\n\\(', tax_otu_label))

###############################################################################
#  plot data
###############################################################################
day_10_hist_genus_corr_plot <- day_10_hist_genus_corr %>% 
    ggplot(aes(x = relative_abundance + 0.001, y = test, color = tax_otu_label)) + 
      geom_point(alpha = 0.3) +
      geom_smooth(method = 'lm', se = F) +
      scale_x_log10() + 
      facet_wrap(tax_otu_label ~ ., scales = 'free_x') + 
      theme_bw() + 
      labs(x = 'Relative Abundance', y = 'Clinical Score') + 
      guides(color = 'none')

day_0_toxin_production_plot <- plot_lefse('day_0_toxin_production_Genus')
day_0_severity_Genus_plot <- plot_lefse('day_0_severity_Genus')

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_4.jpg'),
  cowplot::plot_grid(day_0_severity_Genus_plot + theme(text=element_text(size = 9)),
                     cowplot::plot_grid(day_0_toxin_production_plot + theme(text=element_text(size = 9)),
                                        day_10_hist_genus_corr_plot + theme(text=element_text(size = 9)),
                                        ncol = 1, rel_heights = c(3,2),
                                        labels = c('B', 'C')),
                     nrow = 1, 
                     labels = c('A', NULL)),
  height = 7, width = 6.875, unit = 'in')

###############################################################################
