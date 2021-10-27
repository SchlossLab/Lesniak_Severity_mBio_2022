###############################################################################
#
# Figure 1
#	main idea - germ0free mice colonized with human feces produce diverse communities
# 	plot - relative abundance of communities for each sample and 
#		beta diversity to compare distance from mice recieving same donor to others
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
source(here::here('code/utilities.R'))

###############################################################################
# laod data
###############################################################################
metadata <- read_metadata()

taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
  select(-label, -numOtus) %>% 
  pivot_longer(-Group, names_to = 'otu', values_to = 'counts') %>% 
  mutate(relative_abundance = counts/2107 * 100)
beta_df <- read_beta()

###############################################################################
#  setup data
###############################################################################
day_0_samples <- metadata %>% 
  filter(day == 0) %>% 
  select(group, human_source)

corecipient <- day_0_samples %>% 
  left_join(day_0_samples, by = c('human_source')) %>% 
  filter(group.x != group.y) %>% 
  inner_join(beta_df, by = c('group.x' = 'rows', 'group.y' = 'columns')) %>% 
  mutate(comparison = 'Recipients compared to\nothers w/same donor')

other_recipients <- beta_df %>% 
  filter(grepl('D0', rows) & grepl('D0', columns)) %>% 
  anti_join(corecipient, by = c('rows' = 'group.x', 'columns' = 'group.y')) %>% 
  mutate(comparison = 'Recipients compared to\nothers w/different donor')

###############################################################################
#  analyze data
###############################################################################
wilcox.test(corecipient$distances, other_recipients$distances)

###############################################################################
#  plot data
###############################################################################
day_0_plot <- metadata %>% 
  filter(day == 0,
         cdiff_strain == 431) %>% 
  select(group, human_source, mouse_id) %>% 
  left_join(filter(shared, counts > 0), by = c('group' = 'Group')) %>% 
  left_join(taxonomy, by = c('otu' = 'OTU')) %>% 
  rename(taxa_level = Class) %>% 
  group_by(group, taxa_level, human_source) %>% 
  summarise(relative_abundance = log10(sum(relative_abundance) + 0.0001)) %>% 
  group_by(taxa_level) %>% 
  filter(!is.na(taxa_level)) %>% 
  mutate(median_ra = mean(relative_abundance),
         human_source = as.character(as.numeric(gsub('DA', '', human_source))),
         human_source = factor(human_source, levels =
                                 c('10093', '10148', '431', '884', '1245',
                                   '581', '1324', '953', '1134',
                                   '369', '430', '10027', '10034',
                                   '10082', '578')),
         taxa_level = gsub('_unclassified', '', taxa_level),
         taxa_level = gsub('_', ' ', taxa_level),
         taxa_level = paste0('*', taxa_level, '*')) %>% 
  ggplot(aes(x = group, 
             y = reorder(taxa_level, median_ra), 
             fill = relative_abundance)) + 
    geom_tile() + 
    scale_fill_gradient2(low="white", mid='#0B775E', high = 'black',
		  	limits = c(-2,2), na.value = NA, midpoint = .3,
			  breaks = c(-2.5, -1, 0, 1, 2), labels = c('', '0.1', '1', '10', '100')) + 
	  theme_bw() + 
    labs(x = NULL, y = NULL, fill = NULL) + 
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    facet_grid(.~human_source, scales = 'free_x') + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'bottom',
          axis.text.y = ggtext::element_markdown(angle = 45),
          strip.text.x = element_text(angle = 45),
          strip.background = element_blank(),
          legend.margin=margin(t=-0.1, r=0, b=-0.1, l=0, unit="cm"),
          legend.key.height = unit(0.3, 'cm'))

refined_beta_plot <- bind_rows(corecipient, other_recipients) %>% 
  mutate(comparison = ifelse(comparison == 'Recipients compared to\nothers w/different donor',
                             'Different donor', 'Same donor'),
    comparison = factor(comparison, levels = 
                              c('Same donor',
                                'Different donor'))) %>% 
  ggplot(aes(x = comparison, y = distances)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    #geom_jitter(alpha = 0.1, width = 0.3) + 
    labs(x = NULL, y = 'Similarity') +
    scale_y_continuous(breaks = c(0, 1), limits = c(0,1),
                       labels = c('Complete\nsimilarity', 'No similarity')) + 
    theme_bw()

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_1.jpg'),
  cowplot::plot_grid(day_0_plot + theme(text=element_text(size = 9)), 
                     NULL,
                     beta_plot + theme(text=element_text(size = 9)), 
                     ncol = 1, 
                     rel_heights = c(6, 0.3, 3),
                     labels = c('A', 'B', NULL)),
  height = 9.0625, width = 6.875, unit = 'in')
###############################################################################
