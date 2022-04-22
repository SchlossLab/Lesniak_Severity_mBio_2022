###############################################################################
#
# Figure 1
#	main idea - germ-free mice colonized with human feces produce diverse communities
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

donor_aes <- donor_df

###############################################################################
#  setup data
###############################################################################
day_0_samples <- metadata %>% 
  filter(day == 0) %>% 
  select(group, human_source)

donor <- bind_rows(
  day_0_samples %>% 
    inner_join(beta_df, by = c('group' = 'columns', 'human_source' = 'rows')) %>% 
    mutate(comparison = 'Donor to\nrecipient'),
  day_0_samples %>% 
    inner_join(beta_df, by = c('group' = 'rows', 'human_source' = 'columns')) %>% 
    mutate(comparison = 'Donor to\nrecipient'))

corecipient <- day_0_samples %>% 
  left_join(day_0_samples, by = c('human_source')) %>% 
  filter(group.x != group.y) %>% 
  inner_join(beta_df, by = c('group.x' = 'rows', 'group.y' = 'columns')) %>% 
  mutate(comparison = 'Recipients of\nsame donor')

other_recipients <- beta_df %>% 
  filter(grepl('D0', rows) & grepl('D0', columns)) %>% 
  anti_join(corecipient, by = c('rows' = 'group.x', 'columns' = 'group.y')) %>% 
  mutate(comparison = 'Recipients of\ndifferent donor')

beta_data <- bind_rows(donor, corecipient, other_recipients) %>% 
  mutate(comparison = factor(comparison, levels = 
                              c('Recipients of\ndifferent donor', 
                                'Recipients of\nsame donor',
                                'Donor to\nrecipient')))

relative_abundance_data <- metadata %>% 
  filter(day == 0,
         cdiff_strain == 431) %>% 
  select(group, human_source, mouse_id) %>% 
  left_join(filter(shared, counts > 0), by = c('group' = 'Group')) %>% 
  left_join(taxonomy, by = c('otu' = 'OTU')) %>% 
  rename(taxa_level = Class) %>% 
  group_by(group, taxa_level, human_source, mouse_id) %>% 
  summarise(relative_abundance = log10(sum(relative_abundance) + 0.0001)) %>% 
  group_by(taxa_level) %>% 
  filter(!is.na(taxa_level)) %>% 
  left_join(donor_aes, by = c('human_source')) %>% 
  mutate(median_ra = mean(relative_abundance),
         taxa_level = gsub('_', ' ', taxa_level),
         taxa_level = paste0('*', taxa_level, '*'),
         taxa_level = ifelse(grepl('unclassified', taxa_level), 
         		paste('Unclassified', gsub(' unclassified\\*', '*', taxa_level)),
         		taxa_level),
         mouse_id = ifelse(mouse_id %in% c('NP1_2597', 'NP1_2725'), '+', '')) %>% 
  filter(median_ra > -0.677)


###############################################################################
#  analyze data
###############################################################################
#wilcox.test(corecipient$distances, other_recipients$distances)
# p < 2.2e-16
#wilcox.test(corecipient$distances, donor$distances)
# p < 2.2e-16
#wilcox.test(donor$distances, other_recipients$distances)
# p = 0.004092
###############################################################################
#  plot data
###############################################################################
day_0_plot <- relative_abundance_data %>% # only plot top 10
  ggplot(aes(x = group, 
             y = reorder(taxa_level, median_ra), 
             fill = relative_abundance)) + 
    geom_tile() + 
    scale_fill_gradient2(low="white", mid='#0B775E', high = 'black',
		  	limits = c(-2,2), na.value = NA, midpoint = .3,
			  breaks = c(-2.5, -1, 0, 1, 2), labels = c('', '0.1', '1', '10', '100')) + 
    scale_x_discrete(labels = c(rep('', 8), rep('+', 2), rep('', 41))) + 
	  theme_bw() + 
    labs(x = NULL, y = NULL, fill = NULL) + 
    facet_grid(.~donor_labels, scales = 'free_x', space = 'free_x') + 
    theme(axis.ticks.x = element_blank(),
          legend.position = 'bottom',
          axis.text.x = element_text(margin = margin(t = -3)),
          axis.text.y = ggtext::element_markdown(),
          strip.background = element_blank(),
          legend.margin=margin(t=-0.1, r=0, b=-0.1, l=0, unit="cm"),
          legend.key.height = unit(0.3, 'cm'))

beta_plot <- beta_data %>% 
  ggplot(aes(x = comparison, y = distances)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
    #geom_jitter(alpha = 0.05, width = 0.3) + 
    labs(x = NULL, y = expression(theta['YC'])) +
    coord_flip() + 
    theme_bw()

###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_1.tiff'),
  cowplot::plot_grid(cowplot::plot_grid(
                      day_0_plot + theme(text=element_text(size = 9)),
                      labels = c('       Donor'), label_size = 9, label_fontface = 'plain'), 
                     NULL,
                     beta_plot + theme(text=element_text(size = 9)), 
                     ncol = 1, 
                     rel_heights = c(6, 0.4, 2.2),
                     labels = c('A', 'B', NULL)),
  height = 5.5, width = 5.1, unit = 'in',
  compression = 'lzw')
###############################################################################
