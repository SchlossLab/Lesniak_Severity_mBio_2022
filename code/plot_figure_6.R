###############################################################################
#
# Figure 6
#	main idea - Is our observation conserved across isolates?
#   plot - CFU and histopathology score with different isolates
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
source(here::here('code/utilities.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata() %>% 
  mutate(isolate_431 = ifelse(cdiff_strain == 431, '431', 'Other'),
         isolate = factor(cdiff_strain, levels = c('431', '299', '395','458')))
  
toxin <- read_toxin()

histology <- read_histology() %>% 
  select(summary_score, mouse_id)

donor_aes <- donor_df

###############################################################################
#  setup data
###############################################################################
non_431_donors <- metadata %>% 
  filter(cdiff_strain != 431) %>% 
  pull(human_source) %>% 
  unique

cfu_data <- metadata %>% 
  filter(human_source %in% non_431_donors,
         day == 10) %>% 
  add_row(human_source = 'DA00578', cdiff_cfu = NA, isolate_431 = 'Other') %>% 
  left_join(toxin, by = 'group') %>% 
  left_join(donor_aes, by = c('human_source')) %>% 
  mutate(cdiff_cfu = ifelse(cdiff_cfu == 0, log10(60), log10(cdiff_cfu)),
         donor_labels = factor(donor_labels, levels = c('G' ,'F', 'M')))

histology_data <- metadata %>% 
  filter(human_source %in% non_431_donors) %>% 
  left_join(histology, by = 'mouse_id') %>% 
  select(cdiff_strain, summary_score, human_source, 
		mouse_id, isolate_431, isolate) %>% 
  unique %>% 
  left_join(donor_aes, by = 'human_source') %>% 
  mutate(donor_labels = factor(donor_labels, levels = c('G' ,'F', 'M')))

###############################################################################
#  plot data
###############################################################################
cfu_plot <- cfu_data %>% 
  ggplot(aes(x = donor_labels, y = cdiff_cfu, color = isolate_431)) + 
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) + 
    geom_hline(yintercept = 2, linetype = 'dashed', color = 'lightgray') + 
    geom_label(data = data.frame(x = 3, y = 2), aes(x = x, y = y), label = "LOD", 
    	fill = 'white', color = 'white', inherit.aes = F) + 
    geom_text(data = data.frame(x = 3, y = 2), aes(x = x, y = y), label = "LOD", 
    	color = 'lightgray', size = 12/.pt, inherit.aes = F) + 
    geom_text(data = data.frame(x = 2.75, y = 1.65), aes(x = x, y = y), 
    	label = "DECD", color = "#3588d1", size = 8/.pt, inherit.aes = F) + 
    geom_text(data = data.frame(x = 3.25, y = 1.65), aes(x = x, y = y), 
    	label = "DECD", color = "#76406b", size = 8/.pt, inherit.aes = F) + 
    scale_color_manual(values = c("#3588d1", "#76406b")) + 
    scale_y_continuous(limits = c(1.5, 8.5),
                       breaks = c(2,4,6,8), 
                       labels = c('10^2', '10^4', '10^6', '10^8')) + 
    theme_bw() + 
    theme(legend.position = 'right',
          axis.text.y = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          legend.title = ggtext::element_markdown(),
          plot.title = ggtext::element_markdown(hjust = 0.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) + 
    labs(x = 'Donor Source', y = '*C. difficile* CFU', 
    	title = 'Day 10', color = '*C. difficile*  
isolate')

histology_score_plot <- histology_data %>% 
  ggplot(aes(x = isolate_431, y = summary_score, 
             color = isolate_431, fill = isolate_431)) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5) + 
    facet_grid(.~donor_labels, switch = 'both') + 
    theme_bw() + 
    scale_color_manual(values = c("#3588d1", "#76406b")) + 
    scale_fill_manual(values = c("#3588d1", "#76406b")) + 
    scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    labs(x = 'Donor Source', y = 'Histopathologic score', title = 'Endpoint')

###############################################################################
#  plot data
###############################################################################
ggsave(here('results/figures/Figure_6.jpg'),
       cowplot::plot_grid(cfu_plot, histology_score_plot, 
                          nrow = 1, rel_widths = c(6,4)),
              height = 3, width = 6, unit = 'in')
###############################################################################
