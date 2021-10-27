###############################################################################
#
# Figure 
# 	plot -
#	main idea - 
#
# need files:
# 	data/process/
#	data/mothur/
#	data/process/
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
histology <- read_histology()
# outcome from run_lefse_temporal_trend.R
lefse_hilo <- c("Anaeromassilibacillus", "Anaeroplasma", 
	"Bacteroidales_unclassified", "Bifidobacterium", "Blautia", "Clostridioides", 
	"Clostridium_IV", "Clostridium_XlVa", "Duncaniella", "Escherichia/Shigella", 
	"Faecalibacterium", "Fournierella", "Helicobacter", "Hungatella", "Ihubacter", 
	"Intestinimonas", "Klebsiella", "Kosakonia", "Lactobacillales_unclassified", 
	"Legionella", "Muribaculaceae_unclassified", "Neglecta", "Odoribacter", 
	"Oscillibacter", "Paraclostridium", "Paramuribaculum", 
	"Prevotellaceae_unclassified", "Ruthenibacterium")

###############################################################################
#  setup data
###############################################################################
temporal_lefse_df <- taxonomy_df %>% 
  filter(Genus %in% lefse_hilo) %>% 
  left_join(pivot_longer(shared, 
                         cols = -Group, names_to = 'OTU', values_to = 'counts')) %>% 
  left_join(metadata, by = c('Group' = 'group')) %>%
  group_by(mouse_id, day, Genus, early_euth) %>% 
  summarise(counts = sum(counts)) %>% 
  inner_join(select(histology, mouse_id, hist_score), by = 'mouse_id') %>% 
  mutate(outcome = ifelse(early_euth, 'severe', hist_score),
         outcome = factor(outcome, levels = c('low', 'high', 'severe')),
         rel_abun = counts/2107 * 100) %>% 
  filter(outcome != 'mid') %>% 
  group_by(Genus) %>% 
  mutate(present = sum(counts > 0)) %>% 
  filter(present > 60)

###############################################################################
#  plot data
###############################################################################
#lefse_ temporal trend hi/low day 0-10
temporal_lefse_plot <- temporal_lefse_df %>% 
  ggplot(aes(x = day, y = rel_abun, group = interaction(day, outcome), color = outcome)) + 
    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
                 position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 1/21.07, linetype = 'dashed', size = 0.25) + 
    facet_wrap(Genus~., scales = 'free_y') + 
    scale_y_log10() + 
    theme_bw()
###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_XX_temporal_lefse.jpg'),
       temporal_lefse_plot,
       height = 9, width = 6)
###############################################################################
