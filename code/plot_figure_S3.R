###############################################################################
#
# Figure 
#   plot -
# main idea - 
#
# need files:
#   data/process/
# data/mothur/
# data/process/
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
library(tidyverse)
library(here)
source(here('code/Sum_OTU_by_Tax.R'))

###############################################################################
# load and steup data
###############################################################################
metadata <- read_tsv(here('data/process/metadata_tidy.tsv'),
                     col_type = 'cccccDddddDcdddddcdl') %>% 
  filter(human_source %in% c('DA00430', 'DA00578', 'DA00369'),
         day == 0) %>% 
  mutate(human_source = as.character(as.numeric(gsub('DA', '', human_source))),
         human_source = factor(human_source, levels = c('430','369', '578')))
         
shared <- read_tsv(here('data/mothur/sample.final.0.03.subsample.shared'),
                   col_types = cols(.default = col_double(),
                                    Group = col_character())) %>% 
  mutate(Group = ifelse(grepl('_D|DA', Group), Group,
                        gsub('D', '_D', Group))) %>% 
  select(-label, -numOtus) %>% 
  pivot_longer(-Group, names_to = 'otu', values_to = 'counts') %>% 
  filter(counts > 0) %>% 
  mutate(relative_abundance = counts/2107 * 100)

taxonomy <- read_tsv(here('data/process/final.taxonomy.tidy.tsv'),
                     col_types = cols(.default = col_character()))

features_df <- read_tsv(here('data/process/ml/ml_feature_imp.tsv'),
                        col_type = 'ddcccddcc') %>% 
  filter(dataset != 'day_0_histology_50') %>% 
  mutate(dataset = gsub('_80', '', dataset))

day_0_moribund <- metadata %>% 
  select(group, mouse_id, day, early_euth) %>% 
  filter(day == 0) %>% 
  inner_join(shared, by = c('group' = 'Group')) %>% 
  mutate(early_euth = ifelse(early_euth, 'severe', 'mild'))

###############################################################################
# create relative abundance plot
###############################################################################
day_0_plot <- metadata %>% 
  filter(day == 0,
         cdiff_strain == 431) %>% 
  select(group, human_source, mouse_id) %>% 
  left_join(shared, by = c('group' = 'Group')) %>% 
  left_join(taxonomy, by = c('otu' = 'OTU')) %>% 
  rename(taxa_level = Genus) %>% 
  group_by(group, taxa_level, human_source) %>% 
  summarise(relative_abundance = log10(sum(relative_abundance) + 0.0001)) %>% 
  group_by(taxa_level) %>% 
  filter(!is.na(taxa_level)) %>% 
  mutate(taxa_level = gsub('_unclassified', '', taxa_level),
         taxa_level = gsub('_', ' ', taxa_level),
         taxa_level = paste0('*', taxa_level, '*'),
         present = length(relative_abundance) > 3) %>% 
  filter(present) %>% 
  ggplot(aes(x = group, 
             y = taxa_level, 
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

###############################################################################
# create genus rank RF plot
###############################################################################
day_0_moribund_plot_df <- features_df %>% 
  filter(dataset == 'day_0_moribund',
         method == 'rf',
         taxonomic_level == 'Genus') %>% 
  group_by(names) %>% 
  mutate(median_diff = median(perf_metric_diff)) %>% 
  ungroup() %>% 
  filter(median_diff > 0) 
day_0_moribund_feature_plot <- day_0_moribund_plot_df %>% 
  ggplot(aes(x = reorder(names,median_diff), y = perf_metric_diff)) + 
  stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
  #geom_jitter(alpha = 0.1, width = 0.1) + 
  coord_flip() + 
  facet_grid(method~dataset) + 
  labs(x = NULL, y = 'AUC Diff with feature permuted') +
  theme_bw()

day_0_moribund_features <- day_0_moribund_plot_df %>% 
  select(names, median_diff) %>% 
  unique %>% 
  arrange(median_diff) %>% 
  pull(names)

day_0_moribund_feature_abundance_plot <- day_0_moribund %>% 
  left_join(taxonomy, by = c('otu' = "OTU")) %>% 
  filter(Genus %in% day_0_moribund_features) %>% 
  group_by(group, early_euth, Genus) %>% 
  summarise(counts = sum(counts)/21.07) %>% 
  mutate(taxon = factor(Genus, day_0_moribund_features)) %>% 
  ggplot(aes(x = taxon, y = counts, color = early_euth)) + 
  stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
               position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3), 
              alpha = 0.2) + 
  scale_y_log10() + 
  theme_bw() + 
  coord_flip() + 
  labs(x = NULL, y = 'Relative Abundance (%)', color  = 'Disease Severity') + 
  theme(legend.position = 'bottom')

###############################################################################
# output figure
###############################################################################
ggsave(here('results/figures/Figure_S3.jpg'),
       cowplot::plot_grid(
         day_0_plot + theme(text=element_text(size = 7)),
         cowplot::plot_grid(
           day_0_moribund_feature_plot + theme(text=element_text(size = 7)),
           day_0_moribund_feature_abundance_plot + theme(text=element_text(size = 7)),
           ncol = 1),
         nrow = 1),
       height = 7, width = 7, unit = 'in')
###############################################################################