library(tidyverse)
library(here)


source(here('code/Sum_OTU_by_Tax.R'))

# load in data
metadata <- read_tsv(here('data/process/metadata_tidy.tsv'),
                     col_type = 'cccccDddddDcdddddcdl') %>% 
  mutate(isolate_431 = ifelse(cdiff_strain == 431, '431', 'Other'),
         isolate = factor(cdiff_strain, levels = c('431', '299', '395','458')))
  
toxin <- read_tsv(here('data/process/toxin_tidy.tsv'),
                  col_type = 'cccccdc')
histology <- read_tsv(here('data/process/histology_tidy.tsv'),
                      col_type = 'ccccdddcdcc') %>% 
  select(summary_score, mouse_id)

donor_colors <- c("#41bbc5", "#7e2640", "#9ad859", "#682dbd", "#728e24", "#e26df8", "#23980d", "#fd048f", "#46e33e", "#fa1bfc", "#155126", "#f9b2d1", "#154975", "#e2c627", "#6294ce")


non_431_donors <- metadata %>% 
  filter(cdiff_strain != 431) %>% 
  pull(human_source) %>% 
  unique

cfu_plot <- metadata %>% 
  filter(human_source %in% non_431_donors,
         day == 10) %>% 
  add_row(human_source = 'DA00578', cdiff_cfu = NA, isolate_431 = 'Other') %>% 
  left_join(toxin, by = 'group') %>% 
  mutate(cdiff_cfu = ifelse(cdiff_cfu == 0, log10(60), log10(cdiff_cfu)),
         donor_source = case_when(human_source == 'DA00369' ~ 'F',
                                  human_source == 'DA00430' ~ 'G',
                                  human_source == 'DA00578' ~ 'M'),
         donor_source = factor(donor_source, levels = c('G' ,'F', 'M'))) %>% 
  ggplot(aes(x = donor_source, y = cdiff_cfu, color = isolate_431)) + 
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) + 
    geom_hline(yintercept = 2, linetype = 'dashed', color = 'lightgray') + 
    geom_label(x = 3, y = 2, label = "LOD", fill = 'white', color = 'white') + 
    geom_text(x = 3, y = 2, label = "LOD", color = 'lightgray', size = 12/.pt) + 
    geom_text(x = 2.75, y = 1.65, label = "DECD", color = "#3588d1", size = 10/.pt) + 
    geom_text(x = 3.25, y = 1.65, label = "DECD", color = "#76406b", size = 10/.pt) + 
    scale_color_manual(values = c("#3588d1", "#76406b")) + 
    scale_y_continuous(limits = c(1.5, 8.5),
                       breaks = c(2,4,6,8), 
                       labels = c('10^2', '10^4', '10^6', '10^8')) + 
    theme_bw() + 
    theme(legend.position = 'right',
          axis.text.y = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          plot.title = ggtext::element_markdown(hjust = 0.5),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) + 
    labs(x = 'Donor Source', y = NULL, title = 'Day 10 *C. difficile* CFU', color = NULL)

#metadata %>% 
#  filter(human_source %in% non_431_donors) %>% 
#  left_join(toxin, by = 'group') %>% 
#  ggplot(aes(x = day, y = Log_repiricoal_dilution, color = cdiff_strain == 431)) + 
#  geom_point() + 
#  facet_wrap(human_source~.) + 
#  theme_bw()

clinical_score_plot <- metadata %>% 
  filter(human_source %in% non_431_donors) %>% 
  left_join(histology, by = 'mouse_id') %>% 
  select(cdiff_strain, summary_score, human_source, mouse_id, isolate_431, isolate) %>% 
  unique %>% 
  mutate(donor_source = case_when(human_source == 'DA00369' ~ 'F',
                                  human_source == 'DA00430' ~ 'G',
                                  human_source == 'DA00578' ~ 'M'),
         donor_source = factor(donor_source, levels = c('G' ,'F', 'M'))) %>% 
  ggplot(aes(x = isolate_431, y = summary_score, 
             color = isolate_431, fill = isolate_431)) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5) + 
    facet_grid(.~donor_source, switch = 'both') + 
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
    labs(x = 'Donor Source', y = NULL, title = 'Endpoint Clinical Score')

ggsave(here('results/figures/Figure_6.jpg'),
       cowplot::plot_grid(cfu_plot, clinical_score_plot, 
                          nrow = 1),
              height = 4, width = 8, unit = 'in')
