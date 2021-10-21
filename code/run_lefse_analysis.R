

library(tidyverse)
library(here)

metadata <- read_tsv(here('data/process/metadata_tidy.tsv'),
		col_type = 'cccccDddddDcdddddcdl') 
shared <- read_tsv(here('data/mothur/sample.final.0.03.subsample.shared'),
		col_types = cols(.default = col_double(),
			Group = col_character())) %>% 
	mutate(Group = ifelse(grepl('_D|DA', Group), Group,
		gsub('D', '_D', Group)))
toxin <- read_tsv(here('data/process/toxin_tidy.tsv'),
                  col_type = 'cccccdc') %>% 
  filter(mouse_id %in% metadata$mouse_id)
histology <- read_tsv(here('data/process/histology_tidy.tsv'),
                      col_type = 'ccccdddcdcc') %>% 
	             mutate(mouse_id = paste0(cage_id, '_', ear_tag))
taxonomy <- read_tsv(here('data/process/final.taxonomy.tidy.tsv'),
                     col_types = cols(.default = col_character()))

run_lefse <- function(sample_df_name, tax_level){
	i <- sample_df_name
	current_df <- get(i)

	cfu_column <- metadata %>% 
		select(Group = group, cdiff_cfu) %>% 
		mutate(cdiff_cfu = ifelse(cdiff_cfu == 0, log10(60), log10(cdiff_cfu)),
			cdiff_cfu = round(2107 * ((cdiff_cfu - min(cdiff_cfu, na.rm = T)) / 
											(max(cdiff_cfu, na.rm = T) - min(cdiff_cfu, na.rm = T))), 0)) %>% 
		rename(C_difficile_CFU = cdiff_cfu)

	current_shared <- shared %>% 
		filter(Group %in% current_df$Group) %>% 
		pivot_longer(cols = -c('label', 'Group', 'numOtus'),
			names_to = "OTU", values_to = 'value') %>% 
		left_join(taxonomy, by = c('OTU')) %>% 
		group_by(label, Group, numOtus, .data[[tax_level]]) %>% 
		summarise(value = sum(value)) %>% 
		pivot_wider(names_from = .data[[tax_level]], values_from = value) %>% 
		ungroup %>% 
	# add cdifficile cfu
		left_join(cfu_column, by = c('Group'))

	# remove otus that either have 0 or only present in 5 or fewer samples
	present_otus <- current_shared %>% 
		select(-label, -Group, -numOtus) %>% 
		map_dbl(~ sum(. > 0, na.rm = T)) %>% 
		which(x = (. > 5)) %>% 
		names
	current_shared <- current_shared %>% 
		select(label, Group, numOtus, one_of(present_otus)) %>% 
		mutate(numOtus = length(present_otus))

	# write files to be used in mothur for lefse analysis
	write_tsv(file = paste0('data/process/lefse/', i, '_', tax_level, '.shared'), 
		x = filter(current_shared, Group %in% current_df$Group))
	write_tsv(file = paste0('data/process/lefse/', i, '_', tax_level, '.design'), 
		x = filter(current_df, Group %in% current_shared$Group))
	print(tax_level)
	# run lefse
	system(paste0('code/mothur/mothur "#set.dir(input=data/process/lefse, output=data/process/lefse);
		lefse(shared=', i, '_', tax_level, '.shared, design=', i, '_', tax_level, '.design)"'))
}

# lefse - severe disease by OTU
severe_disease <- metadata %>% 
	filter(cdiff_strain == 431,
		day == 0) %>% 
	select(Group = group, early_euth)

toxin_presence <- toxin %>% 
  mutate(toxin = Log_repiricoal_dilution >1) %>%
  right_join(filter(metadata, cdiff_strain == 431), 
             by = 'group', suffix = c('.toxin', '')) %>% 
  group_by(mouse_id) %>% 
  mutate(toxin_production = any(toxin, na.rm =T)) %>% 
  ungroup()

toxin_severity <- toxin_presence %>% 
  filter(day == 0) %>% 
  mutate(toxin_outcome = case_when(early_euth == T & toxin_production == T ~ 'toxin_moribund',
                                   early_euth == T & toxin_production == F ~ 'no_toxin_moribund',
                                   early_euth == F & toxin_production == T ~ 'toxin_mild',
                                   early_euth == F & toxin_production == F ~ 'no_toxin_mild')) %>% 
  select(Group = group, toxin_outcome) %>% 
  data.frame

day_0_toxin_production <- toxin %>%
  group_by(mouse_id) %>% 
  summarise(toxin_presence = max(Log_repiricoal_dilution) > 1) %>%
  right_join(filter(metadata, cdiff_strain == 431), 
             by = 'mouse_id', suffix = c('.toxin', '')) %>% 
  filter(day == 0,
         !is.na(toxin_presence)) %>% 
  select(Group = group, toxin_presence)

toxin_presence <- toxin_presence %>% 
  filter(toxin_sample_type == 'stool') %>% 
	select(Group = group, toxin)

summary_score_hilo_day10 <- metadata %>% 
	filter(cdiff_strain == 431,
		day == 10) %>% 
	inner_join(histology, by = c('mouse_id')) %>% 
	mutate(summary_score = case_when(summary_score > 5 ~ 'high',
									 summary_score < 5 ~ 'low',	
									 T ~ 'NA')) %>% 
	filter(summary_score != 'NA') %>% 
	select(Group = group, summary_score)

summary_score_hilo_day0 <- metadata %>% 
  filter(cdiff_strain == 431,
         day == 0) %>% 
  inner_join(histology, by = c('mouse_id')) %>% 
  mutate(summary_score = case_when(summary_score > 5 ~ 'high',
                                   summary_score < 5 ~ 'low',	
                                   T ~ 'NA')) %>% 
  filter(summary_score != 'NA') %>% 
  select(Group = group, summary_score)

day_0_severity <- metadata %>% 
  filter(cdiff_strain == 431,
         day == 0) %>% 
  inner_join(histology, by = c('mouse_id')) %>% 
  mutate(summary_score = case_when(early_euth == T ~ 'severe',
                                   summary_score > 5 ~ 'moderate',
                                   summary_score < 5 ~ 'mild',	
                                   T ~ 'NA')) %>% 
  filter(summary_score != 'NA') %>% 
  select(Group = group, summary_score)



run_lefse('severe_disease', 'OTU')
run_lefse('severe_disease', 'Genus')
run_lefse('toxin_presence', 'Genus')
run_lefse('toxin_severity', 'Genus')
run_lefse('summary_score_hilo_day10', 'Genus')
run_lefse('summary_score_hilo_day0', 'Genus')
run_lefse('day_0_severity', 'Genus')
run_lefse('day_0_toxin_production', 'Genus')