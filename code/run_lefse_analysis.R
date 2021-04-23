

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

run_lefse <- function(sample_df_name, tax_level){
	i <- sample_df_name
	current_df <- get(i)

	current_shared <- shared %>% 
		filter(Group %in% current_df$Group) %>% 
		pivot_longer(cols = -c('label', 'Group', 'numOtus'),
			names_to = "OTU", values_to = 'value') %>% 
		left_join(taxonomy, by = c('OTU')) %>% 
		group_by(label, Group, numOtus, .data[[tax_level]]) %>% 
		summarise(value = sum(value)) %>% 
		pivot_wider(names_from = .data[[tax_level]], values_from = value) %>% 
		ungroup
	# remove otus that either have 0 or only present in 5 or fewer samples
	present_otus <- current_shared %>% 
		select(-label, -Group, -numOtus) %>% 
		map_dbl(~ sum(. > 0)) %>% 
		which(x = (. > 5)) %>% 
		names
	current_shared <- current_shared %>% 
		select(label, Group, numOtus, one_of(present_otus)) %>% 
		mutate(numOtus = length(present_otus))

	# write files to be used in mothur for lefse analysis
	write_tsv(path = paste0('data/process/', i, '_', tax_level, '.shared'), 
		x = filter(current_shared, Group %in% current_df$Group))
	write_tsv(path = paste0('data/process/', i, '_', tax_level, '.design'), 
		x = filter(current_df, Group %in% current_shared$Group))
	print(tax_level)
	# run lefse
	system(paste0('code/mothur/mothur "#set.dir(input=data/process, output=data/process);
		lefse(shared=', i, '_', tax_level, '.shared, design=', i, '_', tax_level, '.design)"'))
}

# lefse - severe disease by OTU
severe_disease <- metadata %>% 
	filter(cdiff_strain == 431,
		day == 0) %>% 
	select(Group = group, early_euth)

run_lefse('severe_disease', 'OTU')
run_lefse('severe_disease', 'Genus')

