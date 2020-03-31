library(tidyverse)
library(RColorBrewer)
library(cowplot)

# get error rate of mock communities and plot them

get_error_rate <- function(data_set=T){

	if(data_set == 'sample'){
		seq_length <- 253
	} else if(data_set == 'donor'){
		seq_length <- 257
	} else {
		stop('Input either "donor" or "sample" to select a dataset')
	}

	error_df <- read.table(paste0('data/process/', data_set, '.error.count'), header = T)
	total_seqs <- sum(error_df$Sequences)
	total_errors <- sum(error_df$Errors * error_df$Sequences)
	error_rate <- total_errors / (total_seqs * seq_length) * 100

	# get mock community abundances
	mock_file <- read.table(paste0("data/process/mock.", data_set, ".final.shared"), sep='\t', header = TRUE, stringsAsFactors = F) %>% 
		mutate(Group = gsub('(MOCK|mockv)', '', Group))
	
	if(data_set == 'sample'){
		mock_file <- filter(mock_file, Group %in% c('S282','S343','S344','S96')) 
	} 	

	tax_file <- read.table(paste0('data/process/', data_set, '.final.taxonomy.clean.tsv'), sep='\t',header = T, stringsAsFactors = F)
	
	# group OTUs by genus
	source('code/Sum_OTU_by_Tax.R')
	taxonomy_genus <- sum_otu_by_taxa(tax_file, mock_file, taxa_level = 'Genus')
	# convert data to long format
	mock_file <- mock_file %>% 
		select(-numOtus, -label) %>% 
		gather(otu, abundance, -Group)
	# subset the OTUs to just the most abundant
	top_otus <- mock_file %>% 
		group_by(otu) %>% 
		summarise(total_abundance = sum(abundance)) %>% 
		top_n(20, total_abundance) %>% 
		pull(otu)
	# determine the total abundance for each sample
	total_abundance <- mock_file %>% 
		group_by(Group) %>% 
		summarise(total_abundance = sum(abundance)) 

	# Species in the HMP mock
	# Acinetobacter baumannii, Actinomyces odontolyticus, Bacillus cereus, Bacteroides vulgatus, Candida albicans, Clostridium beijerinckii, Deinococcus radiodurans, Enterococcus faecalis, Escherichia coli, Helicobacter pylori, Lactobacillus gasseri, Listeria monocytogenes, Methanobrevibacter smithii, Neisseria meningitidis, Propionibacterium acnes, Pseudomonas aeruginosa, Rhodobacter sphaeroides, Staphylococcus aureus, Staphylococcus epidermidis, Streptococcus agalactiae, Streptococcus mutans, Streptococcus pneumoniae
	HMP_mock <- c('Acinetobacter', 'Actinomyces', 'Bacillus', 'Bacteroides', 'Candida', 'Clostridium', 'Deinococcus', 'Enterococcus', 'Escherichia', 'Helicobacter', 'Lactobacillus', 'Listeria', 'Methanobrevibacter', 'Neisseria', 'Propionibacterium', 'Pseudomonas', 'Rhodobacter', 'Staphylococcus', 'Staphylococcus', 'Streptococcus', 'Streptococcus', 'Streptococcus')

	brewer_colors <- c(brewer.pal(n = 12, name = "Paired"), brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Dark2"))

	plot_by_otu <- mock_file %>% 
		mutate(otu = case_when(otu %in% top_otus ~ otu,
			T ~ 'Other')) %>% 
		group_by(Group, otu) %>% 
		summarise(absolute = sum(abundance)) %>% 
		full_join(total_abundance, by = c('Group')) %>% 
		mutate(relative = absolute/total_abundance * 100) %>% 	
		ggplot(aes(x = Group, y = relative, fill = otu)) + 
			geom_bar(stat = 'identity', position = 'stack') +
			scale_fill_manual(values = brewer_colors[1:(length(top_otus)+1)]) +
			theme_bw() + 
			labs(x = NULL, y = 'Relative Abundance', title = 'Top 20 OTUs in Mock',
				subtitle = paste0('All others grouped in Other. Error Rate = ', 
					round(error_rate, 2), '%')) + 
			guides(fill = 'none')

	plot_by_genus <- taxonomy_genus %>% 
		mutate(taxonomy = case_when(taxa %in% HMP_mock ~ taxa,
			T ~ 'Other')) %>% 
		group_by(Group, taxonomy) %>% 
		summarise(absolute = sum(abundance)) %>% 
		full_join(total_abundance, by = c('Group')) %>% 
		mutate(relative = absolute/total_abundance * 100) %>% 	
		ggplot(aes(x = Group, y = relative, fill = taxonomy)) + 
			geom_bar(stat = 'identity', position = 'stack') +
			scale_fill_manual(values = c(brewer.pal(n = 8, name = "Paired"), brewer.pal(n = 8, name = "Set1"))) +
			theme_bw() + 
			labs(x = NULL, y = NULL, title = 'Genus level grouping of Species in Mock',
				subtitle = 'Color legend only applies to genus plot')

	ggsave(paste0('exploratory/', data_set, '_mock_by_genus.jpg'), 
		plot_grid(plot_by_otu, plot_by_genus), width = 12, height = 7)
}

get_error_rate('sample')
get_error_rate('donor')