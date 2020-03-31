
# get error rate
error_df <- read.table('data/mothur/humanGF_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.error.count', header = T)
total_seqs <- sum(error_df$Sequences)
total_errors <- sum(error_df$Errors * error_df$Sequences)
total_errors / (total_seqs*252) * 100

# get mock community abundances
mock_file <- read.table("data/process/mock.final.shared", sep='\t', header = TRUE, stringsAsFactors = F)
tax_file <- read.table('data/process/final.taxonomy.clean.tsv', sep='\t',header = T, stringsAsFactors = F)
source('code/Sum_OTU_by_Tax.R')

taxonomy_genus <- sum_otu_by_taxa(tax_file, mock_file, taxa_level = 'Genus')

mock_file %<>% 
	select(-numOtus, -label) %>% 
	gather(otu, abundance, -Group)
top_otus <- mock_file %>% 
	group_by(otu) %>% 
	summarise(total_abundance = sum(abundance)) %>% 
	top_n(17, total_abundance) %>% 
	pull(otu)
total_abundance <- mock_file %>% 
	group_by(Group) %>% 
	summarise(total_abundance = sum(abundance)) 

# Species in the mock
# Bacillus_subtilis, Enterococcus_faecalis, Escherichia_coli, Lactobacillus_fermentum, Listeria_monocytogenes, Pseudomonas_aeruginosa, Salmonella_enterica, Staphylococcus_aureus
zymo_mock <- c('Bacillus', 'Enterococcus', 'Escherichia/Shigella', 'Lactobacillus', 'Listeria', 'Pseudomonas', 'Salmonella', 'Staphylococcus')


plot_by_otu <- mock_file %>% 
	mutate(otu = case_when(otu %in% top_otus ~ otu,
		T ~ 'Other')) %>% 
	group_by(Group, otu) %>% 
	summarise(absolute = sum(abundance)) %>% 
	full_join(total_abundance, by = c('Group')) %>% 
	mutate(relative = absolute/total_abundance * 100) %>% 	
	gather(abundance, value, absolute, relative) %>% 
	ggplot(aes(x = Group, y = value, fill = otu)) + 
		geom_bar(stat = 'identity', position = 'stack') +
		scale_fill_manual(values = c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 9, name = "Set1"))) +
		facet_wrap(.~abundance, scales = 'free_y') + 
		theme_bw()
ggsave('~/Desktop/mock_by_otu.jpg', plot_by_otu, width = 10, height = 7)

plot_by_genus <- taxonomy_genus %>% 
	mutate(taxonomy = case_when(taxa %in% zymo_mock ~ taxa,
		T ~ 'Other')) %>% 
	group_by(Group, taxonomy) %>% 
	summarise(absolute = sum(abundance)) %>% 
	full_join(total_abundance, by = c('Group')) %>% 
	mutate(relative = absolute/total_abundance * 100) %>% 	
	gather(abundance, value, absolute, relative) %>% 
	ggplot(aes(x = Group, y = value, fill = taxonomy)) + 
		geom_bar(stat = 'identity', position = 'stack') +
		scale_fill_manual(values = brewer.pal(n = 8, name = "Paired")) +
		facet_wrap(.~abundance, scales = 'free_y') + 
		theme_bw() + 
		labs(title = 'Genus level mock')
ggsave('~/Desktop/mock_by_genus.jpg', plot_by_genus, width = 10, height = 7)

