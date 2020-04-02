#####################
#
#	Fig 1 - Compare Donor and Initial Communities
#
#	Need files
#		data/process	
#			sample_shared
#		data/mothur
#			donor_shared
#			beta_diversity
#			alpha_diversity
#
#####################
library(tidyverse)
library(cowplot)

donor_shared <- read_tsv('data/mothur/donor/donor.final.0.03.subsample.shared')
donor_taxonomy <- read_tsv('data/process/donor.final.taxonomy.clean.tsv')
donor_beta <- read_tsv('data/mothur/donor/donor.final.thetayc.0.03.lt.ave.nmds.axes')

humanGF_shared <- read_tsv('data/process/human_CdGF.subsample.shared')%>% 
	filter(grepl('D0', sample_id)) %>% 
	rename(Group=sample_id)
humanGF_taxonomy <- read_tsv('data/process/sample.final.taxonomy.clean.tsv')
humanGF_beta <- read_tsv('data/mothur/sample.final.thetayc.0.03.lt.ave.nmds.axes')


source('code/Sum_OTU_by_Tax.R')
# function sums OTUs by taxonomic level
# sum_otu_by_taxa(taxonomy_df, otu_df, taxa_level, top_n=0, silent=T)

donor_abundance_plot <- sum_otu_by_taxa(donor_taxonomy, donor_shared, 'Genus', top_n=20)  %>% 
	group_by(Group) %>% 
	mutate(total = sum(abundance),
		relative_abundance = abundance/total * 100) %>% 
	ggplot(aes(y = Group, x =taxa, fill = relative_abundance)) + 
		geom_tile() +
		scale_fill_gradient(low="white", high="black", limits = c(0,100)) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
		labs(x = NULL, y = NULL, title = 'Donor Community') #+ 
		#guides(fill='none')

day0_abundance_plot <- sum_otu_by_taxa(humanGF_taxonomy, humanGF_shared, 'Genus', top_n=20)  %>% 
	group_by(Group) %>% 
	mutate(total = sum(abundance),
		relative_abundance = abundance/total * 100) %>% 
	ggplot(aes(y = Group, x =taxa, fill = relative_abundance)) + 
		geom_tile() +
		scale_fill_gradient(low="white", high="black", limits = c(0,100)) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
		labs(x = NULL, y = NULL, title = 'Murine Community - Day 0 (~2 weeks post inoculation)')
ggsave('~/Desktop/ra_donor.jpg', donor_abundance_plot)
ggsave('~/Desktop/ra_humanGF.jpg', day0_abundance_plot)


donor_beta %>% 
	ggplot(aes(x = axis1, y = axis2, color = group)) +
	geom_point() + 
	theme_bw()

humanGF_beta %>% 
	ggplot(aes(x = axis1, y = axis2, color = group)) +
	geom_point() + 
	theme_bw()


plot_grid(donor_abundance_plot, day0_abundance_plot)