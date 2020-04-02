###################
#
# sum_otu_by_taxa.R
#
# created by Nicholas Lesniak in 2016
# modified by Nicholas Lesniak in 2019
#    
# function outputs a dataframe with OTUs summed by taxa group per sample
# 
# taxonomy_df - dataframe formatted with convert_OTU_labels.R
# otu_df - dataframe with samples as rownames and otus in columns
# taxa_level - genus, family, order, class, phylum, kingdom
# 
#
###################

levels <- c('Kingdom','Phylum','Class','Order','Family','Genus', 'OTU', 'tax_otu_label', 'otu_label')

sum_otu_by_taxa <- function(taxonomy_df, otu_df, taxa_level = 'NA', top_n = 0, silent = T){
	
	if(!all(is.data.frame(taxonomy_df), is.data.frame(otu_df),
		taxa_level %in% levels)) {stop(paste0(
			'Check to make sure you have entered a data frame for taxonomy and shared 
			and you have selected a classification level - ',
			paste0(levels, collapse = ', ')))}
	if(length(taxonomy_df$OTU) < 1){stop(paste0(
			'Check to make sure you have formatted the taxonomy file with convert_OTU_labels.R'))
	}
	if(silent == F){
		print(paste0('Summing shared by ', taxa_level))
	}

	output_dataframe <- otu_df %>% 
		gather(OTU, abundance, starts_with('Otu')) %>% 
		inner_join(select(taxonomy_df, OTU, taxa = one_of(taxa_level)), by = 'OTU') %>% 
		group_by(Group, taxa) %>%
		summarize(abundance = sum(abundance))

	if(top_n > 0){
		if(silent == F){
			print(paste0('Returning the top ', top_n, ' groups, all others are summed in Other'))
		}

		top_taxa <- output_dataframe %>% 
			group_by(taxa) %>% 
			summarise(median_abundance = mean(abundance, na.rm = T)) %>% 
			top_n(top_n, median_abundance) %>% 
			select(taxa) %>% 
			mutate(top_taxa = taxa)

		output_dataframe <- output_dataframe %>%
			full_join(top_taxa, by = 'taxa') %>%  
			mutate(taxa = ifelse(is.na(top_taxa), 'Other', taxa)) %>% 
			group_by(Group, taxa) %>% 
			summarise(abundance = sum(abundance)) %>% 
			mutate(dataset = paste0('top_', top_n, '_by_', taxa_level)) %>% 
			ungroup 

		} else {
			if(silent == F){
				print('Returning all groups')
			}
		}

	return(output_dataframe)
	}####
# function to output a dataframe with OTU rel abund summed by tax group per sample
#
# input file with tax level class and file with samples OTU rel abund
####

sum_OTU_by_tax_level <- function(TAX_LEVEL, OTU_DF, tax_df){
		 source('code/tax_level.R')
		 TAX_DF <-get_tax(TAX_LEVEL,names(OTU_DF), tax_df)
		 tax_levels <- as.character(unique(TAX_DF$tax))
		 OUTPUT_DF <- data.frame(rep(0,length(rownames(OTU_DF))), row.names=rownames(OTU_DF))
		 for (i in 1:length(tax_levels)){
					OTU_by_level <- rownames(TAX_DF)[TAX_DF$tax %in% tax_levels[i]]
					if (length(OTU_by_level)>1){
							 level_column <- apply(OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]],1,sum)
					} else {
							 level_column <- OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]]
					}     
					OUTPUT_DF[,i] <- level_column
					colnames(OUTPUT_DF)[i] <- tax_levels[i]
		 }
		 return(OUTPUT_DF)
}
