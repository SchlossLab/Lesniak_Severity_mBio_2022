############
#
# determine variation 16S V4 of C difficile
#
############

library(tidyverse)
library(insect) # use rc() to generate reverse complement of 16S V4
library(msa) # for sequence alignment - msa() and sequence identity ?pid()


sequence_df <- read_tsv('data/process/cdiff_16S/all_cdifficile_16S_V4.txt')
 
genome_n <- length(unique(sequence_df$cdiff_genome)) # number of C difficile genomes analyzed
n_16S <- nrow(sequence_df) # number of 16S V4 regions present in all the genomes
V4_distribution_plot <- sequence_df %>% # distribution of 16S genes per genome
	count(cdiff_genome) %>% 
	ggplot(aes(n)) + geom_histogram() + 
		labs(x = 'Number of 16S rRNA genes', y = 'Count', 
			title = 'Number of 16S rRNA genes per C. difficile genome') +
		theme_bw()
# convert reversed sequences to forward and trim primers
cleaned_seqs <- sequence_df %>% 
	mutate(reverse_seq = rc(sequence_16S)) %>% 
	gather(direction, sequence, sequence_16S, reverse_seq) %>% 
	filter(grepl('^GTG', sequence)) %>% 
	mutate(sequence = gsub('^GTGCCAGCAGCCGCGGTAA|ATTAGATACCCTGGTAGTCC$', '', sequence)) %>% 
	select(-direction)

## create fasta of 16S V4
## convert dataframe columns of sequence, name to long format of sequence below sequence name
#fasta_frmt_data <- as.vector(t(data.frame(
#	cleaned_seqs$cdiff_genome, cleaned_seqs$sequence)))
#write_lines(fasta_frmt_data, 'data/process/cdiff_16S/cdifficile_16S_V4.fasta')

unique_seqs <- unique(cleaned_seqs$sequence)
unique_seqs_aligned <- msa(unique_seqs, type = 'dna')
# for alignment and logo plot use following code
# msaPrettyPrint(unique_seqs_aligned, output = 'pdf', showNames ='left')
# determine the adundance of each sequence in the genomes
sequence_distribution <- cleaned_seqs %>% 
	count(sequence) %>% 
	mutate(percent_of_16S_genes = paste(round(n/n_16S*100, 1), '%'))

# get  sequence in most (96.7%) genomes
concensus_seq <- sequence_distribution %>% 
	filter(percent_of_16S_genes == max(percent_of_16S_genes)) %>% 
	pull(sequence)
# compare each other sequence to the consensus sequence
sequence_identity_df <- map_dfr(unique_seqs, function(i){ 
		p_aln <- pairwiseAlignment(concensus_seq, i)
		percent_alignment <- pid(p_aln)	
		return(tibble(sequence = i, percent_identity = percent_alignment))
	})

write_lines(concensus_seq, 'data/process/cdiff_16S/consensus_cdifficile_16S_V4.txt')

# select sequences below 97% similar to consensus sequence
outlier_sequence_16S <- sequence_identity_df %>% 
	filter(percent_identity < 97) %>% 
	pull(sequence)
# compare sequence identity to determine if above 97% similar to any other sequence
outlier_seq_identity_df <-  map_dfr(unique_seqs, function(i){ 
		p_aln <- pairwiseAlignment(outlier_sequence_16S, i)
		percent_alignment <- pid(p_aln)	
		return(tibble(sequence = i, percent_identity = percent_alignment))
	})
# get the 16S V4 seqeunce that doesn't cluster with the other sequences
outlier_seq_identity_df %>% 
	filter(sequence != unique_sequence_16S) %>% 
	filter(percent_identity > 97) %>% 
	nrow
# get the genome that 16S came from
outlier_16S_genome <- cleaned_seqs %>% 
	filter(grepl(outlier_sequence_16S, sequence_16S)) %>% 
	pull(cdiff_genome)
# are there other 16S sequences in this genome?
cleaned_seqs %>% 
	filter(grepl(outlier_16S_genome, cdiff_genome)) %>% 
	nrow
# this genome has 15 concensus sequences and on non-consensus sequence
# and according to https://www.frontiersin.org/articles/10.3389/fmicb.2018.02994/full
# Strain M68 is ribotype 017
# for more information on the C difficile genome containing the outlier 16S sequence, use following code
# outlier_C_difficile <- read_delim('data/process/cdiff_16S/assembly_summary.txt', delim = '\t',
#		skip = 1, comment = '')  %>% 
#	filter(organism_name == gsub('.*\\d |, .*', '', outlier_16S_genome)) %>% 
#	data.frame

