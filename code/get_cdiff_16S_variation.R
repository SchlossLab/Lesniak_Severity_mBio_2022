############
#
# determine variation 16S V4 of C difficile
#
############

library(tidyverse)
library(insect) # use rc() to generate reverse complement of 16S V4
library(msa)

sequence_df <- read_tsv('data/process/cdiff_16S/cdifficile_16S_V4.txt')
 
genome_n <- length(unique(sequence_df$cdiff_genome)) # number of C difficile genomes analyzed
n_16S <- nrow(sequence_df) # number of 16S V4 regions present in all the genomes
16S distribution_plot <- sequence_df %>% # distribution of 16S genes per genome
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

# create fasta of 16S V4
# convert dataframe columns of sequence, name to long format of sequence below sequence name
fasta_frmt_data <- as.vector(t(data.frame(
	cleaned_seqs$cdiff_genome, cleaned_seqs$sequence)))
write_lines(fasta_frmt_data, 'data/process/cdiff_16S/cdifficile_16S_V4.fasta')

unique_seqs <- unique(cleaned_seqs$sequence)
unique_seqs_aligned <- msa(unique_seqs, type = 'dna')
msaPrettyPrint(unique_seqs_aligned, output = 'pdf', showNames ='left')
# determine the adundance of each sequence in the genomes
sequence_distribution <- cleaned_seqs %>% 
	count(sequence) %>% 
	mutate(percent_of_16S_genes = paste(round(n/n_16S*100, 1), '%')) %>% 
	select(-n)
# get  seqeunce in most (96.7%) genomes
concensus_seq <- sequence_distribution %>% 
	filter(percent_of_16S_genes == max(percent_of_16S_genes)) %>% 
	pull(sequence)
# compare each other sequence to the consensus sequence
lapply(unique_seqs, function(i){
		p_aln <- pairwiseAlignment(concensus_seq, i)
		percent_alignment <- pid(p_aln)	
		return(c(i, percent_alignment))
	})

msa_unq_aln_nmd <- msa(inputSeqs = seqs_unique_aligned_named$sequence, 
	names = paste(seqs_unique_aligned_named$names, letters[1:11]), type = 'dna')
msaPrettyPrint(msa_unq_aln_nmd, output = 'pdf')

	mutate(sequence_set = str_extract(sequence, '^...'),
		length_16S = nchar(sequence))  %>% 
	filter(n < 500) %>% 
	left_join(sequence_df, by = c('sequence' = 'sequence_16S')) %>% 
	data.frame

library("ggseqlogo")
data(ggseqlogo_sample)
ggplot() + geom_logo( seqs_dna$MA0001.1 ) + theme_logo()

ggplot() + geom_logo( seqs_unique_aligned ) + theme_logo()



