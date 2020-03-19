#######
#
# Idenitify all available full Cdiff complete genomes
#
#######

library(tidyverse)

# read in table with all available C difficile genomic sequences
cdiff_genome_df <- read_delim('data/process/cdiff_16S/assembly_summary.txt', delim = '\t',
	skip = 1, comment = '') 

# filter to include only complete genomes
cdiff_genome_paths <- cdiff_genome_df %>% 
	filter(assembly_level %in% c('Chromosome', 'Complete Genome')) %>% 
	# convert to use rsync instead of ftp and select just the genomic.fna.gz
	mutate(rsync_path = paste0(gsub('ftp:', 'rsync:', ftp_path), 
		'/', `# assembly_accession`, '_', asm_name, '_genomic.fna.gz')) %>% 
	pull(rsync_path)

write.table(cdiff_genome_paths, file = "data/process/cdiff_16S/cdiff_genome_paths.txt", 
	quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)