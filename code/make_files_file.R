################################################################################
#
# make_files_file.R
#
# This script creates the files file for mothur 
# 	using the names of fastq sequences and the meta data
#
# Dependencies...
# * data/raw/*fastq
# * data/process/human_source_tidy.tsv
#
# Output...
# * data/raw/humanGF_cdiff.files
#
################################################################################

# load packages
library(tidyverse)

# read in the names of all fastq available from project
fastq_list <- list.files(path = 'data/raw/', pattern = 'fastq')

# separate fastqs by sample type - mock, inoculum, experiment since each has a different name structure
fastq_df <- tibble(fastq = fastq_list) %>% 
	mutate(sample_names = gsub('_S.*', '', fastq), # remove sequencing file naming suffixes
		sample_type = case_when(grepl('[Mm]ock', sample_names) ~ 'mock',
			grepl('D', sample_names) ~ 'experiment', # all experiment samples have day denoted with a D
			T ~ 'inoculum')) # remaining samples will be inocula 

# format experiment samples into files format
experiment_fastqs <- fastq_df %>% 
	filter(sample_type %in% 'experiment') %>% 
	select(sample_names, fastq) %>% 
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')
# format inoculum samples into files format and add S## suffix to differentiate mocks
mock_fastqs <- fastq_df %>% 
	filter(sample_names %in% c('mock', 'mock1', 'mock2', 'Mock1')) %>% 
	mutate(sample_names = gsub('.*(S\\d{2,3}).*', 'mock_\\1', fastq)) %>% 
	select(sample_names, fastq) %>% 
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')

# replace all - with _  since mothur uses - to split sames
files_file <- bind_rows(experiment_fastqs, mock_fastqs) %>% 
	mutate(sample_names = gsub('-', '_', sample_names))

# output files file
write_tsv(files_file, 'data/mothur/humanGF_cdiff.files', col_names = FALSE)