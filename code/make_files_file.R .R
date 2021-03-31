# make_files_file.R
# create the files file for mothur using the names of fastq sequences and the meta data

library(tidyverse)

human_GF_mouse <- read_tsv('data/process/metadata_tidy.tsv',
	col_type = cols(.default = col_character())) # load metadata
fastq_list <- read_csv('data/process/list_seqs.txt', col_names = F, col_type = 'c') # read in the names of all fastq available from project
#fastq_list <- list.files(path = 'data/raw/', pattern = 'fastq') # read in the names of all fastq available from project

# create vector with all cage names used in experiment
cages <- human_GF_mouse %>% 
	pull(cage_id) %>% 
	unique

# create vector with all human source numbers used in experiment
sources <- read_tsv('data/process/human_source_tidy.tsv',
	col_type = cols(.default = col_character())) %>% 
	pull(sample_id)
source_number_meta <- as.numeric(gsub('DA0*','',sources))

# separate fastqs by sample type - mock, inoculum, experiment since each has a generally unique name structure
fastq_df <- fastq_list %>% 
	rename(fastq = X1) %>% 
	mutate(sample_names = gsub('_S.*', '', fastq), # remove sequencing file naming suffixes
		sample_type = case_when(grepl('[Mm]ock', sample_names) ~ 'mock',
			grepl('DA', sample_names) ~ 'inoculum', # inocula all have a DA number, so select most with this
			!grepl('D', sample_names) ~ 'inoculum', # remaining samples without a D will be inocula since all experiment samples have day denoted with a D
			T ~ 'experiment')) 

# subset fastqs to the cages/sources that are used in the analysis
# subset by cage
experiment_fastqs <- fastq_df %>% 
	filter(sample_type == 'experiment') %>% 
	separate(sample_names, c('cage', 'day'), 
		sep = '[^MID](?<!^)D', remove = F) %>% # ignore cage MID, separate at any D that doesnt start the string (which should split only on the day and not any cage name)
	mutate(cage = gsub('-.*', '', cage)) %>% # remove mouse tag from sample name to leave cage name
	filter(cage %in% cages) %>% # select only samples from cages used in the analysis/metadata
	select(sample_names, fastq) %>% 
	filter(!grepl('-D-7', sample_names)) %>% # remove day -7 fastq samples
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')
# subset inoculum
inoculum_fastqs <- fastq_df %>% 
	filter(sample_type == 'inoculum',
		!is.na(as.double(sample_names))) %>% 
	filter(sample_names %in% source_number_meta) %>% # select only samples from sources used on analysis
	# convert to files format - sample, forward read, reverse read columns
	select(sample_names, fastq) %>% 
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')
# subset mock
mock_fastqs <- fastq_df %>% 
	filter(sample_type == 'mock') %>% 
	mutate(sample_names = gsub('.*(S\\d{2,3}).*', 'mock_\\1', fastq)) %>% 
	select(sample_names, fastq) %>% 
	filter(sample_names %in% c('mock_S282', 'mock_S343', 'mock_S344', 'mock_S102', 'mock_S96')) %>% # remove water samples labeled as mock 
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')

files_file <- bind_rows(experiment_fastqs, inoculum_fastqs, mock_fastqs) %>% 
	# since there is - in the sample nanmes and this is the symbol mothur uses to split sames
	# mothur attempts to split all the file names containing - 
	# replace all - with another character ('v')
	mutate(sample_names = gsub('-', '_', sample_names))

write_tsv(files_file, 'data/raw/humanGF_cdiff.files', col_names = FALSE)

# fastqs to remove (unused donors, mock, water, days)
unused_fastqs <- fastq_list  %>% 
	anti_join(pivot_longer(files_file, cols = c('R1', 'R2')), by = c('X1' = 'value'))

write_tsv(unused_fastqs, 'data/process/unused_fastqs.txt', col_names = F)