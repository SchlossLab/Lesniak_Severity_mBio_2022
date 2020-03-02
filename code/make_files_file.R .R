# make_files_file.R
# create the files file for mothur using the names of fastq sequences and the meta data

library(tidyverse)
library(readxl)

input_mouse <- 'data/raw/humanGF_ids.xlsx'
human_GF_mouse <- read_excel(input_mouse, sheet = 'complete metadata') # load metadata
fastq_list <- list.files(path = 'data/mothur/', pattern = 'fastq') # read in the names of all fastq available from project

# create vector with all cage names used in experiment
cages <- human_GF_mouse$cage_id%>%unique%>%as.character
cages <- cages[!cages %in% c('N', 'NC', 'Y')]

# create vector with all human source numbers used in experiment
sources <- human_GF_mouse %>% 
	filter(cage_id != 'inoculum') %>% 
	filter(human_source != 'AMS') %>% 
	pull(human_source) %>% unique%>%as.character
source_number_meta <- as.numeric(gsub('DA0*','',sources))

# separate fastqs by sample type - mock, inoculum, experiment since each has a generally unique name structure
fastq_df <- data.frame(fastq = fastq_list, stringsAsFactors = F) %>% 
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
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')
# subset inoculum
inoculum_fastqs <- fastq_df %>% 
	filter(sample_type == 'inoculum') %>% 
	filter(!(sample_names == '581-incolum' | sample_names == 'DA581')) %>% # remove extra 581 samples from separate run than rest of inocula samples
	mutate(source_number = as.numeric(regmatches(sample_names, gregexpr('\\d{3,}', sample_names)))) %>% # remove any extra text/symbols other than the source number
	filter(source_number %in% source_number_meta) %>% # select only samples from sources used on analysis
	# convert to files format - sample, forward read, reverse read columns
	select(sample_names, fastq) %>% 
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')
# subset mock
mock_fastqs <- fastq_df %>% 
	filter(sample_type == 'mock') %>% 
	mutate(sample_names = gsub('.*(S\\d{2,3}).*', 'mock_\\1', fastq)) %>% 
	select(sample_names, fastq) %>% 
	mutate(read = str_extract(fastq, 'R[12]')) %>%
	spread(read, fastq, fill = 'NA')

files_file <- bind_rows(experiment_fastqs, inoculum_fastqs, mock_fastqs)	

write.table(files_file, 'data/mothur/humanGF_cdiff.files', sep = '\t',
	quote = FALSE, row.names = FALSE, col.names = FALSE)