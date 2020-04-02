###################
#
# make_Humanhuman_GF_mouse_file.R
#
#Create file with experiment data with toxin and metadata
#
#	Need files:
#		data/raw/humanGF_ids.xlsx
#		data/process/sample.final.0.03.subsample.shared
#		data/process/sample.final.shared
#		data/raw/meta_humanGF_cdiff.files
#		data/raw/humanGF_cdiff.files
#
#	Outputs files:
#		data/process/human_CdGF_metadata.txt
#		data/process/human_CdGF.an.unique_list.0.03.subsample.shared
#
###################
pack_used <- c('tidyverse', 'readxl')
for (dep in pack_used){
	if (dep %in% installed.packages()[,"Package"] == FALSE){
		install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
					quiet=TRUE);
	}
	library(dep, verbose=FALSE, character.only=TRUE)
}

input_mouse <- 'data/raw/humanGF_ids.xlsx'
input_shared <- 'data/mothur/sample.final.0.03.subsample.shared'
output_metadata <- 'data/process/human_CdGF_metadata.txt'
output_shared <- 'data/process/human_CdGF.subsample.shared'
input_full_shared <- 'data/mothur/sample.final.shared'
output_full_shared <- 'data/process/humanCdGF.full.shared'

# files file used for metadata file
input_meta_sample_ids <- 'data/raw/meta_humanGF_cdiff.files' 
# files file used for mothur processing
input_shared_sample_ids <- 'data/raw/humanGF_cdiff.files'

#read in data files
human_GF_mouse <- read_xlsx(input_mouse, sheet = 'complete metadata')
shared_file <- read_tsv(input_shared)
meta_sample_ids <- read_tsv(input_meta_sample_ids, col_names = FALSE)
shared_sample_ids <- read_tsv(input_shared_sample_ids, col_names = FALSE)
shared_full <- read_tsv(input_full_shared)

# convert sample ids of shared to match metafile
colnames(meta_sample_ids) <- c('meta_sample_id','fastq_1', 'fastq_2')
colnames(shared_sample_ids) <- c('shared_sample_id','fastq_1', 'fastq_2')
sample_id_key <- merge(meta_sample_ids, shared_sample_ids, 
					   by= 'fastq_1', all = T)[, c('meta_sample_id', 'shared_sample_id')]
# add labels to human source samples
sample_id_key[is.na(sample_id_key$meta_sample_id), 'meta_sample_id'] <- 
  paste('DA', str_pad(sample_id_key[is.na(sample_id_key$meta_sample_id), 
									'shared_sample_id'], 5, pad = '0'), sep = '')
# change labels of group/sample_ids in shared to match metafile
shared_file <- merge(sample_id_key, shared_file, by.x = 'shared_sample_id', 
						by.y = 'Group', all.y = TRUE)
shared_file <- rename(shared_file, sample_id = meta_sample_id)
human_GF_mouse <- rename(human_GF_mouse, sample_id = group)

# do the above for the full shared file too 
shared_full <- merge(sample_id_key, shared_full, by.x = 'shared_sample_id', by.y = 'Group', all.y = TRUE)
shared_full <- rename(shared_full, sample_id = meta_sample_id)

# subset df to include only sample_id and otus
shared_file <- select(shared_file, sample_id, starts_with('Otu'))
shared_full <- select(shared_full, sample_id, starts_with('Otu'))

# Remove AMS human sources, cecum and day -7 samples, and file col
human_GF_mouse <- human_GF_mouse %>% 
	filter(!human_source %in% c('AMS'),
		!sample_id %in% c('581-inoculum', 'DA581'),
		!human_GF_mouse$day %in% -7) %>% 
	select(-file)

# add unique mouse identify column to differentiate NT mice
human_GF_mouse <- rename(human_GF_mouse, ear_tag = mouse_id)
human_GF_mouse$mouse_id <- paste(human_GF_mouse$cage_id, human_GF_mouse$ear_tag,
								 sep = '_')

# add T/F col for mice euthanized early
Euthanized <- human_GF_mouse %>% 
  select(mouse_id, day) %>% 
  group_by(mouse_id) %>% 
  summarise(last_day_by_mouse = max(day))
Euthanized$Early_Euth <- ifelse(Euthanized$last_day_by_mouse < 10, T, F)
human_GF_mouse <- merge(human_GF_mouse, Euthanized, by='mouse_id')

# include only samples with both metadata and OTU data
shared_full <- shared_full[shared_full$sample_id %in% human_GF_mouse$sample_id, ]
shared_file <- shared_file[shared_file$sample_id %in% human_GF_mouse$sample_id, ]

human_GF_mouse <- human_GF_mouse[human_GF_mouse$sample_id %in% shared_file$sample_id, ]
human_GF_mouse <- human_GF_mouse[human_GF_mouse$sample_id %in% shared_full$sample_id, ]

# add log transformed CFU
human_GF_mouse <- human_GF_mouse %>% 
	mutate(log_cfu = log10(as.numeric(cdiff_cfu) + 1))

# output metadata file
write_tsv(human_GF_mouse, path = output_metadata)
write_tsv(shared_file, path = output_shared)
write_tsv(shared_full, path = output_full_shared)

