###################
#
# make_Humanhuman_GF_mouse_file.R
#
#Create file with experiment data with toxin and metadata
#
#    Need files:
#         data/raw/humanGF_ids.xlsx
#         data/mothur/gf_all.an.0.03.subsample.shared
#
#    Outputs files:
#         data/process/human_CdGF_metadata.txt
<<<<<<< HEAD
#         data/process/human_CdGF.an.unique_list.0.03.subsample.shared
=======
#         data/process/human_CdGF_all.0.03.subsample.shared
>>>>>>> master
#
###################
input_mouse <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/raw/humanGF_ids.xlsx'
input_shared <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/mothur/gf_all.an.0.03.subsample.shared'
output_metadata <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/process/human_CdGF_metadata.txt'
output_shared <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/process/human_CdGF.an.unique_list.0.03.subsample.shared'

# files file used for metadata file
input_meta_sample_ids <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/raw/gf_cdiff.files'
# files file used for mothur processing
input_shared_sample_ids <- '~/Documents/Github/Schubert_humanCdGF_XXXX_2016/data/raw/gf_all.files'

# gdata loads excel files, 
pack_used <- c('gdata', 'dplyr', 'stringr')
for (dep in pack_used){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
                     quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

#read in data files
human_GF_mouse <- read.xls(input_mouse, sheet = 'complete metadata')
shared_file <- read.table(input_shared, sep = '\t',header = T)
meta_sample_ids <- read.table(input_meta_sample_ids, sep = '\t', header = F, stringsAsFactors = FALSE)
shared_sample_ids <- read.table(input_shared_sample_ids, sep = '\t', header = F, stringsAsFactors = FALSE)

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

# subset df to include only sample_id and otus
shared_file <- select(shared_file, sample_id, contains('Otu0'))

# Remove AMS and DA00810 human sources, cecum and day -7 samples, and file col
human_GF_mouse <- subset(human_GF_mouse, !human_source %in% c('AMS', 'DA00810') &
                                   sample_type == 'stool' &
                                   !human_GF_mouse$day %in% -7
                                 ,select = -file)

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
shared_file <- shared_file[shared_file$sample_id %in% human_GF_mouse$sample_id, ]
human_GF_mouse <- human_GF_mouse[human_GF_mouse$sample_id %in% shared_file$sample_id, ]

# add log transformed CFU
human_GF_mouse$log_cfu <- log10(human_GF_mouse$cdiff_cfu + 1)

# output metadata file
write.table(human_GF_mouse, file = output_metadata, quote = F, sep = '\t', row.names = F)
write.table(shared_file, file = output_shared, quote = F, sep = '\t', row.names = F)
