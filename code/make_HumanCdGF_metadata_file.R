###################
#
# make_HumanCdGF_metadata_file.R
#
#Create file with experiment data with toxin and metadata
#
#
#
###################
setwd("~/Documents/SchlossLab_NL/git-Projects/Schubert_humanCdGF_XXXX_2016/code")

human_GF_exp <- read.table('../data/raw/humanGF_ids.txt', sep='\t',header = T)
human_GF_tox <- read.table('../data/raw/Alyx_Humice_toxinassay_results.txt',sep='\t',header=T)
human_GF_IDs <- read.table('../data/raw/MIMARKS_cdclinical.txt',sep='\t',header=T)

#Create list of human sources used in experiment
Sources <- levels(human_GF_exp$human_source)

#Subset the human source metadata dataframe to a data frame of only samples used in study
human_source_data <- human_GF_IDs[human_GF_IDs$sample_id %in% Sources,c(1,8:9,13,17,32:45)]
#Combine Toxin, Exp, and metadata
human_GF_exp_metadata <- merge(human_GF_exp, human_GF_tox,by.x='group',by.y='Cage_Mouse',
                               incomparables=NA, all.x=T)
human_GF_metadata <- merge(human_GF_exp_metadata,human_source_data,by.x='human_source',by.y='sample_id',
                           invcomparables=NA,all.x=T)
#output metadata file
write.table(human_GF_metadata, file='../data/raw/human_GF_metadata.txt',quote=F,sep='\t')
