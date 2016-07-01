###################
#
# make_HumanCdGF_metadata_file.R
#
#Create file with experiment data with toxin and metadata
#
#    Need files:
#         data/raw/humanGF_ids.xlsx
#         data/raw/Alyx_Humice_toxinassay_results.xlsx
#         data/raw/MIMARKS_cdclinical.xlsx
#         data/raw/humanGF_ids.xlsx
#         data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#
#    Outputs files:
#         data/process/human_CdGF_metadata.txt
#         data/process/human_CdGF.an.unique_list.0.03.subsample.shared
#
###################
setwd("~/Documents/Github/Schubert_humanCdGF_XXXX_2016")

#need gdata to read excel files
install.packages('gdata')
library(gdata)

#read in data files
human_GF_mouse <- read.xls('data/raw/humanGF_ids.xlsx',sheet='complete metadata')
human_GF_toxin <- read.xls('data/raw/Alyx_Humice_toxinassay_results.xlsx', sheet = 'to metadata')
human_GF_clinical <- read.xls('data/raw/MIMARKS_cdclinical.xlsx',sheet='Sheet1')
shared_file <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared', sep='\t',header = T)

#subset df to include only sample and otus
shared_file <- shared_file[,!colnames(shared_file) %in% c("label",'numOtus')]

#Subset the human source metadata dataframe to a data frame of only samples used in study
#   Subset sources to only ones used in study
human_GF_clinical <- human_GF_clinical[human_GF_clinical$sample_id %in% human_GF_mouse$human_source,]
#   Subset data columns to remove columns with same information
human_GF_clinical <- human_GF_clinical[,sapply(human_GF_clinical, function(col) length(unique(col)))>1]
#   Remove information from clinical seq and redundant columns
human_GF_clinical <- human_GF_clinical[,!colnames(human_GF_clinical) %in% 
                                         c("feature",'lib_reads_seqd','sfffile_id','mid')]
# Remove file names
human_GF_mouse <- human_GF_mouse[,!colnames(human_GF_mouse) %in% 'file']
# Remove AMS samples
human_GF_mouse <- human_GF_mouse[!human_GF_mouse$human_source=='AMS',]

#Combine Toxin, Exp, and metadata
human_CdGF_metadata <- merge(human_GF_mouse, human_GF_toxin,by.x='group',by.y='Cage_Mouse',
                               incomparables=NA, all.x=T)
human_CdGF_metadata <- merge(human_CdGF_metadata,human_GF_clinical,by.x='human_source',by.y='sample_id',
                           invcomparables=NA,all.x=T)

#include only samples with both metadata and OTU data
shared_file <- shared_file[shared_file$Group %in% human_CdGF_metadata$group,]
human_CdGF_metadata <- human_CdGF_metadata[human_CdGF_metadata$group %in% shared_file$Group,]
human_CdGF_metadata <- human_CdGF_metadata[order(human_CdGF_metadata$group),]

# - add mouse name column to differentiate NT mice
human_CdGF_metadata <- cbind(human_CdGF_metadata,Mouse_ID=paste(human_CdGF_metadata$mouse_id,human_CdGF_metadata$cage_id))

# - fill in all rows with end day data
#         create data frame with day euthanized and mouse name
Euthanized <- unique(human_CdGF_metadata[!is.na(human_CdGF_metadata$day_end),c('day_end','Mouse_ID')])
Euthanized <- rbind(Euthanized,unique(human_CdGF_metadata[!human_CdGF_metadata$Mouse_ID %in% Euthanized$Mouse_ID,c('day_end','Mouse_ID')]))
#         add day mouse NT LINE was euthanized, missing day but has all samples upto end of exp at day 10
Euthanized[Euthanized$Mouse_ID=='NT LINE','day_end'] <- 10
# add day euthanized data to metadata
human_CdGF_metadata <- merge(human_CdGF_metadata, Euthanized, by='Mouse_ID')
names(human_CdGF_metadata)[names(human_CdGF_metadata) %in% c('mouse_id','day_end.x','day_end.y')] <- c('mouse_eartag','day_end','Euthanized')

#output metadata file
write.table(human_CdGF_metadata, file='data/process/human_CdGF_metadata.txt',quote=F,sep='\t',row.names=F)
write.table(shared_file, file='data/process/human_CdGF.an.unique_list.0.03.subsample.shared',quote=F,sep='\t',row.names=F)
