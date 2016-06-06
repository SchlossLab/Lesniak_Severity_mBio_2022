###################
#
# HumanCdGF_microbiota_summary.R
#
#Create file with microbiota community summary data
#
#
#
###################
setwd("~/Documents/Github/Schubert_humanCdGF_XXXX_2016")

#read in files
meta_file <- read.table('data/process/human_GF_metadata.txt', sep='\t',header = T)
shared_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header = T, row.names = 2)
tax_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)


#Eliminate columns with label and numOtus
shared_file <- shared_file[,!colnames(shared_file) %in% c("label",'numOtus')]
#Create vector with number of seq per sample
num_seq <- apply(shared_file, 1, sum)
#Create df with relative abundances
rel_abund <- 100*shared_file/num_seq

#remove samples not followed in this study or sequenced
rel_abund <- rel_abund[rownames(rel_abund) %in% meta_file$group,]
meta_file <- meta_file[meta_file$group %in% rownames(rel_abund),]
rownames(meta_file) <- meta_file[,2]
meta_file <- meta_file[order(rownames(meta_file)),-2]

#determine median abundances >1%
med_rel_abund_cage <- aggregate(rel_abund, by=list(meta_file$cage_id),median)
cage_IDs <- med_rel_abund_cage[,1]
med_rel_abund_cage <- med_rel_abund_cage[,-1]
otus_1 <- apply(med_rel_abund_cage, 2, max) > 1
med_rel_abund_cage_1 <- med_rel_abund_cage[!cage_IDs=='inoculum',otus_1]
num_otus <- sum(otus_1)
#abundances of inoculum in OTUs >1%
inoculum_rel_abund <- rel_abund[meta_file$cage_id=='inoculum',otus_1]

#get taxonomy
taxonomy <- tax_file$Taxonomy
  taxonomy <- gsub("\\(\\d*\\)",'',taxonomy)#elimnate values
  taxonomy <- gsub(";unclassified", "", taxonomy) #eliminate unclassified
  taxonomy <- gsub("/.*", "", taxonomy) #eliminate 
  taxonomy <- gsub(";$", "", taxonomy) #eliminate end ;
  taxonomy <- gsub("\"$", "", taxonomy) #eliminate end parenthesis
  taxonomy <- gsub(".*;", "", taxonomy) #eliminate preceding text, except last
  taxonomy <- gsub("^\"", "", taxonomy) #eliminate beginning parenthesis
  names(taxonomy) <- rownames(tax_file)
  taxonomy <- taxonomy[colnames(med_rel_abund_cage_1)] #only want phyla for OTUs >1%

