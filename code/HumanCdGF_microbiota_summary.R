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
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header = T, row.names=1)
tax_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)

#Create df with relative abundances
rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))

#Create vector of OTUs with median abundances >1%
med_rel_abund_cage <- aggregate(rel_abund, by=list(meta_file$cage_id),median)
cage_IDs <- med_rel_abund_cage[,"Group.1"]
med_rel_abund_cage <- med_rel_abund_cage[,!colnames(med_rel_abund_cage) %in% 'Group.1']
OTUs_1 <- apply(med_rel_abund_cage, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]

#df of OTUs with abundances >1% - by cage and inoculum
rel_abund_cage <- rel_abund[!meta_file$cage_IDs=='inoculum',OTUs_1]
inoculum_rel_abund <- rel_abund[meta_file$cage_id=='inoculum',OTUs_1]

#get taxonomy - df with columns for each level - 
# level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom)
## - convert taxonomy text list to df, remove percentages
## - subset df by desired tax level, then replace any unclassified with next level up
get_tax <- function(tax_level=1){
    if (tax_level %in% c(1:5)){
      taxonomy <- tax_file[OTU_list,]
      taxonomy <-  data.frame(do.call('rbind', strsplit(as.character(taxonomy$Taxonomy),';',fixed=TRUE)))
      taxonomy <- data.frame(sapply(taxonomy,gsub,pattern="\\(\\d*\\)",replacement=""))
              level <- 7-tax_level
      tax_out <- as.character(taxonomy[,level])
      for (i in level:2){
        next_level <- i-1
        tax_out[tax_out=='unclassified'] <- 
          as.character(taxonomy[tax_out=='unclassified',next_level])
      }
      return(data.frame(tax= tax_out, row.names=OTU_list))
    } else {print(
      'Error: Level must be 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum)')
    }
}

taxonomy <- get_tax(1)

#Create median rel abundance and IQRs
median_rel_abund <- aggregate()

low_QTR

up_QTR