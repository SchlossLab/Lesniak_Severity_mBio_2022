##GET NEW SHARED FILE !!!

##what bacteria were present in the mice on day 0?

setwd("~/Documents/Schloss_Lab/Schubert_humanCdGF_XXXX_2016")

#import Nick's awesome combined metadata file
meta_file <- read.table("data/process/human_CdGF_metadata.txt", header = TRUE, sep='\t', fill = TRUE, row.names=2)
shared_file <- read.table("data/process/human_CdGF.an.unique_list.0.03.subsample.shared", sep='\t', header = TRUE, row.names=1)
tax_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)

####Use Nick's awesome code to get the OTU abundance table sorted out right###

#Create df with relative abundances
rel_abund <- 100*shared_file/unique(apply(shared_file, 1, sum))

#Create vector of OTUs with median abundances >1%
med_rel_abund_cage <- aggregate(rel_abund, by=list(meta_file$cage_id),median)
cage_IDs <- med_rel_abund_cage[,"Group.1"]
med_rel_abund_cage <- med_rel_abund_cage[,!colnames(med_rel_abund_cage) %in% 'Group.1']
OTUs_1 <- apply(med_rel_abund_cage, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]

#df of OTUs with abundances >1% - by cage and inoculum
rel_abund_cage <- rel_abund[!meta_file$cage_id=='inoculum',OTUs_1]
inoculum_rel_abund <- rel_abund[meta_file$cage_id=='inoculum',OTUs_1]

#df of OTUs w abundances >1% 
rel_abund_d0 <- rel_abund[meta_file$day == 0, OTUs_1]
#it's adding weird NAs, idk why but remove them:
rel_d0 <- na.omit(rel_abund_d0)

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

taxonomy_genus <- get_tax(1)
taxonomy_phylum <- get_tax(5)

sum_OTU_by_tax_level <- function(TAX_DF,OTU_DF){
  tax_levels <- as.character(unique(TAX_DF$tax))
  OUTPUT_DF <- data.frame(rep(0,length(rownames(OTU_DF))), row.names=rownames(OTU_DF))
  for (i in 1:length(tax_levels)){
    OTU_by_level <- rownames(TAX_DF)[TAX_DF$tax %in% tax_levels[i]]
    if (length(OTU_by_level)>1){
      level_column <- apply(OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]],1,sum)
    } else {
      level_column <- OTU_DF[,names(OTU_DF)[names(OTU_DF) %in% OTU_by_level]]
    }     
    OUTPUT_DF[,i] <- level_column
    colnames(OUTPUT_DF)[i] <- tax_levels[i]
  }
  return(OUTPUT_DF)
}

#get phylum level for D0 only

phylum_d0 <- sum_OTU_by_tax_level(taxonomy_phylum, rel_d0)

#plot each mouse individually

colors <- c('red', 'blue', 'green', 'orange', 'purple', 'pink')
legend_labels <- as.character(taxonomy_phylum$tax)
legend_labels[5] <- "Other"

barplot(t(phylum_d0), ylab='Relative Abundance', main="Taxonomic composition of mice on D0", 
        col = colors, ylim=c(0,100), cex.lab=0.9, cex.axis=0.7, cex.names=0.7)

#combine by cages/donors/ median per cage 

#add a column for cages to the df
cages <- meta_file$cage_id[meta_file$day==0]
cages <- na.omit(cages)
rel_d0[90] <- cages
colnames(rel_d0)[90] <- "cage"

d0_med <- aggregate(rel_d0[, 1:89], list(rel_d0$cage), median)
d0_med_phy <- sum_OTU_by_tax_level(taxonomy_phylum, d0_med)
cage_only <- cage_IDs[-15]
d0_med_phy[7] <- cage_only

d0_med_phy <- sum_OTU_by_tax_level(taxonomy_phylum, d0_med)
barplot(t(d0_med_phy), ylab='Relative Abundance', main="Taxonomic composition of mice on D0 by cage", 
        col = colors, ylim=c(0,100), cex.lab=0.9, cex.axis=0.7, cex.names=0.7, xaxt='n', xlab="cage")
axis(1, at=1:26, labels=cage_only, cex.axis=0.6, cex.names=0.7)

#need to add a column of cages to the rel_d0 df

#what kind of plots do we want? line plots of diversity? 

