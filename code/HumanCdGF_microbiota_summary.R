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
rel_abund_cage <- rel_abund[!meta_file$cage_id=='inoculum',OTUs_1]
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

taxonomy_genus <- get_tax(1)
taxonomy_phylum <- get_tax(5)

#missing inoculum

shared_presample <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared', sep='\t',header = T)

needed_inoculum <- as.character(unique(meta_file[meta_file$cage_id!='inoculum','human_source']))
sequenced_inoculum <- as.character(unique(meta_file[meta_file$cage_id=='inoculum','human_source']))
needed_inoculum <- needed_inoculum[!needed_inoculum %in% sequenced_inoculum]
inoculum_presample <- shared_presample[,'Group']
inoculum_presample <- inoculum_presample[grep('*inoc*',inoculum_presample)]

#create function to output a dataframe with OTU rel abund summed by tax group per sample
# input file with tax level classs and file with samples OTU rel abund
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

phylum_level_RA <- sum_OTU_by_tax_level(taxonomy_phylum, inoculum_rel_abund)


# stacked barplot of tax comp (phylum) in inoculum samples
colors <- c('#0066FF','#CC0033','#66CC33','#FF9933','#CC3399','#666666')
legend_labels <- as.character(unique(taxonomy_phylum$tax))
legend_labels[5] <- "Other"
     
#stacked bar plot with base
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(phylum_level_RA), ylab='Relative Abundance', las=2, cex.sub=0.9,
        main="Taxonomic Composition of GF mice inoculum",cex.lab=0.9, cex.axis = 0.7, cex.names = 0.7,
        border= colors, col = colors, ylim=c(0,100))
     mtext('Human Source', side=1, line=5, cex=0.9)
     legend('right',fill=colors,legend_labels, inset=c(-0.2,0), bty='n', border=colors, y.intersp = 0.8, cex=0.9)

#stacked bar plot with ggplot
install.packages('ggplot2')
library(ggplot2)
install.packages("reshape2")
library(reshape2)
ggplot_phylum_RA <- cbind(phylum_level_RA,Names= rownames(phylum_level_RA))
ggplot_phylum_RA <- melt(ggplot_phylum_RA, id.vars='Names')
ggplot(ggplot_phylum_RA, aes(x=Names, y=value, fill=variable)) +
     geom_bar(stat='identity')+
     xlab('Human Source')+
     ylab("Relative Abundance")+
     ggtitle("Taxonomic Composition of GF mice inoculum") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
     scale_fill_manual(values=colors, name='', labels=legend_labels)+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line.y = element_line(colour = "black"),
           axis.ticks.x = element_blank())

#Create median rel abundance and IQRs
median_rel_abund <- aggregate()

low_QTR

up_QTR