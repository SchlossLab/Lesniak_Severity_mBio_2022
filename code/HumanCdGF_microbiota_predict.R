###################
#
# HumanCdGF_microbiota_summary.R
#
# Predict disease outcome from community structure at Day 0
#    Disease outcome - 
#
#
###################
setwd("~/Documents/Github/Schubert_humanCdGF_XXXX_2016")

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header = T, row.names=1)
tax_file <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)

# subset files to day 0 data
day_0_meta <- meta_file[meta_file$day==0,]
day_0_meta <- day_0_meta[complete.cases(day_0_meta$day),-c(16:29)]
day_0_shared <- shared_file[rownames(shared_file) %in% rownames(day_0_meta),]

#Create df with relative abundances
rel_abund_day_0 <- 100*day_0_shared/unique(apply(day_0_shared, 1, sum))

#Create vector of OTUs with median abundances >1%
med_rel_abund_cage <- aggregate(rel_abund_day_0, by=list(day_0_meta$cage_id),median)
cage_IDs <- med_rel_abund_cage[,"Group.1"]
med_rel_abund_cage <- med_rel_abund_cage[,!colnames(med_rel_abund_cage) %in% 'Group.1']
OTUs_1 <- apply(med_rel_abund_cage, 2, max) > 1
OTU_list <- colnames(rel_abund_day_0)[OTUs_1]
#df of OTUs with abundances >1% - by cage and inoculum
rel_abund_day_0 <- rel_abund_day_0[,OTUs_1]

# mouse age of death
TOD <- meta_file[!meta_file$day %in% c(-7, 0,NA),names(meta_file) %in% c('day','mouse_id','cage_id')]
TOD$unique_id <- paste(TOD$mouse_id,TOD$cage_id)
Time_of_Death <- as.data.frame(apply(data.frame(do.call(rbind,split(TOD$day,c(TOD$unique_id)))),1,max))
names(Time_of_Death)[1] <- 'Endpoint'



#Load taxonomy function from tax_level.R
# Create dataframe with tax level as 
#  get_tax(
#    tax_level = level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom),
#    subsets rows by row_list (default = all OTUs),
#    dataframe dataframe (default = tax_file),
#    )
#
source('code/Sum_OTU_by_Tax.R')

family_level_RA <- sum_OTU_by_tax_level(2,rel_abund_day_0,tax_file)

colors <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# stacked barplot of tax comp (phylum) in inoculum samples
legend_labels <- as.character(names(family_level_RA))
#legend_labels[names(phylum_level_RA)=="Bacteria_unclassified"] <- "Other"
     
#stacked bar plot with base
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(family_level_RA), ylab='Relative Abundance', las=2, cex.sub=0.9,
        main="Taxonomic Composition of GF mice inocula",cex.lab=0.9, cex.axis = 0.7, cex.names = 0.7,
        border= colors, col = colors, ylim=c(0,100))
     mtext('Human Source', side=1, line=5, cex=0.9)
     legend('right',fill=colors,legend_labels, inset=c(-0.2,0), bty='n', border=colors, y.intersp = 0.8, cex=0.9)

#Create median rel abundance and IQRs
median_rel_abund <- aggregate()

low_QTR

up_QTR