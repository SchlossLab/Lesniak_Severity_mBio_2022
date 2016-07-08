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
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 'group')
shared_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header = T, row.names='Group')
tax_file <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 'OTU')

     
#### subset files to day 0 data
day_0_meta <- meta_file[which(meta_file$day==0),!is.na(meta_file['430-307-D0',])]
# subset meta to day 1 data
cdiff_effect <- c("cdiff_cfu","nextDayCFU_colon","nextDayCFU_any","weight_grams","percent_weightLoss.x","Number","Experiment_day",
                  "loss_below90","day_end","Recip_Dilution_from_DP2","Recip_Dilution_from_adding_to_cells",
                  "Recip_Dilution_from_Postion_on_PlateA","Recip_Dilution_from_Postion_on_PlateB","Average.",
                  "Log_repiricoal_dilution")
day_1_cdiff <- meta_file[which(meta_file$day==1),c(cdiff_effect, 'Mouse_ID')]
names(day_1_cdiff)[!names(day_1_cdiff) %in% 'Mouse_ID'] <- 
     paste('Day_1_',names(day_1_cdiff)[!names(day_1_cdiff) %in% 'Mouse_ID'], sep='')
#subset meta_file to endpoint data
stool_meta_file <- meta_file[!meta_file$sample_type=='cecum',]
day_euthanized_cdiff <- stool_meta_file[which(stool_meta_file$day == stool_meta_file$Euthanized), c(cdiff_effect,'Mouse_ID')]
names(day_euthanized_cdiff)[!names(day_euthanized_cdiff) %in% 'Mouse_ID'] <- 
     paste('Day_Euth_',names(day_euthanized_cdiff)[!names(day_euthanized_cdiff) %in% 'Mouse_ID'], sep='')

#### subset shared
#Create df with relative abundances
day_0_shared <- shared_file[rownames(day_0_meta),]
rel_abund_day_0 <- 100*day_0_shared/unique(apply(day_0_shared, 1, sum))
#Create vector of OTUs with median abundances >1%
day_0_shared_median <- aggregate(rel_abund_day_0, by=list(day_0_meta$cage_id),median)
#cage_IDs <- day_0_shared_median[,"Group.1"]
day_0_shared_median <- day_0_shared_median[,!colnames(day_0_shared_median) %in% 'Group.1']
OTUs_day_0_01 <- apply(day_0_shared_median, 2, max) > 1
OTU_list <- colnames(rel_abund_day_0)[OTUs_1]
#df of OTUs with abundances >1% - by cage and inoculum
rel_abund_day_0 <- rel_abund_day_0[,OTUs_day_0_01]

#Merge data by mouse
Disease_predict_df <- merge(
          merge(day_1_cdiff[,c('Mouse_ID','Day_1_cdiff_cfu','Day_1_percent_weightLoss.x','Day_1_Average.','Day_1_Log_repiricoal_dilution')],
                day_euthanized_cdiff[,c('Mouse_ID','Day_Euth_cdiff_cfu','Day_Euth_percent_weightLoss.x','Day_Euth_Average.','Day_Euth_Log_repiricoal_dilution')], by='Mouse_ID',all.x=T),
          merge(day_0_meta[,c('Mouse_ID','cage_id','Euthanized')],rel_abund_day_0, by='row.names'),
          by='Mouse_ID',all.y=T)
Predict_df <- Disease_predict_df[,!names(Disease_predict_df) %in% c('Mouse_ID','Row.names','cage_id')]
Predict_median_cage <- aggregate(Predict_df, by=list(Disease_predict_df$cage_id), median, na.rm=TRUE)
Predict_median_cage <- data.frame(Predict_median_cage[,!names(Predict_median_cage) %in% 'Group.1'], row.names = Predict_median_cage$Group.1)

## plot cfu
library(ggplot2)
CFU_data <- stool_meta_file[-which(stool_meta_file$day==-7),c('cage_id','day','human_source','cdiff_cfu','Log_repiricoal_dilution','percent_weightLoss.x','Euthanized')]
CFU_431_data <- CFU_data[!CFU_data$cage_id %in% c('CON','CON2','CON3'),]
CFU_data_median <- aggregate(CFU_431_data[,c(4:7)], by=list(CFU_431_data$cage_id,CFU_431_data$day), median, na.rm=TRUE)
names(CFU_data_median)[names(CFU_data_median) %in% c('Group.1','Group.2')] <- c('cage_id','day')
Cage_Inoc <- unique(meta_file[!meta_file$cage_id=='inoculum',c('cage_id','human_source')])
CFU_data_median <- merge(CFU_data_median,Cage_Inoc, by='cage_id')

cfu <- cbind(CFU_data_median[,c(1,2,7,3)],data_type=rep('cfu',length(CFU_data_median$human_source)))
toxin<- cbind(CFU_data_median[,c(1,2,7,4)],data_type=rep('toxin',length(CFU_data_median$human_source)))
wt_loss <- cbind(CFU_data_median[,c(1,2,7,5)],data_type=rep('wt_loss',length(CFU_data_median$human_source)))
names(cfu)[names(cfu) %in% 'cdiff_cfu'] <- 'data'
names(toxin)[names(toxin) %in% 'Log_repiricoal_dilution'] <- 'data'
names(wt_loss)[names(wt_loss) %in% 'percent_weightLoss.x'] <- 'data'

Outcome_cat_df <- merge(unique(CFU_data_median[,c('cage_id','Euthanized')]),CFU_data_median[CFU_data_median$day==1,c('cdiff_cfu','cage_id')],by='cage_id')
Outcome_cat_df$Euthanized[Outcome_cat_df$Euthanized %in% 10] <- 'Day 10'
Outcome_cat_df$Euthanized[!Outcome_cat_df$Euthanized %in% 'Day 10'] <- 'Early'
Outcome_cat_df$cdiff_cfu[Outcome_cat_df$cdiff_cfu <= 10000000] <- 'Low CFU'
Outcome_cat_df$cdiff_cfu[!Outcome_cat_df$cdiff_cfu == 'Low CFU'] <- 'High CFU'
Outcome_cat_df$cdiff_cfu[Outcome_cat_df$Euthanized == 'Early'] <- 'Euthanized'

facet_inoc_df <- rbind(rbind(cfu,toxin),wt_loss)
facet_inoc_df <- merge(facet_inoc_df,Outcome_cat_df,by='cage_id')
facet_df <- cbind(cbind(cage_id=paste(facet_inoc_df$cage_id,facet_inoc_df$human_source, sep = '_'),facet_inoc_df[,c('day','data','data_type','Euthanized')]),
                       outcome=paste(facet_inoc_df$Euthanized,facet_inoc_df$cdiff_cfu, sep = '_'))

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

# by cage - faceted by datatype
# not by human_source because human source has different treatments and day 0 communities for same treatments (DA00581)
ggplot(data=facet_df, aes( x= day, y=data, color=cage_id)) + geom_point() + geom_line() + scale_y_log10() + 
     facet_grid(data_type ~ Euthanized, scales='free_y') +
     theme_bw() + xlab('Day') + ylab('') + scale_x_continuous(breaks=c(0:10)) + scale_color_manual(values=colors)



library(randomForest)

important <- importance(randomForest(Euthanized ~ ., data=Predict_median_cage[,c(9:117)], na.action=na.omit))
important <- importance(randomForest(Day_Euth_cdiff_cfu ~ ., data=Predict_median_cage[,c(5,10:117)], na.action=na.omit))
important <- importance(randomForest(Day_Euth_Log_repiricoal_dilution ~ ., data=Predict_median_cage[,c(8,10:117)], na.action=na.omit))

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