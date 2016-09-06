###################
#
# HumanCdGF_microbiota_inocula.R
#
#Create file with microbiota community summary data + plot by family
#
#    analyze initial shared file from inocula - processed by mothur alone
#
###################

install.packages("RColorBrewer")
library(RColorBrewer)

setwd("~/Documents/Github/Schubert_humanCdGF_XXXX_2016")

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared_file <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header = T, row.names=1)
tax_file <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)

#subset shared to be inoculum only
inocula_df <- shared_file[meta_file$cage_id=='inoculum',]

#Create df with relative abundances
rel_abund_inoc <- 100*inocula_df/unique(apply(inocula_df, 1, sum))

#Create vector of OTUs with max abundances >1%
OTUs_1_inoc <- apply(rel_abund_inoc, 2, max) > 1
OTU_list_inoc <- colnames(rel_abund_inoc)[OTUs_1_inoc]

#df of OTUs with abundances >1%
rel_abund_0.01 <- rel_abund_inoc[,OTUs_1_inoc]


#Load taxonomy function from tax_level.R
# Create dataframe with tax level as 
#  get_tax(
#    tax_level = level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom),
#    subsets rows by row_list (default = all OTUs),
#    dataframe dataframe (default = tax_file),
#    )
#
source('code/Sum_OTU_by_Tax.R')

family_level_RA <- sum_OTU_by_tax_level(2,rel_abund_0.01,tax_file)

#reorganize the file a bit 
inoc <- meta_file$human_source[meta_file$cage_id=='inoculum']
inoc <- as.vector(inoc)
family_level_RA[25] <- inoc
colnames(family_level_RA)[25] <- "donor"

#remove donor 810 and duplicates of others

family_level_RA <- family_level_RA[-8,]
family_level_RA <- family_level_RA[-1,]
family_level_RA <- family_level_RA[-4,]
family_level_RA <- family_level_RA[-17,]

inoc_labels <- colnames(family_level_RA[1:24])
inoc_list <- as.vector(inoc_labels)

#make multiplot of family plots for each donor
#run this entire code to the end to plot the multipanel plot for all donors
#set up layout
par(mar=c(2,1.5,1.1,1), oma=c(6,7,6,0.5))
layout(matrix(1:16, nrow=4))

#loop to plot all donors 
for(d in family_level_RA$donor){
  one_d <- subset(family_level_RA, donor == d)
  bug_levels <- t(one_d[,1:24])
  #make barplot
  z <- barplot(t(one_d[,1:24]), beside = TRUE, xaxt='n', col="white", ylim=c(0,max(one_d[,1:24])+5), ylab = "Relative Abundance (%)", main=paste("Donor", d))
  axis(1, at=seq(1,24, by=1), labels=FALSE, tick=FALSE)
  #text(seq(1,24, by=1), par("usr")[3] - 11, labels=inoc_list, srt = 65, pos = 1, xpd = TRUE, cex=0.6)
  box()
}

#tweak final plot so it looks nice 
#text(x=z+0.8, y=-3, xpd=NA, label=parse(text=inoc_list), pos=2, srt=70, cex=0.8)
#text(x=z-28.5, y=-100, xpd=NA, label=parse(text=inoc_list), pos=2, srt=70, cex=0.8)
#text(x=z-57.5, y=-100, xpd=NA, label=parse(text=inoc_list), pos=2, srt=70, cex=0.8)
text(x=z+0.8, y=-5, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-29, y=-5, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-58, y=-5, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-88, y=-5, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
mtext("Relative Abundance (%)", side = 2, line =2, las = 3, cex = 1, adj=-1, padj=-73)

#mtext("Relative Abundance (%)", side = 2, line =2, las = 3, cex = 1, adj=0, padj=-55)
#mtext("Donor Inocula Communities", side = 3, cex = 1, las = 1, adj =6.5, padj=-48)

dev.off()

#barplot of donors for comparison of figure types 
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p = 17
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
x <- barplot(t(family_level_RA[,1:24]), ylab= "Relative Abundance %", ylim=c(0,100), xaxt='n', col = getPalette(p), names.arg=family_level_RA$donor, cex.names = 0.9, cex.axis=0.9, cex.lab=0.9)
legend_labels <- colnames(family_level_RA[,1:24])
legend_labels <- as.character(legend_labels)
legend("right",fill=getPalette(p),legend_labels, bty='n', cex=0.6, y.intersp = 0.7, inset=c(-0.25,0))
box()
text(cex=0.8, x=x+0.5, y=-4, label=family_level_RA$donor, xpd=TRUE, srt=45, pos=2)

