#############
#
# mouse_d0_community.R
#
# Create microbiota files for day 0 and plot different ways
#
#
#
#############

install.packages("RColorBrewer")
library(RColorBrewer)

#import Nick's awesome combined metadata file
meta_file <- read.table("data/process/human_CdGF_metadata.txt", header = TRUE, sep='\t', fill = TRUE, row.names=3)
shared_file <- read.table("data/process/human_CdGF.an.unique_list.0.03.subsample.shared", sep='\t', header = TRUE, row.names=1)
tax_file <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)

#make OTU abundance file

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


# Create dataframe with tax level as 
#  get_tax(
#    tax_level = level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom),
#    subsets rows by row_list (default = all OTUs),
#    dataframe dataframe (default = tax_file),
#    )
#
source('code/Sum_OTU_by_Tax.R')

#ok but this stuff is giving me more than just d0
taxonomy_genus <- sum_OTU_by_tax_level(1, rel_d0, tax_file)
taxonomy_family <- sum_OTU_by_tax_level(2, rel_d0, tax_file)
taxonomy_phylum <- sum_OTU_by_tax_level(5, rel_d0, tax_file)

#plot each mouse individually

colors <- c('red', 'blue', 'green', 'orange', 'purple', 'pink')
legend_labels <- as.character(taxonomy_phylum$tax)
legend_labels[5] <- "Other"

barplot(t(taxonomy_phylum), ylab='Relative Abundance', main="Taxonomic composition of mice on D0", 
        col = colors, ylim=c(0,100), cex.lab=0.9, cex.axis=0.7, cex.names=0.7)

#combine by cages/donors/ median per cage 

#add a column for cages to the df
cages <- meta_file$cage_id[meta_file$day==0]
cages <- na.omit(cages)
rel_d0[82] <- cages
colnames(rel_d0)[82] <- "cage"

#add a column of donors to the df
donors <- meta_file$human_source[meta_file$day==0]
donors <- na.omit(donors)
rel_d0[83] <- donors
colnames(rel_d0)[83] <- "donor"
 
d0_med <- aggregate(rel_d0[, 1:81], list(rel_d0$cage), median)
d0_med_phy <- sum_OTU_by_tax_level(5, d0_med, tax_file)
#removes inoculum cage
cage_only <- cage_IDs[-15]
cage_only_char <- as.character(cage_only)

#use these to change dataframe rownames, cleans up for plotting later
rownames(d0_med_phy) <- cage_only_char
d0_med_phy <- d0_med_phy[,-1]

#run this whole command to make the figure of phylum by cage/donor
barplot(t(d0_med_phy), ylab='Relative Abundance', main="Taxonomic composition of mice on D0 by cage", 
        col = colors, ylim=c(0,100), cex.lab=0.9, cex.axis=0.7, cex.names=0.7, xlab="cage")
legend_labels <- colnames(d0_med_phy)
legend('right',fill=colors,legend_labels, bty='n', inset=c(-0.1,0), border=colors, y.intersp = 0.8, cex=0.9)

#make dataframe of median family abundances and organize it 
d0_med_fam <- sum_OTU_by_tax_level(2, d0_med, tax_file)
rownames(d0_med_fam) <- cage_only_char
colnames(d0_med_fam)[1] <- "cage"

#just make table of human source and cage unique and then merge by cage, take out unused ones after 

dc <- meta_file[!meta_file$cage_id=='inoculum',2:3]
uc <- unique(dc)
d0_med_fam <- merge(d0_med_fam, uc, by.x='cage', by.y='cage_id')

#get rid of samples/cages not used here: 
d0_med_fam <- d0_med_fam[-23,]
d0_med_fam <- d0_med_fam[-4,]
d0_med_fam <- d0_med_fam[-3,]
d0_med_fam <- d0_med_fam[-6,]
d0_med_fam <- d0_med_fam[-4,]
d0_med_fam <- d0_med_fam[-4,]
d0_med_fam <- d0_med_fam[-12,]



#colorbrewer
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#run this whole command to make figure of family by cage/donor - add colorbrewer!
#main="Taxonomic composition of mice on D0 by donor, family level"
n=16
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(d0_med_fam[2:20]), ylab='Relative Abundance', 
        col = getPalette(n), ylim=c(0,100), cex.lab=0.9, cex.axis=0.8, cex.names=0.8, xaxt = 'n')
fam_labels <- colnames(d0_med_fam[2:20])
legend('right', fill = getPalette(n), fam_labels, inset=c(-0.15,0), cex=0.6)
text(cex=0.8, x=x+0.5, y=-4, label=d0_med_fam$human_source, xpd=TRUE, srt=45, pos=2)
box()


#barplots for individual families
#combine into loop and output plots to a main figure looking file 

#calculates upper and lower quantiles, then combines by family. 
d0_upper <- aggregate(rel_d0[1:81], list(rel_d0$cage), FUN= quantile, probs =0.75)
d0_lower <- aggregate(rel_d0[1:81], list(rel_d0$cage), FUN= quantile, probs =0.25)
d0_upper_fam <- sum_OTU_by_tax_level(2, d0_upper, tax_file)
d0_lower_fam <- sum_OTU_by_tax_level(2, d0_lower, tax_file)
colnames(d0_upper_fam)[1] <- "donor"
colnames(d0_lower_fam)[1] <- "donor"
fam_list <- as.vector(fam_labels)

#to save as pdf to results, can't get the size right to render but fine for now
pdf(file="results/figures/figure2.pdf", width=25, height=15)

#run this entire code to the end to plot the multipanel plot for all donors
#set up layout
par(mar=c(2,1.5,1.1,1), oma=c(6,7,6,0.5))
layout(matrix(1:30, nrow=5))

#loop to plot all donors 
for(d in d0_med_fam$cage){
  one_d <- subset(d0_med_fam, cage == d)
  one_uci <- subset(d0_upper_fam, cage == d)
  one_lci <- subset(d0_lower_fam, cage == d) 
  medians <- t(one_d[,2:20])
  uci <- t(one_uci[,2:20])
  lci <- t(one_lci[,2:20])
  #make barplot
  z <- barplot(t(one_d[,2:20]), beside = TRUE, xaxt='n', col="white", ylim=c(0,max(uci)+5), ylab = "Relative Abundance (%)", main=paste("Cage", d))
  #axis(1, at=seq(1,19, by=1), labels=FALSE, tick=FALSE)
  #text(seq(1,19, by=1), par("usr")[3] - 11, labels=fam_list, srt = 65, pos = 1, xpd = TRUE, cex=0.6)
  box()
  #plot error bars
  arrows(x0=z, y0=medians, z, y1=uci, angle=90, length=0.05)
  arrows(x0=z, y0=medians, z, y1=lci, angle=90, length=0.05)
}

#tweak final plot so it looks nice 
text(x=z+0.8, y=-5, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-95, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-71, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-47, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-23, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
mtext("Relative Abundance (%)", side = 2, line =2, las = 3, cex = 1, adj=1, padj=-78)
mtext("Mouse community on day 0, by cage", side = 3, cex = 1, adj =7, padj=-22)

dev.off()