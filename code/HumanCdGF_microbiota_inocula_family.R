###################
#
# HumanCdGF_microbiota_inocula.R
#
#Create file with microbiota community summary data + plot by family
#
#    analyze initial shared file from inocula - processed by mothur alone
#
###################

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
family_level_RA[25] <- inoc
colnames(family_level_RA)[25] <- "donor"

#remove donor 810 because we can't use it anymore

family_level_RA <- family_level_RA[-8,]

inoc_labels <- colnames(family_level_RA[1:24])
inoc_list <- as.vector(inoc_labels)

#make multiplot of family plots for each donor
#run this entire code to the end to plot the multipanel plot for all donors
#set up layout
par(mar=c(2,1.5,1.1,1), oma=c(6,7,6,0.5))
layout(matrix(1:30, nrow=5))

#loop to plot all donors 
for(d in family_level_RA$donor){
  one_d <- subset(family_level_RA, donor == d)
  bug_levels <- t(one_d[,1:24])
  #make barplot
  z <- barplot(t(one_d[,1:24]), beside = TRUE, xaxt='n', col="white", ylim=c(0,max(one_d[,1:24])+5), ylab = "Relative Abundance (%)", main=paste("Donor", d))
  axis(1, at=seq(1,24, by=1), labels=FALSE, tick=FALSE)
  text(seq(1,24, by=1), par("usr")[3] - 11, labels=inoc_list, srt = 65, pos = 1, xpd = TRUE, cex=0.6)
  box()
}

#tweak final plot so it looks nice 
text(x=z+0.8, y=-5, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-95, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-71, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-47, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
text(x=z-23, y=-225, xpd=NA, label=parse(text=fam_list), pos=2, srt=70, cex=0.8)
mtext("Relative Abundance (%)", side = 2, line =2, las = 3, cex = 1, adj=1, padj=-78)
mtext("Mouse community on day 0, by cage", side = 3, cex = 1, adj =7, padj=-22)







# stacked barplot of tax comp (family) in inoculum samples
colors <- primary.colors(length(names(phylum_level_RA)))
legend_labels <- names(phylum_level_RA)
legend_labels[names(phylum_level_RA)=="Bacteria_unclassified"] <- "Other"
legend_labels[names(phylum_level_RA)=="Candidatus_Saccharibacteria"] <- "C. Saccharibacteria"
legend_labels[names(phylum_level_RA)=="Deinococcus-Thermus"] <- "D. Thermus"

#stacked bar plot with base
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(phylum_level_RA), ylab='Relative Abundance', las=2, cex.sub=0.9,
        main="Taxonomic Composition of GF mice inocula",cex.lab=0.9, cex.axis = 0.7, cex.names = 0.7,
        border= colors, col = colors, ylim=c(0,100))
     mtext('Human Source', side=1, line=5, cex=0.9)
     legend('right',fill=colors,legend_labels, inset=c(-0.16,0), bty='n', border=colors, y.intersp = 0.8, cex=0.9)

# community structure from 454 data
shared_454_file <- read.table('../untracked_Schubert_humanCdGF_XXXX_2016/DA_v53_454/final.an.0.03.subsample.shared', sep='\t',header = T, row.names=2)
tax_454_file <- read.table('../untracked_Schubert_humanCdGF_XXXX_2016/DA_v53_454/final.an.0.03.cons.taxonomy', sep='\t',header = T, row.names = 1)
shared_454_file <- shared_454_file[,-c(1,2)]
rel_454_abund <- 100*shared_454_file/unique(apply(shared_454_file, 1, sum))

phylum_454_level_RA <- sum_OTU_by_tax_level(5,rel_454_abund,tax_454_file)

legend_454 <- names(phylum_454_level_RA)
legend_454[names(phylum_454_level_RA)=="Bacteria_unclassified"] <- "Other"

#stacked bar plot with base of 454
par(mar=c(7.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(phylum_454_level_RA), ylab='Relative Abundance', las=2, cex.sub=0.9,
        main="Taxonomic Composition of GF mice inocula",cex.lab=0.9, cex.axis = 0.7, cex.names = 0.7,
        border= colors, col = colors, ylim=c(0,100))
mtext('Human Source', side=1, line=5, cex=0.9)
legend('right',fill=colors,legend_454, inset=c(-0.16,0), bty='n', border=colors, y.intersp = 0.8, cex=0.9)


# merge all samples to compare differences between miseq v4 and 454 v53 - at the family level
family_level_RA <- sum_OTU_by_tax_level(2,rel_abund,tax_file)
family_454_level_RA <- sum_OTU_by_tax_level(2,rel_454_abund,tax_454_file)

v4_16s <- cbind(as.data.frame(t(family_level_RA)),tax=names(family_level_RA))
v53_454 <- cbind(as.data.frame(t(family_454_level_RA)),tax=names(family_454_level_RA))

all_inocula <- merge(v53_454,v4_16s,by='tax')
rownames(all_inocula) <- all_inocula[,'tax']
all_inocula <- all_inocula[,-1]
all_inocula_OTUs_1 <- apply(all_inocula, 1, mean) > 1
less <- all_inocula[!all_inocula_OTUs_1,]
less_abund <- apply(less, 2, sum)
family_inocula <- rbind(all_inocula[all_inocula_OTUs_1,],'Less Abundant Taxa'=less_abund)

names(family_inocula)[names(family_inocula)=='581-inoculum'] <- 'DA00581_1_inoc'
names(family_inocula)[names(family_inocula)=='DA581'] <- 'DA00581_2_inoc'
sort_all_inocula <- family_inocula[,order(names(family_inocula))]

par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(as.matrix(sort_all_inocula), ylab='Relative Abundance', las=2, cex.sub=0.9,
        main="Comparison of Taxonomic Composition of GF mice inocula",cex.lab=0.9, cex.axis = 0.7, cex.names = 0.7,
        border= brewer.pal(length(rownames(sort_all_inocula)), 'Paired'), col = brewer.pal(length(rownames(sort_all_inocula)), 'Paired'), ylim=c(0,100),
        space=c(0.5,0.1,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5
                ,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1))
mtext('DA##### = v53 454, all others v4 MiSeq\nAll v4 samples from sequencing of inocula on 6/5/15 (AMSno, DA00581_1_inoc and DA00581_2_inoc from original seq set) ', side=1, line=7, cex=0.8)
legend('right',fill=brewer.pal(length(rownames(sort_all_inocula)), 'Paired'),legend=rownames(sort_all_inocula), inset=c(-0.16,0), bty='n', 
       border=brewer.pal(length(rownames(sort_all_inocula)), 'Paired'), y.intersp = 0.8, cex=0.9)

# - which genera are increased in the june 2015 inoculum samples?
family_OTU_list <-get_tax(2,rownames(tax_file),tax_file)
genus_OTU_list <- get_tax(1,rownames(tax_file),tax_file)
entero <- genus_OTU_list[family_OTU_list$tax=='Enterobacteriaceae',]
pseudo <- genus_OTU_list[family_OTU_list$tax=='Pseudomonadaceae',]
pseudo_list <- unique(pseudo)
entero_list <- unique(entero)
genus_level_RA <- sum_OTU_by_tax_level(1,rel_abund,tax_file)
genus_level_RA[,names(genus_level_RA) %in% pseudo_list]
genus_level_RA[,names(genus_level_RA) %in% entero_list]
