###################
#
# HumanCdGF_microbiota_inocula.R
#
#Create file with microbiota community summary data
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
rel_abund <- 100*inocula_df/unique(apply(inocula_df, 1, sum))

#Create vector of OTUs with max abundances >1%
OTUs_1 <- apply(rel_abund, 2, max) > 1
OTU_list <- colnames(rel_abund)[OTUs_1]

#df of OTUs with abundances >1%
rel_abund_0.01 <- rel_abund[,OTUs_1]


#Load taxonomy function from tax_level.R
# Create dataframe with tax level as 
#  get_tax(
#    tax_level = level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom),
#    subsets rows by row_list (default = all OTUs),
#    dataframe dataframe (default = tax_file),
#    )
#
source('code/Sum_OTU_by_Tax.R')

phylum_level_RA <- sum_OTU_by_tax_level(5,rel_abund,tax_file)

# stacked barplot of tax comp (phylum) in inoculum samples
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

# merge all samples to compare differences between miseq v4 and 454 v53
v4_16s <- cbind(t(phylum_level_RA),names(phylum_level_RA))
colnames(v4_16s)[24] <- 'tax'
v53_454 <- cbind(t(phylum_454_level_RA),names(phylum_454_level_RA))
colnames(v53_454)[18] <- 'tax'

all_inocula <- merge(v53_454,v4_16s,by='tax')
all_inocula <- t(all_inocula)
colnames(all_inocula) <- all_inocula[1,]
all_inocula <- all_inocula[-1,]
rownames(all_inocula)[18] <- 'DA00581_1_inoc'
rownames(all_inocula)[40] <- 'DA00581_2_inoc'
sort_all_inocula <- all_inocula[order(rownames(all_inocula)),]

par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(sort_all_inocula), ylab='Relative Abundance', las=2, cex.sub=0.9,
        main="Comparison of Taxonomic Composition of GF mice inocula",cex.lab=0.9, cex.axis = 0.7, cex.names = 0.7,
        border= brewer.pal(8, 'Set2'), col = brewer.pal(8, 'Set2'), ylim=c(0,100),
        space=c(0.5,0.1,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.1,0.5,0.1,0.1,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5
                ,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5,0.1))
mtext('DA##### = v53 454, all others v4 MiSeq\nAll v4 samples from sequencing of inocula on 6/5/15 (AMSno, DA00581_1_inoc and DA00581_2_inoc from original seq set) ', side=1, line=7, cex=0.8)
legend('right',fill=brewer.pal(8, 'Set2'),legend=colnames(sort_all_inocula), inset=c(-0.16,0), bty='n', 
       border=brewer.pal(8, 'Set2'), y.intersp = 0.8, cex=0.9)
