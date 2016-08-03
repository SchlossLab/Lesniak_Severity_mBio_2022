##########
#
#
# This script builds a plot of abundances of each phyla or family over time
# this script is incomplete as of 8/3/16 KF
#
#
# Dependencies: 
#     * data/process/human_CdGF_metadata.txt
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#
#
# Output: * results/figures/phyla_over_time.pdf

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

#then sum OTU by tax level to get species names
source('code/Sum_OTU_by_Tax.R')
cage_fam <- sum_OTU_by_tax_level(2, rel_abund_cage, tax_file)
cage_phy <- sum_OTU_by_tax_level(5, rel_abund_cage, tax_file)

#does anything else need to go here? 


#get median and IQR by cage and day
relabund_cageday <- na.omit(relabund_cageday)

# transform <- t(relabund_cageday), try a solution without transforming
lci <- aggregate(rel_abund, by = list(meta_file$cage_id, meta_file$day), function(x){quantile(x, probs=c(0.25), na.rm=TRUE)})
uci <- aggregate(rel_abund, by = list(meta_file$cage_id, meta_file$day), function(x){quantile(x, probs=c(0.75), na.rm=TRUE)})

