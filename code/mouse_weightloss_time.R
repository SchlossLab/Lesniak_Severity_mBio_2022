##########
#
#
# This script builds the mouse weight loss figure
#
#
# Dependencies: 
#     * data/process/human_CdGF_metadata.txt
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#
#
# Output: * results/figures/weight_loss_time.pdf

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