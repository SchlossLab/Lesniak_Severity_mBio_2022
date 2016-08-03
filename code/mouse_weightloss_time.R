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
