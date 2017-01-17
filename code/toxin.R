# C diff toxin code
#
#
# Dependencies:
# 
#
#

library('tidyr')
library('ggplot2')

meta_file   <- 'data/process/human_CdGF_metadata.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, row.names = 'sample_id')
toxin_file <- 'data/process/toxin_results.txt'
toxin_file <- read.table(toxin_file, sep = '\t', header = T)
sampleID_bins <- read.table('data/process/sampleID_bins.design', header=T)

toxin_file <- separate(toxin_file, Cage_Mouse, sep="-D", into = c("mouse_ID", "day"))
sampleID_bins <- separate(sampleID_bins, sample_ID, sep='-D', into = c("mouse_ID", "day"))
mouse_toxin <- merge(toxin_file, sampleID_bins, by.x="mouse_ID", by.y="mouse_ID")