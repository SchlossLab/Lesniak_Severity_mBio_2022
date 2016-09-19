###################
#
# thetayc_time.R
#
#Plot thetaYC over time distance from day 0 
#
#    
#
###################

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header=T)