###################
#
# thetayc_time.R
#
#Plot thetaYC over time distance from day 0 
#thetaYC distances cant have sample names that are just numbers, this script wont work 
#    
#
###################

source('code/read.dist.R')

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header=T)
theta <- read.dist(file='data/process/shared_group.thetayc.0.03.lt.dist', input="lt", make.square=T, diag=0)

#then maybe write it as a function 

#################THE RIGHT ONE!!!#######
#compares day0 to every subsequent day. this asks, how different are the communities from the start? can this be linked to outcome?

mouseIDVector <- as.character(unique(meta_file$mouse_id))
mouseIDVector <- mouseIDVector[mouseIDVector != "inoculum_NA"]
stability <- as.data.frame(matrix(ncol=10, nrow=length(mouseIDVector)))
row.names(stability) <- mouseIDVector
for (k in 1:10){
  sampA <- rownames(meta_file)[meta_file$day %in% 0]
  sampB <- rownames(meta_file)[meta_file$day %in% k]
  dist <- theta[sampA, sampB]
  row.names(dist) <- meta_file$mouse_id[rownames(meta_file) %in% rownames(dist)]
  names(dist) <- meta_file$mouse_id[rownames(meta_file) %in% names(dist)]
  for (i in mouseIDVector){
    coord <- dist[i,i]
    try(
    stability[i,k] <- as.numeric(coord), silent = TRUE
    )
  }
}

#plot





#build as function 


#now write version where you compare each day to previous day. this asks question: how stable is community over time?



