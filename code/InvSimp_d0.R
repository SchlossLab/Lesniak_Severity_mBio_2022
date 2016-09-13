###################
#
# InvSimp_d0.R
#
#Plot Inverse Simpson diversity metric for day 0 communities +/or over time
#
#    
#
###################

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header=T)
simp <- read.table('data/mothur/gf_all.an.0.03.subsample.groups.summary', sep='\t',header=T)

#subset simp to be day 0 only
day0list <- rownames(meta_file)[meta_file$day==0]
simp[2] <- gsub("-", "", simp$group)
day0list <- gsub("-", "", day0list)

#for the below to work, need to do some sort of gsubbing to get the dashes fixed 
day0_simp <- simp[simp$group %in% day0list,]
day0_meta <- subset(meta_file, day==0)

day0_cages <- day0_meta$cage_id
day0_simp$cage <- day0_meta$cage_id

fullsimp <- cbind(day0_simp[2:3], day0_simp[7])

plot(fullsimp$cage, fullsimp$invsimpson)


#plot day0 inv simpson with cage along bottom inv simpson on top
#need to merge w metadata file to get column of cages 

#but would be better to show over time like i did with the cdiff time paper plot 0-10 days, by donor? cage? but dont want to do average, too small n 
