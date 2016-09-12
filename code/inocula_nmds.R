###################
#
# Inocula_nmds.R
#
#Create shared file of inocula only and plot NMDS and color points by outcome
#
#    
#
###################

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header=T)

#orig_shared <- read.table('data/mothur/gf_all.an.0.03.subsample.shared', sep='\t', header = T)
#remove dashes from group names so the next stuff will work
#orig_shared[2] <- gsub("-", "", orig_shared$Group)

#subset shared to be inoculum only
inocula <- rownames(meta_file)[meta_file$cage_id=='inoculum']
inocula_shared <- shared[shared$sample_id %in% inocula,]

#subset shared to be day 0 only- will use later for comparison 
day0 <- rownames(meta_file)[meta_file$day==0]
day0_shared <- shared[shared$sample_id %in% day0,]

#add back label and numOTUs column to subsetted shared file-- need this to run dist.shared
inocula_shared <- cbind(label=0.03, inocula_shared[1], numOTUs=2255, inocula_shared[2:ncol(inocula_shared)])
colnames(inocula_shared)[2] <- 'Group'

#for day 0 also
day0_shared <- cbind(label=0.03, day0_shared[1], numOTUs=2255, day0_shared[2:ncol(day0_shared)])
colnames(day0_shared)[2] <- 'Group'

#combine these all to run as one file duhhhh
inputd0_shared <- rbind(inocula_shared, day0_shared)

#write file out so we can use it in mothur
write.table(inocula_shared, file='data/process/inocula.shared', quote=F,sep='\t',row.names=F)
write.table(day0_shared, file='data/process/day0.shared', quote=F,sep='\t',row.names=F)
write.table(inputd0_shared, file='data/process/inputd0.shared', quote=F, sep='\t', row.names=F)

#run dist.shared, nmds in mothur to get axes file
#mothur ouput stats for combined file: 
#Number of dimensions:	2
#Lowest stress :	0.362683
#R-squared for configuration:	0.380491

combo_nmds <- read.table(file='data/process/inputd0.thetayc.0.03.lt.nmds.axes', header = T)
inoc_nmds <- read.table(file='data/process/inocula.thetayc.0.03.lt.nmds.axes', header = T)
day0_nmds <- read.table(file='data/process/day0.thetayc.0.03.lt.nmds.axes', header = T)




#need to color points by type 
inocula_design <- read.table(file='data/raw/donor.design', header = F)
inocula_nmds <- merge(nmds, inocula_design, by.x='group', by.y='V1')

severe <- inocula_nmds[grep('Severe', inocula_nmds$V2), c(2,3)]
asymptomatic <- inocula_nmds[grep('Asymptomatic', inocula_nmds$V2), c(2,3)]

#make nmds inocula plot
plot(inocula_nmds$axis1, inocula_nmds$axis2, main="Similarity of donor inocula")
points(severe, pch=16, col = "red")
points(asymptomatic, pch=16, col = "black")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)








