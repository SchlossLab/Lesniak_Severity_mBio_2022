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
shared <- read.table('data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared', sep='\t',header=T)

#subset shared to be inoculum only
inocula_df <- shared_file[meta_file$cage_id=='inoculum',]
inocula_shared <- shared[shared$Group %in% rownames(inocula_df),]

#take out samples we cant use
inocula_shared <- inocula_shared[-1,]
inocula_shared <- inocula_shared[-4,]
inocula_shared <- inocula_shared[-6,]
inocula_shared <- inocula_shared[-17,]

#write file out so we can use it in mothur
write.table(inocula_shared, file='data/process/inocula.shared', quote=F,sep='\t',row.names=F)

#run dist.shared, nmds in mothur to get axes file

nmds <- read.table(file='data/process/inocula.thetayc.0.03.lt.ave.nmds.axes', header = T)

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








