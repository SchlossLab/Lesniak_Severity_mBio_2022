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


#merge design and nmds files to prep for plotting 
inocula_design <- read.table(file='data/raw/donor.design', header = F)
inocula_nmds <- merge(inoc_nmds, inocula_design, by.x='group', by.y='V1')

combo_design <- read.table(file='data/raw/d0_outcome.design', header = F)
combo_nmds <- merge(combo_nmds, combo_design, by.x='group', by.y='V1')

d0_donor <- read.table(file='data/raw/sample_donor.design', header = F)
day0_nmds <- merge(day0_nmds, d0_donor, by.x='group', by.y='V1')


#make nmds inocula plot
severe <- inocula_nmds[grep('Severe', inocula_nmds$V2), c(2,3)]
asymptomatic <- inocula_nmds[grep('Asymptomatic', inocula_nmds$V2), c(2,3)]
plot(inocula_nmds$axis1, inocula_nmds$axis2, main="Similarity of donor inocula")
points(severe, pch=16, col = "red")
points(asymptomatic, pch=16, col = "black")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)
#amova in mothur output for asymptomatic vs severe: p-value: 0.508

#donor and day 0 by outcome plot
severe2 <- combo_nmds[grep('Severe', combo_nmds$V2), c(2,3)]
asymp2 <- combo_nmds[grep('Asymptomatic', combo_nmds$V2), c(2,3)]
plot(combo_nmds$axis1, combo_nmds$axis2, main="Similarity of donor and day0 by outcome")
points(severe2, pch=16, col = "red")
points(asymp2, pch=16, col = "black")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)

#day 0 by outcome plot, no donors, need to make separate nmds for this one... ugh. 



#nmds colored by donor plot
d369 <- day0_nmds[grep('DA00369', day0_nmds$V2), c(2,3)]
d430 <- day0_nmds[grep('DA00430', day0_nmds$V2), c(2,3)]
d431 <- day0_nmds[grep('DA00431', day0_nmds$V2), c(2,3)]
d578 <- day0_nmds[grep('DA00578', day0_nmds$V2), c(2,3)]
d581 <- day0_nmds[grep('DA00581', day0_nmds$V2), c(2,3)]
d884 <- day0_nmds[grep('DA00884', day0_nmds$V2), c(2,3)]
d953 <- day0_nmds[grep('DA00953', day0_nmds$V2), c(2,3)]
d1134 <- day0_nmds[grep('DA01134', day0_nmds$V2), c(2,3)]
d1146 <- day0_nmds[grep('DA01146', day0_nmds$V2), c(2,3)]
d1245 <- day0_nmds[grep('DA01245', day0_nmds$V2), c(2,3)]
d1324 <- day0_nmds[grep('DA01324', day0_nmds$V2), c(2,3)]
d10027 <- day0_nmds[grep('DA10027', day0_nmds$V2), c(2,3)]
d10034 <- day0_nmds[grep('DA10034', day0_nmds$V2), c(2,3)]
d10082 <- day0_nmds[grep('DA10082', day0_nmds$V2), c(2,3)]
d10093 <- day0_nmds[grep('DA10093', day0_nmds$V2), c(2,3)]
d10148 <- day0_nmds[grep('DA10148', day0_nmds$V2), c(2,3)]

plot(day0_nmds$axis1, day0_nmds$axis2, main="Similarity of day 0 communities, colored by donor")

points(d369, pch=16, col = 'red')
points(d430, pch=16, col = 'orange')
points(d431, pch=1, col = 'yellow')
points(d578, pch=16, col = 'green')
points(d581, pch=16, col = 'blue')
points(d884, pch=1, col = 'purple')
points(d953, pch=16, col = 'black')
points(d1134, pch=16, col = 'brown')
points(d1146, pch=16, col = 'aquamarine')
points(d1245, pch=1, col = 'cyan')
points(d1324, pch=16, col = 'darkgoldenrod')
points(d10027, pch=16, col = 'burlywood')
points(d10034, pch=1, col = 'coral1')
points(d10082, pch=1, col = 'darkolivegreen')
points(d10093, pch=16, col = 'darkgray')
points(d10148, pch=16, col = 'violetred2')

par(mar=c(5,5,5,5), xpd=TRUE)
legend <- c("d369", "d430", "d431", "d578", "d581", "d884", "d953", "d1134", "d1146", "d1245", "d1324", "d10027", "d10034", "d10082", "d10093", "d10148")
legend(x="topright", legend, col = c("red", "orange", "yellow", "green", "blue", "purple", "black", "brown", "aquamarine", "cyan", "darkgoldenrod", "burlywood", "coral1", "darkolivegreen", "darkgray", "violetred2"), 
        pch=16, bty='n', inset=c(-0.1,0), cex=0.9)

#Anova for comparing severe vs asymptomatic
#can do this in mothur very easily with distance matrix and design file 






