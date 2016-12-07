###################
#
# Inocula_nmds.R
#
#Create shared file of inocula only and plot NMDS and color points by outcome
#
#    
#
###################
install.packages("rgl")
library("rgl")


#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header=T)
full_shared <- read.table('data/process/humanCdGF_full.shared', sep='\t', header = T)

#shared files are processed in a make_meta script to remove erroneous samples 
#so we need format processed shared to look like a true shared file:
full_shared_group <- cbind(label=0.03, full_shared[1], numOtus=5671, full_shared[2:ncol(full_shared)])
colnames(full_shared_group)[2] <- 'Group'

#subset shared to be inoculum only
inocula <- rownames(meta_file)[meta_file$cage_id=='inoculum']
inocula_shared <- full_shared[full_shared$sample_id %in% inocula,]

#subset shared to be day 0 only 
day0 <- rownames(meta_file)[meta_file$day==0]
day0_shared <- full_shared[full_shared$sample_id %in% day0,]

#add back label and numOTUs column to subsetted shared file-- need this to run dist.shared
inocula_shared <- cbind(label=0.03, inocula_shared[1], numOtus=5671, inocula_shared[2:ncol(inocula_shared)])
colnames(inocula_shared)[2] <- 'Group'

#for day 0 also.. probably could make this a function
day0_shared <- cbind(label=0.03, day0_shared[1], numOtus=5671, day0_shared[2:ncol(day0_shared)])
colnames(day0_shared)[2] <- 'Group'

#combine these all to run as one file duhhhh
inputd0_shared <- rbind(inocula_shared, day0_shared)

#write file out so we can use it in mothur
write.table(inocula_shared, file='data/process/inocula.shared', quote=F,sep='\t',row.names=F)
write.table(day0_shared, file='data/process/day0.shared', quote=F,sep='\t',row.names=F)
write.table(inputd0_shared, file='data/process/inputd0.shared', quote=F, sep='\t', row.names=F)
write.table(full_shared_group, file='data/process/full_shared.shared', quote=F, sep='\t', row.names=F)

#run dist.shared with iters and subsamp in TYC form, nmds in mothur to get axes file
#OLD mothur ouput stats for OLD SUBSAMP combined file: 
#Number of dimensions:	2
#Lowest stress :	0.362683
#R-squared for configuration:	0.380491

#mothur output stats for new combined file
#for combo
#Number of dimensions:	2
#Lowest stress :	0.364672
#R-squared for configuration:	0.372893

#for day 0 
#Number of dimensions:	2
#Lowest stress :	0.348974
#R-squared for configuration:	0.407514

combo_nmds <- read.table(file='data/process/inputd0.thetayc.0.03.lt.nmds.axes', header = T)

#need to get these other nmds's done to get the right plot 
inoc_nmds <- read.table(file='data/process/inocula.thetayc.0.03.lt.nmds.axes', header = T)
day0_nmds <- read.table(file='data/process/day0.thetayc.0.03.lt.nmds.axes', header = T)

day0_3D <- read.table(file='data/process/day0.thetayc.3Dnmds.axes', header = T)


#merge design and nmds files to prep for plotting 
inocula_design <- read.table(file='data/raw/donor.design', header = F)
inocula_nmds <- merge(inoc_nmds, inocula_design, by.x='group', by.y='V1')

inocula_meta_design <- read.table(file='data/raw/donor_meta.design', header = F)
inocula_meta_nmds <- merge(inoc_nmds, inocula_meta_design, by.x='group', by.y='V1')
colnames(inocula_meta_nmds)[5] <- "Disease_status"

combo_design <- read.table(file='data/raw/d0_outcome.design', header = F)
combo_nmds <- merge(combo_nmds, combo_design, by.x='group', by.y='V1')

#need 3D version of these 
d0_donor <- read.table(file='data/raw/sample_donor.design', header = F)
day0_donor_nmds <- merge(day0_nmds, d0_donor, by.x='group', by.y='V1')
day0_3D_nmds <- merge(day0_3D, d0_donor, by.x='group', by.y='V1')

d0_outcome <- read.table(file='data/raw/d0_only_outcome.design', header=F)
day0_outcome_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')
day0_outcome_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')

day0_donor_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')

#use this one
day0_3D_outcome <- merge(day0_3D_nmds, d0_outcome, by.x='group', by.y='V1')

#make nmds inocula plot
severe <- inocula_nmds[grep('Severe', inocula_nmds$V2), c(2,3)]
asymptomatic <- inocula_nmds[grep('Asymptomatic', inocula_nmds$V2), c(2,3)]
plot(inocula_nmds$axis1, inocula_nmds$axis2, main="Similarity of donor inocula")
points(severe, pch=16, col = "red")
points(asymptomatic, pch=16, col = "black")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)
#amova in mothur output for asymptomatic vs severe DONORS: p-value: 0.508

#nmds inocula plot with clinical metadata, use this one, except need just 2D nmds file***

diarrhea <- inocula_meta_nmds[grep('Diarrhea')]
plot(inocula_meta_nmds$axis1, inocula_meta_nmds$axis2)

pdf("fig1B_actual.pdf", width = 8, height = 5)

ggplot(inocula_meta_nmds, aes(axis1, axis2, group = Disease_status, color = Disease_status)) + geom_point(size = 3) + 
  theme_bw() + xlab("NMDS Axis 1") + ylab("NMDS Axis 2") + theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linetype='solid')) +
  theme(legend.background = element_rect(color='black', size=.5, linetype="solid")) + theme(legend.key = element_blank()) +
  scale_fill_discrete(name="Clinical status", breaks=c("C._difficile", "Diarrhea", "No_diarrhea"), labels=c("C. difficile", "Diarrhea", "Healthy"))

dev.off()

plot3d(inocula_meta_nmds$axis1, inocula_meta_nmds$axis2, inocula_meta_nmds$axis3, type = "p", col = as.numeric(inocula_meta_nmds$Disease_status), size = 5, xlab='NMDS axis 1', ylab= 'NMDS axis 2', zlab= 'NMDS axis 3')
legend3d("topleft", legend= paste(c("C. difficile", "No diarrhea", "Diarrhea")), pch=16, col=c('black', 'green', 'red'), inset=c(0.08))

plot3d(day0_3D_outcome$axis1, day0_3D_outcome$axis2, day0_3D_outcome$axis3, type = "p", col = as.numeric(day0_3D_outcome$V2.x), pch= as.numeric(day0_3D_outcome$V2.y), size = 5)
legend3d("topright", legend= unique(day0_3D_outcome$V2.x), pch=16, col=unique(as.numeric(day0_3D_outcome$V2.x)), inset=c(0.08))

#donor and day 0 by outcome plot
severe2 <- combo_nmds[grep('Severe', combo_nmds$V2), c(2,3)]
asymp2 <- combo_nmds[grep('Asymptomatic', combo_nmds$V2), c(2,3)]
plot(combo_nmds$axis1, combo_nmds$axis2, main="Similarity of donor and day0 by outcome")
points(severe2, pch=16, col = "red")
points(asymp2, pch=16, col = "black")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)

#day 0 by outcome plot, no donors. 
#amova output for that is P <0.001! 

severe4 <- day0_outcome_nmds[grep('Severe', day0_outcome_nmds$V2), c(2,3)]
asymp4 <- day0_outcome_nmds[grep('Asymptomatic', day0_outcome_nmds$V2), c(2,3)]
plot(day0_outcome_nmds$axis1, day0_outcome_nmds$axis2, main="Similarity day0 by outcome", col = day0_outcome_nmds$V2, pch=as.numeric(day0_outcome_nmds$V2))
points(severe4, pch=16, col = "red")
points(asymp4, pch=16, col = "black")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)



#nmds colored by donor plot
library("RColorBrewer")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


day0_donor_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')
names(day0_donor_nmds)[4] <- "Donor"
names(day0_donor_nmds)[5] <- "Outcome"

plot(day0_donor_nmds$axis1, day0_donor_nmds$axis2, main="Similarity of day 0 communities, colored by donor", col = day0_donor_nmds$V2, pch=16)

colors <- c("dodgerblue2","#E31A1C", # red
             "green4",
             "#6A3D9A", # purple
             "#FF7F00", # orange
             "black","gold1",
             "skyblue2","#FB9A99", # lt pink
             "palegreen2",
             "#CAB2D6", # lt purple
             "#FDBF6F", # lt orange
             "gray70", "khaki2",
             "maroon","orchid1","deeppink1","blue1","steelblue4",
             "darkturquoise","green1","yellow4","yellow3",
             "darkorange4","brown”")


#use this plot
pdf("fig1B.pdf", width = 13, height = 5)

ggplot(day0_donor_nmds, aes(axis1, axis2, group = Donor, color = Donor, shape = Outcome)) + geom_point(aes(fill=Donor), size = 3) +
   theme_bw() + ggtitle("Similarity of day 0 communities") + scale_shape_discrete(name = "Outcome") + labs(x='NDMS Axis 1', y = 'NDMS axis 2') + scale_color_manual(values = colors)




### Plotting function to plot convex hulls
### Filename: Plot_ConvexHull.R
### Notes:
############################################################################

# INPUTS:
# xcoords: x-coordinates of point data
# ycoords: y-coordinates of point data
# lcolor: line color

# OUTPUTS:
# convex hull around data points in a particular color (specified by lcolor)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
} 


store <- Plot_ConvexHull(xcoord=day0_donor_nmds$axis1, ycoord = day0_donor_nmds$axis2, lcolor=day0_donor_nmds$V2.x)





#geom_polygon(aes(fill = V2.x), alpha = 0.3) makes shitty polygons
#stat_ellipse()

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





points(d369, pch=16, col = 'red')
points(d430, pch=16, col = 'orange')
points(d431, pch=1, col = 'yellow')
points(d578, pch=16, col = 'green')
points(d581, pch=16, col = 'blue')
points(d884, pch=1, col = 'purple')
points(d953, pch=16, col = 'black')
points(d1134, pch=16, col = 'brown')
points(d1146, pch=16, col = 'aquamarine')
points(d1245, pch=15, col = 'cyan')
points(d1324, pch=16, col = 'darkgoldenrod')
points(d10027, pch=16, col = 'burlywood')
points(d10034, pch=1, col = 'coral1')
points(d10082, pch=1, col = 'darkolivegreen')
points(d10093, pch=16, col = 'darkgray')
points(d10148, pch=16, col = 'violetred2')

par(mar=c(5,5,5,5), xpd=TRUE)
legend <- c("d369", "d430", "d431", "d578", "d581", "d884", "d953", "d1134", "d1146", "d1245", "d1324", "d10027", "d10034", "d10082", "d10093", "d10148")
legend(x="topright", legend, col = c("red", "orange", "yellow", "green", "blue", "purple", "black", "brown", "aquamarine", "cyan", "darkgoldenrod", "burlywood", "coral1", "darkolivegreen", "darkgray", "violetred2"), 
        pch=16, bty='n', inset=c(-0.07,0), cex=0.9)

#Anova for comparing severe vs asymptomatic
#can do this in mothur very easily with distance matrix and design file 


#Dendrogram code
install.packages("ape")

library(ape)

inoc_tree <- read.tree(file='data/process/inocula_2194.thetayc.0.03.tre')
day0_tree <- read.tree(file='data/process/day0.thetayc.0.03.ave.tre')

#for plotting donors by outcome 
dead <- inoc_tree$tip.label[c(3,6,7,8,15)]
plot(inoc_tree, type='radial', tip.color=ifelse(inoc_tree$tip.label %in% dead, 'red', 'black'))
plot(inoc_tree, tip.color=ifelse(inoc_tree$tip.label %in% dead, 'red', 'black'))

#for plotting donors by initial disease state
diah <- inoc_tree$tip.label[c(2,5,7,11)]
cdiff <- inoc_tree$tip.label[8]
plot(inoc_tree, type='radial', tip.color=ifelse(inoc_tree$tip.label %in% diah, 'blue', ifelse(inoc_tree$tip.label %in% cdiff, 'red', 'black')))
plot(inoc_tree, tip.color=ifelse(inoc_tree$tip.label %in% diah, 'blue', ifelse(inoc_tree$tip.label %in% cdiff, 'red', 'black')))
legend("bottomleft", c('diarrhea', 'cdiff infected', 'no diarrhea'), col = c('blue', 'red', 'black'), pch =16)

#for plotting mice at day 0 by outcome 

mild <- subset(d0_outcome, d0_outcome$V2 == 'Asymptomatic')
milder <- mild[1]

severe9 <- subset(d0_outcome, d0_outcome$V2 == 'Severe')
severer <- severe9[1]

#day0 mouse by outcome
plot(day0_tree, tip.color=ifelse(day0_tree$tip.label %in% severer$V1, 'red', 'black'),  cex=0.7, main = "Day 0 dendrogram")
legend("bottomleft", c("mild", "severe"), col = c('black', 'red'), pch = 16)


