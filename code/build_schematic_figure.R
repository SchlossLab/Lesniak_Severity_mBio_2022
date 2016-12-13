# Builds figure 1
# Experimental timeline code from Matt Jenior

# Load dependencies 
deps <- c('shape', 'plotrix', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Documents/Schloss_Lab/Schubert_humanCdGF_XXXX_2016/results/figures/figure_1.pdf'
pdf(file=plot_file, width=9, height=9)

# Create layout for multi-plot
layout(mat=matrix(c(1,
                    2,
                    3,
                    4,
                    5), nrow=5, ncol=1, byrow=TRUE))

# Figure 1A - experimental timeline code 
#-------------------------------------------------------------------------------------------------------------------------------------#

#---------Option 1, text based timeline---------#
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot
Arrows(x0=-4, y0=0, x1=3, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.3)
segments(x0=c(-4,-0.5,2), y0=c(0.5,0.5,0.5), x1=c(-4,-0.5,2), y1=c(-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-3.5,-3,-2.5,-2,-1.5, -1), y0=c(0.25,0.25,0.25,0.25,0.25), x1=c(-3.5,-3,-2.5,-2,-1.5, -1), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
segments(x0=c(0,0.5,1,1.5), y0=c(0.25,0.25,0.25,0.25,0.25), x1=c(0,0.5,1,1.5), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
text(x=c(-4,-0.5, 2), y=c(-0.8,-0.8, -0.8), c('Day -14', 'Day 0', 'Day 10'), cex=1)
text(x=c(-3.35,-0.5, 2), y=c(0.7,0.7, 0.7), c('Gavage human stool', 'Gavage C. difficile', 'Euthanize'), cex=1)
text(x=-4.5, y=0, 'A', cex=1.5, font=2)

#---------Option 2, text based timeline---------#


#------------Option 3, rectangles timeline------#

#-----------Option 4, mouse timeline--------#

#----------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------#

#Figure 1B - nmds of donor stool prior to inoculation 

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

#add back label and numOTUs column to subsetted shared file-- need this to run dist.shared
inocula_shared <- cbind(label=0.03, inocula_shared[1], numOtus=5671, inocula_shared[2:ncol(inocula_shared)])
colnames(inocula_shared)[2] <- 'Group'

#write file out so we can use it in mothur
write.table(inocula_shared, file='data/process/inocula.shared', quote=F,sep='\t',row.names=F)

#run dist.shared with iters and subsamp in TYC form, nmds in mothur to get axes file
#mothur commands
#dist.shared(shared=inocula.shared, calc=thetayc, subsample=2000)
#nmds(phylip=inocula.thetayc.0.03.lt.dist)
#mothur NMDS output:
# Number of dimensions:	2
# Lowest stress :	0.318153
# R-squared for configuration:	0.458028

inoc_nmds <- read.table(file='data/process/inocula.thetayc.0.03.lt.nmds.axes', header = T)
#merge design and nmds files to prep for plotting 
inocula_design <- read.table(file='data/raw/donor.design', header = F)
inocula_meta_design <- read.table(file='data/raw/donor_meta.design', header = F)
inocula_meta_nmds <- merge(inoc_nmds, inocula_meta_design, by.x='group', by.y='V1')
colnames(inocula_meta_nmds)[4] <- "Disease_status"

#build NMDS plot
plot(inocula_meta_nmds$axis1, inocula_meta_nmds$axis2, xlab="NMDS Axis 1", ylab= "NMDS Axis 2")
no_d <- inocula_meta_nmds[grep('No_diarrhea', inocula_meta_nmds$Disease_status), c(2,3)]
diarrhea <- inocula_meta_nmds[grep('Diarrhea', inocula_meta_nmds$Disease_status), c(2,3)]
cdiff <- inocula_meta_nmds[grep('C._difficile', inocula_meta_nmds$Disease_status), c(2,3)]

points(no_d, pch=16, col = "blue")
points(diarrhea, pch=16, col = "yellow2")
points(cdiff, pch=16, col = "black")

legend <- c("No diarrhea", "Diarrhea", "C. difficile")
legend(x="topright", legend, col = c("blue", "yellow2", "black"), pch=16, cex=0.8)

#put labels outside of plot
#text(x=-0.3, y=0.7, 'B', cex=1.5, font=2)
#will need to do that ADONIS here for stats

#-----------------------------------------------------------------------------------------------------#
#Figure 1C - NMDS of mice/cages on day 1 

#prepare files and run on mothur
#subset shared to be day 0 only 
day0 <- rownames(meta_file)[meta_file$day==0]
day0_shared <- full_shared[full_shared$sample_id %in% day0,]
day0_shared <- cbind(label=0.03, day0_shared[1], numOtus=5671, day0_shared[2:ncol(day0_shared)])
colnames(day0_shared)[2] <- 'Group'
write.table(day0_shared, file='data/process/day0.shared', quote=F,sep='\t',row.names=F)

#run in mothur
#dist.shared(shared=day0.shared, calc=thetayc, subsample=2000)
#nmds(phylip=day0.thetayc.0.03.lt.dist)
#mothur NMDS output:
#Number of dimensions:	2
#Lowest stress :	0.35381
#R-squared for configuration:	0.385726
day0_nmds <- read.table(file='data/process/day0.thetayc.0.03.lt.nmds.axes', header = T)
#merge with design files to get ready to plot 
d0_donor <- read.table(file='data/raw/sample_donor.design', header = F)
day0_donor_nmds <- merge(day0_nmds, d0_donor, by.x='group', by.y='V1')
d0_outcome <- read.table(file='data/raw/d0_only_outcome.design', header=F)
day0_outcome_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')
#day0_outcome_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')
day0_donor_nmds <- merge(day0_donor_nmds, d0_outcome, by.x='group', by.y='V1')
names(day0_donor_nmds)[4] <- "Donor"
names(day0_donor_nmds)[5] <- "Outcome"

#plot
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
            "maroon","orchid1")
plot(day0_donor_nmds$axis1, day0_donor_nmds$axis2, xlab='NMDS Axis 1', ylab='NMDS Axis 2', col=colors[as.numeric(day0_donor_nmds$Donor)], pch=c(16,17)[as.numeric(day0_donor_nmds$Outcome)])

par(mar=c(5,5,5,5), xpd=TRUE)
#just make the legend mild and severe 
legend <- c("Mild", "Severe")
legend(x="topright", legend, pch=c(16,17), cex=0.8)


#now sew together plots into one figure and save as pdf 

#Clean up
dev.off()
rm(plot_file)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()