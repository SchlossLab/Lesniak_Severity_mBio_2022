################################################################################
#
# build_CageSummaries.R
#
# This script builds a figure with 4 summary graphs for each cage, including
# weight change over time, *toxin activity detected over time*, C. difficile 
# levels observed, phylum changes over time. 
# 
# * Toxin graph to be added
# * this is Alyx's old code, I'm going to try using it instead of my phyla over time graph. let's see...
#
# Dependencies...
#     * data/process/human_CdGF_metadata.txt
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#
#
# Output...
#   * results/figures/CAGENAME_summary.tiff
#
################################################################################

# read in the metadata file
meta_file <- read.table("data/process/human_CdGF_metadata.txt", header = TRUE, sep='\t', fill = TRUE, row.names=3)
numCages <- length(unique(all_metadata$cage_id))
cageIds <- unique(all_metadata$cage_id)
rownames(all_metadata) <- all_metadata$group

# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/humanGF.final.tx.5.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% all_metadata[all_metadata$sample_type=="stool",1])]
rel_abund <- rel_abund[overlap,]
end_metadata <- all_metadata[overlap,]
end_cageIds <- unique(end_metadata$cage_id)


# let's get the relative abundances for those phyla that have at least one
# sample where they are more than 10% o the community
remove(max) #i accidentally changed mine
good_phyla <- apply(rel_abund, 2, max) > 10
rel_abund_good <- rel_abund[,good_phyla]
numTaxa <- sum(good_phyla)

# let's get the taxonomy data for phyla
taxonomy_file <- read.table(file="data/process/humanGF.final.tx.2.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub("Bacteria;", "", taxonomy)
taxonomy <- gsub(";.*", "", taxonomy)
taxonomy <- gsub("unclassified", "Bacteria", taxonomy)
taxonomy <- gsub("\"$", "", taxonomy) #gets rid of end parenthesis
taxonomy <- gsub("^\"", "", taxonomy) #gets rid of beginning parenthesis
names(taxonomy) <- rownames(taxonomy_file)
taxonomy <- taxonomy[colnames(rel_abund_good)] #only want the good phyla


#let's rename the OTUs
colnames(rel_abund_good) <- taxonomy


#Individual graph functions
time_phylum_linechart <- function(cage){
  #let's get the median relative abundances and IQRs
  med_ra <- aggregate(rel_abund_good, by=list(end_metadata$cage_id, end_metadata$day), median)
  med_ra <- t(med_ra)
  
  l_qtr <- aggregate(rel_abund_good, by=list(end_metadata$cage_id, end_metadata$day), function(x){quantile(x, probs=c(0.25), na.rm=TRUE)})
  l_qtr <- t(l_qtr)
  
  u_qtr <- aggregate(rel_abund_good, by=list(end_metadata$cage_id, end_metadata$day), function(x){quantile(x, probs=c(0.75), na.rm=TRUE)})
  u_qtr <- t(u_qtr)
  
  # Which indices correspond to cage?
  cageIndex <- which(med_ra[1,]==cage)
  
  # Get the days with data
  days <- med_ra[2, cageIndex]
  days <- as.numeric(days)
  
  # Don't want negative days right now
  neg_days <- which(days < 0)
  if(length(neg_days)>0) {
    days <- days[-neg_days]
    cageIndex <- cageIndex[-neg_days]
  }
  dayOrder <- order(days)
  
  # Plot + Parameters
  plot(c(0,10), c(0, 100), type="n", xaxt="n", yaxt="n", xlab="", ylab="% Relative Abundance")
  axis(1, las=1, at=c(0:10), label=c(0:10), line=0) 
 # mtext("Time (days)", side=1, line=1, cex=.8)
  axis(2, las=1, at=seq(0, 100, by=20), label=seq(0, 100, by=20))
  abline(h=seq(0, 100, by=10), col="light gray", lty="longdash")
  colors <- rainbow(numTaxa, v=.75)
 
  # Plotting each line
  xvalues <- sort(days)
  for(i in 1:numTaxa){
    yvalues <- med_ra[2+i, cageIndex]
    yvalues <- as.numeric(yvalues[dayOrder])
    upper_values <- as.numeric(u_qtr[2+i,cageIndex])
    upper_values <- upper_values[dayOrder]
    lower_values <- as.numeric(l_qtr[2+i,cageIndex])
    lower_values <- lower_values[dayOrder]
    lines(xvalues, yvalues, lty=1, col=colors[i], type="o")
    arrows(x0=xvalues, x1=xvalues, y0=yvalues, y1=upper_values, angle=90, length=0.1, col=colors[i])
    arrows(x0=xvalues, x1=xvalues, y0=yvalues, y1=lower_values, angle=90, length=0.1, col=colors[i])
  }
  legend("top", colnames(rel_abund_good), col=colors, pch="o", lty=1, horiz=TRUE,  cex=.4)
}

time_weight_linechart <- function(cage){
  cageIndex <- which(end_metadata$cage_id==cage & end_metadata$day>=0)
  yvalues <- as.numeric(end_metadata[cageIndex, "percent_weightLoss"])
  xvalues <- as.numeric(end_metadata[cageIndex, "day"])
  donor <- unique(end_metadata[cageIndex, "human_source"])
  
  #get median weights and upper and lower quartiles for percent weight loss
  med_w <- aggregate(yvalues, by=list(xvalues), median)
  med_w <- t(med_w)
  
  l_qtr <- aggregate(yvalues, by=list(xvalues), function(x){quantile(x, probs=c(0.25), na.rm=TRUE)})
  l_qtr <- t(l_qtr)
  
  u_qtr <- aggregate(yvalues, by=list(xvalues), function(x){quantile(x, probs=c(0.75), na.rm=TRUE)})
  u_qtr <- t(u_qtr)
  
  #plot
  plot(c(0,10), c(75, 125), type="n", xaxt= "n", xaxt="n", yaxt="n", ylab="% Weight Loss")
  axis(1, las=1, at=c(0:10), label=rep("", 11))
  axis(2, las=1, at=seq(65, 125, by = 10), label=seq(65, 125, by = 10))
  abline(h=seq(75, 125, by=5), col="gray", lty="dashed")
  lines(med_w[1,], med_w[2,], lty=1, type="o" )
  arrows(x0=med_w[1,], x1=med_w[1,], y0=med_w[2,], y1=u_qtr[2,], angle=90, length=0.1)
  arrows(x0=med_w[1,], x1=med_w[1,], y0=med_w[2,], y1=l_qtr[2,], angle=90, length=0.1)
  text(0, 120, labels=donor, cex=2, pos=4)
  
}

time_cdiff_linechart <- function(cage) {
  cageIndex <- which(end_metadata$cage_id==cage & end_metadata$day>=0)
  yvalues <- as.numeric(end_metadata[cageIndex, "cdiff_cfu"])
  xvalues <- as.numeric(end_metadata[cageIndex, "day"])
  cdiff_strain <- unique(end_metadata[cageIndex, "cdiff_strain"])
  
  #get median weights and upper and lower quartiles for percent weight loss
  med_cd <- aggregate(yvalues, by=list(xvalues), median)
  med_cd <- t(med_cd)
  
  l_qtr <- aggregate(yvalues, by=list(xvalues), function(x){quantile(x, probs=c(0.25), na.rm=TRUE)})
  l_qtr <- t(l_qtr)
  
  u_qtr <- aggregate(yvalues, by=list(xvalues), function(x){quantile(x, probs=c(0.75), na.rm=TRUE)})
  u_qtr <- t(u_qtr)
  
  #plot
  plot(c(0,10), c(1, 10e10), log="y", type="n", xaxt= "n", yaxt="n", ylab="C. difficile CFU/g Feces", xlab="")
  axis(1, las=1, at=c(0:10), label=c(rep("", 11)))
  axis(2, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8, 1e10), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8), expression(10^10)))      
  abline(h=c(1, 1e2, 1e4, 1e6, 1e8, 1e10), col="gray", lty="longdash")
  abline(h=c(10, 1e3, 1e5, 1e7, 1e9 ), col="light gray", lty="dashed")
  
  for(i in 1:dim(med_cd)[2]){
    if(is.na(med_cd[2,i])){next}
    if(med_cd[2,i]==0){
      med_cd[2,i] <- 1
    } 
  }
  lines(med_cd[1,], med_cd[2,], lty=1, type="o" )
  arrows(x0=med_cd[1,], x1=med_cd[1,], y0=med_cd[2,], y1=u_qtr[2,], angle=90, length=0.1)
  arrows(x0=med_cd[1,], x1=med_cd[1,], y0=med_cd[2,], y1=l_qtr[2,], angle=90, length=0.1)
  text(0, 1e10, labels=cdiff_strain, cex=2, pos=4)
}

time_toxin <- function(cage){
  cageIndex <- which(end_metadata$cage_id==cage & end_metadata$day>=0)
  yvalues <- as.numeric(end_metadata[cageIndex, "toxin_log"])
  xvalues <- as.numeric(end_metadata[cageIndex, "day"])
  
  #get median weights and upper and lower quartiles for percent weight loss
  med_t <- aggregate(yvalues, by=list(xvalues), median)
  med_t <- t(med_t)
  
  l_qtr <- aggregate(yvalues, by=list(xvalues), function(x){quantile(x, probs=c(0.25), na.rm=TRUE)})
  l_qtr <- t(l_qtr)
  
  u_qtr <- aggregate(yvalues, by=list(xvalues), function(x){quantile(x, probs=c(0.75), na.rm=TRUE)})
  u_qtr <- t(u_qtr)
  
  #plot
  plot(c(0,10), c(1, 7), type="n", xaxt= "n", xaxt="n", yaxt="n", ylab="Log Reciprocal Dilution", xlab="")
  axis(1, las=1, at=c(0:10), label=rep("", 11))
  axis(2, las=1, at=c(1:7), label=c(1:7))
  abline(h=c(1:7), col="gray", lty="dashed")
  lines(med_t[1,], med_t[2,], lty=1, type="o" )
  arrows(x0=med_t[1,], x1=med_t[1,], y0=med_t[2,], y1=u_qtr[2,], angle=90, length=0.1)
  arrows(x0=med_t[1,], x1=med_t[1,], y0=med_t[2,], y1=l_qtr[2,], angle=90, length=0.1)
  
}

build_cage_summary_graphs <- function(cage){
  
  png(height=6, width=4, file=paste0("results/figures/", cage, "_summary.png"), unit="in", res=300)
  
  z <- layout(
    matrix( c(  seq(1, 5)), byrow=F, ncol=1), 
    widths=c(1), heights=c(rep(1,4), 0.3)
  ) #there are 4 graphs
  
  par(mar=c(0.5, 5, 1.25, 0.5))
  
  time_weight_linechart(cage)
  time_toxin(cage)
  time_cdiff_linechart(cage)
  time_phylum_linechart(cage)
  plot.new()
  mtext("Time (Days)", cex=.8, side=1, line=-1)
  
  dev.off()
}

#make individual cage summary graphs
numCages <- length(end_cageIds)
for(i in 2:numCages){
  build_cage_summary_graphs(end_cageIds[i])
}

build_cage_summary_graphs(cage)
