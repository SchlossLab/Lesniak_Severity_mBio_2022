################################################################################
#
# build_microbiotaSummaries.R
#
# This script builds a figure with microbiota summary graphs for each cage, 
# including OTUs belonging to Firmicutes, Bacteroidetes, Proteobacteria, and 
# Verrucomicrobia over time. 
# 
# 
#
# Dependencies...
#   * data/process/humanGF.final.an.unique_list.0.03.subsample.shared
#   * data/process/humanGF.final.an.unique_list.0.03.cons.taxonomy
#   * data/process/humanGF_metadata.txt
#
# Output...
#   * results/figures/CAGENAME_MicrobiotaLevel#_CUTOFFpercRA.tiff
#     where CAGENAME is the name of the cage, # is 0, 1, 2 for otu to family, 
#     and CUTOFF is the minimum relabundance for an OTU
#
################################################################################
setwd("~/Documents/Github/Schubert_humanCdGF_2015")

# read in the metadata file
all_metadata <- read.table(file="data/process/humanGF_metadata.txt", header=T, colClasses="character")
rownames(all_metadata) <- all_metadata$group

# SPECIFY THE LEVEL INTERESTED IN
# level 0 (otu), 1 (genus), 2 (family)
level <- 0

# read in the shared file and get the relative abundance
if(level==0){
  shared_file <- read.table(file="data/process/humanGF.final.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
} else if(level==1) {
  shared_file <- read.table(file="data/process/humanGF.final.tx.1.subsample.shared", header=T, row.names=2)
} else if(level==2) {
  shared_file <- read.table(file="data/process/humanGF.final.tx.2.subsample.shared", header=T, row.names=2)
}
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% all_metadata[all_metadata$sample_type=="stool",1])]
rel_abund <- rel_abund[overlap,]
metadata <- all_metadata[overlap,]

# let's get the relative abundances for OTUs that have an median 
# relative abundance over 1% within each cage
remove(max) # I accidentally changed mine

cage_med <- aggregate(rel_abund, by=list(metadata$cage_id), median)
cageIds <- cage_med[,1]
cage_med <- cage_med[,-1]
good_otus <- apply(cage_med, 2, max) > 1.0
rel_abund_good <- rel_abund[,good_otus]
numOtus <- sum(good_otus)

# get_taxonomy function returns a list containing 
# an array of phylum names and given taxonomy level names
# level is 0 (otu), 1 (genus), 2 (family)
get_taxonomy <- function(level){
  if(level==0){
    # let's get the taxonomy data for these OTUs
    taxonomy_file <- read.table(file="data/process/humanGF.final.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
    taxonomy <- taxonomy_file$Taxonomy
    taxonomy <- gsub("\\(\\d*\\)", "", taxonomy) #gets rid of bootstrap values
    taxonomy <- gsub(";unclassified", "", taxonomy) #gets rid of all unclassified
    taxonomy <- gsub("/.*", "", taxonomy) 
    taxonomy <- gsub(";$", "", taxonomy) #gets rid of extra ; at the end
    taxonomy <- gsub("\"$", "", taxonomy) #gets rid of end parenthesis
    taxonomy <- gsub(".*;", "", taxonomy) #eliminates all preceding text, except last classification
    taxonomy <- gsub("^\"", "", taxonomy) #gets rid of beginning parenthesis
    names(taxonomy) <- rownames(taxonomy_file)
    taxonomy <- taxonomy[colnames(rel_abund_good)] #only want the good phyla
  } else if(level==1){
    # let's get the taxonomy data for genera
    taxonomy_file <- read.table(file="data/process/humanGF.final.tx.1.cons.taxonomy", header=T, row.names=1)
    taxonomy <- taxonomy_file$Taxonomy
    taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
    taxonomy <- gsub("Bacteria;", "", taxonomy)
    taxonomy <- gsub(";$", "", taxonomy)
    # remove unclassified
    taxonomy <- gsub(";unclassified", "", taxonomy)
    taxonomy <- gsub("unclassified", "Bacteria", taxonomy)
    # then grab the last classification and clean up
    taxonomy <- gsub(".*;", "", taxonomy)
    taxonomy <- gsub("\"$", "", taxonomy) #gets rid of end quotes
    taxonomy <- gsub("^\"", "", taxonomy) #gets rid of beginning quotes
    names(taxonomy) <- rownames(taxonomy_file)
    taxonomy <- taxonomy[colnames(rel_abund_good)] #only want the good phyla
  } else if(level==2){
    # let's get the taxonomy data for families, edit these
    taxonomy_file <- read.table(file="data/process/humanGF.final.tx.2.cons.taxonomy", header=T, row.names=1)
    taxonomy <- taxonomy_file$Taxonomy
    taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
    taxonomy <- gsub("Bacteria;", "", taxonomy)
    taxonomy <- gsub(";$", "", taxonomy)
    # eliminate the end
    taxonomy <- gsub("[^;]*$", "", taxonomy)
    taxonomy <- gsub(";$", "", taxonomy)
    # then remove all unclassified
    taxonomy <- gsub(";unclassified", "", taxonomy)
    taxonomy <- gsub("unclassified", "Bacteria", taxonomy)
    # then grab the last and clean up
    taxonomy <- gsub(".*;", "", taxonomy)
    taxonomy <- gsub("\"$", "", taxonomy) #gets rid of end parenthesis
    taxonomy <- gsub("^\"", "", taxonomy) #gets rid of beginning parenthesis
    names(taxonomy) <- rownames(taxonomy_file)
    taxonomy <- taxonomy[colnames(rel_abund_good)] #only want the good phyla
  } else {print("Error: Level must be 0, 1, or 2")}
  
  # With the taxonomy file defined, pick out the corresponding phylum names
  phylTax <- taxonomy_file$Taxonomy
  phylTax <- gsub("\\(\\d*\\)", "", phylTax)
  phylTax <- gsub("Bacteria;", "", phylTax)
  phylTax <- gsub(";.*", "", phylTax)
  phylTax <- gsub("unclassified", "Bacteria", phylTax)
  phylTax <- gsub("\"$", "", phylTax) #gets rid of end parenthesis
  phylTax <- gsub("^\"", "", phylTax) #gets rid of beginning parenthesis
  names(phylTax) <- rownames(taxonomy_file)
  phylTax <- phylTax[colnames(rel_abund_good)] #only want the good phyla
  
  # Return a list with the phylum and taxonomy level labels interested in
  taxonomy <- list(phylTax, taxonomy)
}

# Let's get our taxonomy information
taxonomy <- get_taxonomy(level)
phylTax <- taxonomy[[1]]
otuTax <- taxonomy[[2]]


#let's get the median relative abundances and IQRs
med_ra <- aggregate(rel_abund_good, by=list(metadata$cage_id, metadata$day), median)  
med_ra <- t(med_ra)

l_qtr <- aggregate(rel_abund_good, by=list(metadata$cage_id, metadata$day), function(x){quantile(x, probs=c(0.25), na.rm=TRUE)})
l_qtr <- t(l_qtr)

u_qtr <- aggregate(rel_abund_good, by=list(metadata$cage_id, metadata$day), function(x){quantile(x, probs=c(0.75), na.rm=TRUE)})
u_qtr <- t(u_qtr)

# Graphing function
# cutoff is the median abundance % cutoff by cage
# level is 0 (otu), 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum)
time_otu_linechart_byPhylum <- function(cage, cutoff, level){
  
  # Which columns correspond to inputted cage parameter
  # Create that cage's own subset of data
  cageIndex <- which(med_ra[1,]==cage)
  med_ra_cage <- subset(med_ra, select = cageIndex)
  l_qtr_cage <- subset(l_qtr, select = cageIndex)
  u_qtr_cage <- subset(u_qtr, select = cageIndex)
  
  # What are the total number of phyla represented in this cage with OTUs greater than cutoff% abundance
  cage_numPhyla <- as.numeric(apply(med_ra_cage[3:(numOtus+2),], 1, max)) > cutoff
  cage_phyla <- unique(phylTax[cage_numPhyla])
  cage_numPhyla <- length(cage_phyla)
  
  png(height=10, width=6, file=paste0("results/figures/", cage, "_MicrobiotaLevel", level, "_", cutoff, "percRA.png"), unit="in", res=300)
  
  layoutMatrix <- matrix( c(  seq(1, cage_numPhyla+1)), byrow=F, ncol=1)
  layoutMatrix <- cbind( matrix(c(rep(cage_numPhyla+2, cage_numPhyla+1)), ncol=1), layoutMatrix)
  
  z <- layout(
    layoutMatrix, 
    widths=c(.07,  1), heights=c(rep(1,cage_numPhyla), 0.3)
  ) 
  #layout.show(z)
  
  #want to graph by each phylum in its own graph
  #only showing the OTUs in that cage with greater than 1.0
  #so that the graph is hopefully not too cluttered
  for(i in 1:cage_numPhyla){
    
    # Which indices correspond to phylum i for the given med_ra_cage
    phylInd <- which(phylTax==cage_phyla[i]) 
    
    ###! Might consider removing this filter to keep colors/symbols
    ###! consistent across each cage
    # For phylum i, which OTUs over the time course still
    # have a median relative abundance greater than cutoff% for this cage
    if(length(phylInd)==1){
      good_phyl_otus <- as.numeric(max(med_ra_cage[(phylInd+2),])) > cutoff
      names(good_phyl_otus) <- row.names(med_ra_cage)[phylInd+2]
      goodOtuInd <- which(good_phyl_otus)
      num_goodOtus <- sum(good_phyl_otus)
    } else{
      good_phyl_otus <- as.numeric(apply(med_ra_cage[(phylInd+2),], 1, max)) > cutoff
      names(good_phyl_otus) <- row.names(med_ra_cage[(phylInd+2),])
      goodOtuInd <- which(good_phyl_otus)
      num_goodOtus <- sum(good_phyl_otus)
    }
        
    days <- med_ra_cage[2, ]
    days <- as.numeric(days)
   
    # Not looking at negative days right now 
    neg_days <- which(days < 0)
    if(length(neg_days)>0) {
      days <- days[-neg_days]
      cageIndex <- cageIndex[-neg_days]
    }
    
    # What is the maximum relabund (upper limit) for any OTU belonging to phylum i
    maxVal <- max(as.numeric(u_qtr_cage[phylInd[goodOtuInd]+2,]))
      
    
    # Plot a line chart showing relative abundance over time
    # Start with the box and axes
    par(mar=c(0, 1, 1, 10)+.1, xpd=FALSE)
    plot(c(0,10), c(0, ceiling(maxVal)), las=1, type="n", xaxt="n", ylab="")
    axis(1, at=c(0:10), label=c(rep("", 11))) 
    abline(h=seq(0, 100, by=5), col="light gray", lty="dashed")
    
    # If the last graph, add the xvalues on the axes
    if(i==cage_numPhyla){
      mtext(c(0:10), side=1, at=c(0:10), padj=1)      
    }   
   
    # Customize graphs for different taxonomy levels
    if(level==0){
      pch <- gsub("Otu0*", "", names(goodOtuInd))
      legend_names <- paste0("OTU", pch, ": ", otuTax[names(goodOtuInd)])
      level_name <- "OTUs"
    } else{
      pch <- c(1:num_goodOtus)
      legend_names <- otuTax[names(goodOtuInd)]
      level_name <- ifelse(level==1, "Genus Phylotypes", "Family Phylotypes")
    }
    colors <- rainbow(num_goodOtus, v=.75)
    
    # Graphing each line
    dayOrder <- order(days)
    xvalues <- sort(days)
    for(j in 1:num_goodOtus){
      yvalues <- med_ra_cage[phylInd[goodOtuInd[j]]+2,]
      yvalues <- as.numeric(yvalues[dayOrder])
      upper_values <- as.numeric(u_qtr_cage[phylInd[goodOtuInd[j]]+2,])
      upper_values <- upper_values[dayOrder]
      lower_values <- as.numeric(l_qtr_cage[phylInd[goodOtuInd[j]]+2,])
      lower_values <- lower_values[dayOrder]
      arrows(x0=xvalues, x1=xvalues, y0=yvalues, y1=upper_values, angle=90, length=0.05, col=colors[j])
      arrows(x0=xvalues, x1=xvalues, y0=yvalues, y1=lower_values, angle=90, length=0.05, col=colors[j])
      if(level==0){
        lines(xvalues, yvalues, lty=1, col=colors[j], type="o", cex=1.3, pch=16,  lwd=1)
        text(xvalues, yvalues, pch[j], col="black", cex=.6)
      } else{
        lines(xvalues, yvalues, lty=1, col=colors[j], type="o", pch[j], lwd=1, cex=1)
      }
    }
    par(xpd=TRUE)
    if(level==0){
      legend(10.5, ceiling(maxVal), legend_names, col=colors, lwd=1.5, lty=1,  cex=.7, pt.cex=1.5, title=paste0(cage_phyla[i], ": ", level_name))    
    } else{
      legend(10.5, ceiling(maxVal), legend_names, col=colors, pch=pch, lty=1,  cex=.7, pt.cex=1.5, title=paste0(cage_phyla[i], ": ", level_name))    
    }
  }
  # Add remaining labels in their own space on the layout
  par(mar=c(0, 0, 0, 10)+.1)
  plot.new()
  mtext("Time (Days)", side=1, line=-1.2, cex=1.2)
  par(mar=c(0,0,0,0)+.1)
  plot.new()
  mtext("% Relative Abundance", side=2, line=-1.2, cex=1.2)
  
  dev.off()
}

# make individual cage summary graphs
inoculaInd <- which(cageIds=="inoculum")
cageIds <- cageIds[-inoculaInd]
numEndCages <- length(cageIds)
for(i in 1:numEndCages){
  time_otu_linechart_byPhylum(cageIds[i], 5, level)
}

# Testing: 
#cage <- "NP1"
#cutoff <- 1
#time_otu_linechart_byPhylum(cage, cutoff, level)
