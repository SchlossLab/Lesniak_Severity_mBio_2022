################################################################################
#
# build_randomForest.R
#
# This runs and plots randomForest results.
# Want to predict: 
#   * Day 1 cdifficile levels based on day 0 community
#     - Maybe a function to generally take the previous day and subsequent 
#       day's cdiff predictions, or input the y value
#   * Level of toxin based on community
#   * Community to predict death or not--try using just last day of each 
#     or include each day with its outcome of death or not
#     - Categorical random forest to predict low colonization, high/sustained 
#       colonization, or death
#
# Dependencies...
#   * data/process/humanGF.final.an.unique_list.0.03.subsample.shared
#   * data/process/humanGF.final.an.unique_list.0.03.cons.taxonomy
#   * data/process/humanGF_metadata.txt
#
# Output...
#   * 
#
################################################################################
setwd("~/Documents/Github/Schubert_humanCdGF_2015")

# Load the randomForest package
load_package <- function(package){
  if(!(package %in% rownames(installed.packages()))){
    install.packages(package)
  }
  library(package, quietly=TRUE, character.only=TRUE)
}

load_package("randomForest")



# read in the metadata file
all_metadata <- read.table(file="data/process/humanGF_metadata.txt", header=T)
rownames(all_metadata) <- all_metadata$group


# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/humanGF.final.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs

# Remove CONA group, which was a positive control
CONA_ind <- grep("CONA-", rownames(rel_abund))
rel_abund <- rel_abund[-CONA_ind,]

# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% all_metadata[all_metadata$sample_type=="stool" & all_metadata$day >= 0 & all_metadata$cdiff_strain=="431",1 ])]
rel_abund <- rel_abund[overlap,]
end_metadata <- all_metadata[overlap,]
#logCFU <- log10(end_metadata$nextDayCFU_colon+1)

# let's get the taxonomy data so that we have the string from the kingdom to
# the family level name or whatever the next level up is that provided a robust
# classification.
taxonomy_file <- read.table(file="data/process/humanGF.final.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
names(taxonomy) <- rownames(taxonomy_file)
otus <- gsub("tu0*", "TU", names(taxonomy))
names(otus) <- names(taxonomy)

taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub(";unclassified", "", taxonomy)
taxonomy <- gsub("_1", "", taxonomy)
taxonomy <- gsub(";$", "", taxonomy)
taxonomy <- gsub("/.*", "", taxonomy)
taxonomy <- gsub(".*;", "", taxonomy)
taxonomy <- gsub("_sensu_stricto", "", taxonomy)
taxonomy <- gsub("\"$", "", taxonomy) #gets rid of end parenthesis
taxonomy <- gsub("^\"", "", taxonomy) #gets rid of beginning parenthesis

names(taxonomy) <- rownames(taxonomy_file) #the names get deleted for some reason
tax_otu_labels <- paste0("italic('", taxonomy, "')~plain('(", otus, ")')")
names(tax_otu_labels) <- names(taxonomy)
otu_labels <- paste0(otus, "_", taxonomy)
names(otu_labels) <- names(taxonomy)

# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each experimental unit (each cage-day combination)
remove(max) # I accidentally changed mine
cage_med <- aggregate(rel_abund, by=list(end_metadata$cage_id, end_metadata$day), median)
# This function will select a subset of days to include in the model 
select_days <- function(days){
  #day_index <- which(cage_med$Group.2 %in% days)
  cage_med <- cage_med[which(cage_med$Group.2 %in% days),]
  #cageIds <- cage_med[,1]
  cage_med <- cage_med[,-1]
  cage_med <- cage_med[,-1]
  
  good_otus <- apply(cage_med, 2, max) > 1.0
  good_otus["Otu00008"] <- FALSE # Remove the C. difficile OTU
  rel_abund_good <- rel_abund[,good_otus]
  colnames(rel_abund_good) <- otu_labels[colnames(rel_abund_good)] #just for now
  #colnames(rel_abund_good) <- taxonomy[colnames(rel_abund_good)]
  
  # Select days from rel_abund_good
  day_index <- which(as.numeric(gsub(".*-D", "", rownames(rel_abund_good))) %in% days)
  rel_abund_good <- rel_abund_good[day_index,]
  # then properly select logCFU
  end_metadata <- end_metadata[day_index,]
  logCFU <- log10(end_metadata$nextDayCFU_colon+1)
  return(list(rel_abund_good, end_metadata, logCFU))
}

days <- c(0:10)
all_days <- select_days(days)
rel_abund_good <- all_days[[1]]
end_metadata_good <- all_days[[2]]
logCFU <- all_days[[3]]

# see that it actually stabilized by about 2000 trees, but 1e4 sounds more
# impressive :)
n_trees <- 5000

# build the random forest model using log transformed CFU counts.
# Currently this rf_full is over all 10 days
# rf_full.431 is only for cdiff strain 431
set.seed("6201976")
rf_full <- randomForest(logCFU ~ ., data=rel_abund_good, importance=TRUE, ntree=n_trees, na.action=na.omit)

png(height=6, width=10, file=paste0("results/figures/rf_full.png"), unit="in", res=300)
varImpPlot(rf_full, type=1, main=paste0("Cdiff431 CFU predicted over all 10 days, 
                                             rsq: ", round(rf_full$rsq[n_trees], digits=3)))
dev.off()
write(paste("full", ncol(rel_abund_good), rf_full$rsq[n_trees], sep="\t"), file="data/process/random_forest.data")

# sort the features by their effect on % increase in the standard error
importance_sorted <- sort(importance(rf_full)[,1], decreasing=T)

# fit the model back on the data and tack it to the end of the counts file
end_metadata_good$fit_full <- predict(rf_full, rel_abund_good)

# let's plot the top 12 features
n_features <- 12
importance_subset <- importance_sorted[1:n_features]

# Which samples have the top OTUs?
rownames(rel_abund_good[rel_abund_good$OTU110_Coriobacteriaceae > 1, ])
rownames(rel_abund_good[rel_abund_good$OTU42_Lachnospiraceae > 1, ])
rownames(rel_abund_good[rel_abund_good$OTU44_Clostridium > 1, ])
rownames(rel_abund_good[rel_abund_good$OTU17_Clostridium_XlVa > 1, ])
rownames(rel_abund_good[rel_abund_good$OTU5_Akkermansia > 1, ])

# let's get Rsq for the subset
set.seed("6201976")
rf_partial <- randomForest(logCFU ~ ., data=rel_abund_good[,names(importance_sorted)[1:n_features]], importance=TRUE, ntree=n_trees, na.action=na.omit)
png(height=6, width=10, file=paste0("results/figures/rf_partial.png"), unit="in", res=300)
varImpPlot(rf_partial, type=1, main=paste0("Cdiff431 CFU predicted over all 10 days, partial, 
                                             rsq: ", round(rf_partial$rsq[n_trees], digits=3)))
dev.off()
write(paste("partial", n_features, rf_partial$rsq[n_trees], "\t"), file="data/process/random_forest.data", append=T)


# Now let's look at just Day 0 to 1 results
days <- c(0)
day0 <- select_days(days)
rel_abund_good <- day0[[1]]
end_metadata_good <- day0[[2]]
logCFU <- day0[[3]]

set.seed("6201976")
rf_day0 <- randomForest(logCFU ~ ., data=rel_abund_good, importance=TRUE, ntree=n_trees, na.action=na.omit)

png(height=6, width=10, file=paste0("results/figures/rf_day0.png"), unit="in", res=300)
varImpPlot(rf_day0, type=1, main=paste0("Cdiff431 CFU predicted for Day 1, 
                                             rsq: ", round(rf_day0$rsq[n_trees], digits=3)))
dev.off()
write(paste("day0", ncol(rel_abund_good), rf_day0$rsq[n_trees], "\t"), file="data/process/random_forest.data", append=T)

end_metadata_good$fit_full <- predict(rf_day0, rel_abund_good)

# Binary: yes/no colonized
colonized <- rep("yes", length(logCFU))
colonized[logCFU<2] <- "no"
binary_model <- randomForest(x=rel_abund_good, y=as.factor(colonized), importance=TRUE, ntree=n_trees, na.action=na.omit)
png(height=6, width=10, file=paste0("results/figures/binary_modelD0.png"), unit="in", res=300)
varImpPlot(binary_model, type=1, main=paste0("Cdiff431 Colonization predicted for Day 1, 
                                             err: ", round(binary_model$err.rate[n_trees], digits=3)))
dev.off()
write(paste("binary_D1", ncol(rel_abund_good), binary_model$err.rate[n_trees], "\t"), file="data/process/random_forest.data", append=T)

# supplemental figure 4: full feature importance plot
tiff(file="results/figures/figureS4.tiff", width=4, height=5.0, unit="in", res=300)

par(mar=c(3,8,0.5,0.5))
plot(NA, yaxt="n", xlab="", ylab="", xlim=c(min(importance_sorted), 100),
     ylim=c(1, length(importance_sorted)), axes=F)

abline(h=1:length(importance_sorted), lty=3, col="gray")
points(x=rev(importance_sorted), y=1:length(importance_sorted), pch=19, cex=0.8)
axis(1, at=seq(0,100,25), label=c("0", "", "50", "", "100"), cex=0.8)
box()
mtext(side=2, line=7.5, adj=0, at=1:length(importance_sorted),
      text=parse(text=rev(tax_otu_labels[names(importance_sorted)])),
      las=2, cex=0.6,
      col=c(rep("black", length(importance_sorted)-12), rep("red", 12)))
mtext(side=1, text="% Increase in MSE", line=2.0)

dev.off()


# Plot functions
make_gray_plot <- function(cage){
  plot_data <- !grepl(paste0(cage, "-"), rownames(end_metadata_good))
  plot(logCFU[plot_data], end_metadata_good$fit_full[plot_data], xlim=c(0,9),
       ylim=c(0,9), xlab="", ylab="", cex=0.8, axes=F, col="gray", pch= 19)
  box()
}

make_colored_plot <- function(cage){
  plot_data <- grepl(paste0(cage, "-"), rownames(end_metadata_good))
  points(logCFU[plot_data], end_metadata_good$fit_full[plot_data],
         cex=0.8, col=clrs[cage], pch=19) 
  axis(1, labels=rep("", 6), at=seq(0,10,2), tick=T)
  axis(2, labels=rep("", 6), at=seq(0,10,2), tick=T)
}

## GRAPH OVER TIME THE ABS(RESIDUALS), THEN SHOW BY CAGE
## update blog

## MAKE THIS A FUNCTION FOR BOTH INDIVID DAYS AND DAY 0
cageID <- gsub("-.*", "", rownames(end_metadata_good))
cageID <- unique(cageID)
numCages <- length(cageID)
clrs <- rainbow(numCages)
names(clrs) <- cageID

png(file="results/figures/alldays_fit_byCage.png", width=6.0, height=14.0, unit="in", res=300)

design <- matrix(1:21, nrow=7, byrow=T)
design <- cbind(c(rep(22,7)), design)
design <- rbind(design, c(23,23,23, 23))
layout(design, widths=c(0.2,1,1,1), heights=c(1,1,1,1,1,1,1,0.3))
#layout.show(design)

par(mar=c(0.5,1,1.5,0.5))

for(i in 1:numCages){
  make_gray_plot(cageID[i])
  make_colored_plot(cageID[i])
  text(x=-0.25, y=9.75, xpd=TRUE, label=cageID[i], adj=c(0,0), font=2)
  if(i %in% c(1, 4, 7, 10, 13, 16, 19)){
    axis(2, las=2)
  }
  if(i %in% c(19, 20, 21)){
    axis(1)
  } 
}
plot.new()
par(mar=rep(0.1,4))
text(x=0.5,y=0.5, label="Predicted colonization (log CFU)", cex=1.2, srt=90)

plot.new()
par(mar=rep(0.1,4))
text(x=0.5,y=0.5, label="Observed colonization (log CFU)", cex=1.2)
dev.off()


# let's build Figure 6 (w/ color & pch)

correlations <- read.table(file="data/process/correlation_analysis.tsv", header=T, row.names=1)
correlations$sig_corrs <- round(correlations$sig_corrs, 2)

rho <- correlations[names(importance_subset), "sig_corrs"]
rho[is.na(rho)] <- "N.S."

tax_otu_imp_labels <- paste0("bolditalic('", taxonomy[names(importance_subset)],
                             "')~bold('(",
                             otus[names(importance_subset)], "; \u03C1=", rho, ")')")
names(tax_otu_imp_labels) <- names(taxonomy[names(importance_subset)])

tiff(file="results/figures/figure6.tiff", width=6.875, height=7.5, unit="in", res=300)

#want to jitter the relative abundance for those mice that had no Cdiff
#colonization
cd_zeroes <- logCFU == 0
logCFU[cd_zeroes] <- runif(sum(cd_zeroes),0,1)

par(mar=c(0.5,0.5,0.5,0.5))

design <- matrix(1:n_features, nrow=4, byrow=T)
design <- cbind(c(rep(13,4)), design)
design <- rbind(design, c(0,14,14,14))
layout(design, widths=c(0.3,1,1,1), heights=c(1,1,1,1,0.3))

for(i in 1:n_features){
  #get the row and column number for each spot in the layout
  row <- ceiling(i/3)
  column <- ((i-1) %% 3) + 1
  
  #extract the relative abundance data for this OTU
  otu_abund <- rel_abund[,names(importance_subset)[i]]
  
  #want to jitter the number of tumors for those mice that had a zero
  #relative abundance
  ra_zeroes <- otu_abund == 0
  otu_abund[ra_zeroes] <- runif(sum(ra_zeroes),1.0e-2,1.5e-2)
  
  #plot the relative abundance with the number of cdiff for each animal. plot
  #on consistent log scaled x-axis for all OTUs. will throw errors because it
  #can't plot zeroes on a log scale
  plot(otu_abund[grouping == "control-NA-NA"], logCFU[grouping == "control-NA-NA"],
       log="x", pch=19,
       col="black",
       cex=0.8,
       ylab="", xlab="",
       xlim=c(1e-2, 100), ylim=c(0,9),
       yaxt="n", xaxt="n"
  )
  
  points(otu_abund[grouping != "control-NA-NA"], logCFU[grouping != "control-NA-NA"],
         pch=pch[grouping[grouping != "control-NA-NA"]],
         col=clrs[grouping[grouping != "control-NA-NA"]],
         cex=0.8
  )
  
  #create a vertical line to denote the limit of detection
  abline(v=2.2e-2, col="gray")
  
  #create a horizontal line to denote the limit of detection
  abline(h=1.5, col="gray")
  
  #put the OTU label in the upper left corner of the plot
  text(x=0.7e-2, y=8.8, label=parse(text=tax_otu_imp_labels[i]), pos=4, cex=0.9)
  
  #if it's on the bottom row, put a customized axis indicating the % rabund
  if(row == 4){
    axis(1, at=c(1.25e-2, 1e-1,1e0,1e1,1e2),
         label=c("0", "0.1", "1", "10", "100"),
         cex.axis=1.5)
  }
  
  #if it's in the first column turn the axis labels to be horizontal
  if(column == 1){
    axis(2, las=2, cex.axis=1.5)
  }
}

plot.new()
text(x=0.15, y=0.5, label="Observed colonization (log CFU)", cex=1.5, srt=90)

plot.new()
text(x=0.5, y=0.2, label="Relative abundance at Day 0 (%)", cex=1.5)

dev.off()
