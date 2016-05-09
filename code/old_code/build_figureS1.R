################################################################################
#
# build_figure1.R
#
# This script builds Figure 1, which is a barchart of the median relative
# abundance of each phylum found in the mice treated with the top dose of
# antibiotics as well as the untreated control mice. Also included is the
# number of CFU per gram of feces.
#
# Dependencies...
#   * data/process/humanGF.final.tx.5.subsample.shared
#   * data/process/humanGF.final.tx.5.cons.taxonomy
#   * data/process/abxD1.counts
#
# Output...
#   * results/figures/figure1.tiff
#
################################################################################



#from last file... old stuff
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
top_dose <- counts_file[counts_file$experiment=="top_dose" | counts_file$abx=="control",]
top_dose_med <- aggregate(top_dose$CFU, by=list(top_dose$abx), median)
abx_cfu <- format(top_dose_med$x, scientific=T, digits=2)
abx_cfu <- gsub("e\\+0", "E", abx_cfu)
abx_cfu <- gsub("0.0E0", "0", abx_cfu)
names(abx_cfu) <- top_dose_med$Group.1


# read in the metadata file
all_metadata <- read.table(file="data/raw/humanGF_metadata.txt", header=T, colClasses="character")
D0 <- all_metadata[all_metadata$day=="0" & all_metadata$cdiff_strain=="431", ]
D0 <- D0[!is.na(D0$group),]
rownames(D0) <- D0$group

# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/humanGF.final.tx.5.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% all_metadata[,1])]
rel_abund <- rel_abund[overlap,]
D0 <- D0[overlap,]


# let's get the relative abundances for those phyla that have at least one
# sample where they are more than 10% o the community
good_phyla <- apply(rel_abund, 2, max) > 10
rel_abund_good <- rel_abund[,good_phyla]


# let's get the taxonomy data
taxonomy_file <- read.table(file="data/process/humanGF.final.tx.5.cons.taxonomy", header=T, row.names=1)
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

#let's get the medians and IQRs
med_ra <- aggregate(rel_abund_good, by=list(D0$cage_id), median)
cage <- med_ra[,1]
med_ra <- med_ra[,-1]
rownames(med_ra) <- cage
med_ra <- t(med_ra)

l_qtr <- aggregate(rel_abund_good, by=list(D0$cage_id), function(x){quantile(x, probs=c(0.25))})
l_qtr <- l_qtr[,-1]
rownames(l_qtr) <- cage
l_qtr <- t(l_qtr)

u_qtr <- aggregate(rel_abund_good, by=list(D0$cage_id), function(x){quantile(x, probs=c(0.75))})
u_qtr <- u_qtr[,-1]
rownames(u_qtr) <- cage
u_qtr <- t(u_qtr)

D0_phylum_barplot <- function(cage, label){
    pos <- barplot(as.vector(med_ra[,cage]), ylim=c(0,105), axes=F, col="white")
    arrows(x0=pos, x1=pos, y0=med_ra[,cage], y1=u_qtr[,cage], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,cage], y1=l_qtr[,cage], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    box()
    text(x=par("usr")[1], y=par("usr")[4]*1.05, label=label, font=2,
                                                    adj = c(0,0), xpd=TRUE)
    pos
}

#need this to be 2 by 11, loop through each cage
#then another graph add the cdiff CFU

tiff(height=9, width=3.75, file="results/figures/figureS1.tiff", unit="in", res=300)

    z <- layout(
        matrix( c(  seq(1, 24)), byrow=F, ncol=2 ), 
        widths=c(1), heights=c(rep(1,8), 1)
        ) #there are 22 different cages

    par(mar=c(0.5, 5, 1.25, 0.5))
    numCages <- length(colnames(med_ra))
    for (i in 1:numCages){ 
      print(i)
      D0_phylum_barplot(colnames(med_ra)[i], colnames(med_ra)[i])
    }

   # mtext(side=2, line=3, at=0, "Relative Abundance (%)")

    text(x=pos+0.1, y=par("usr")[3]-10, labels=rownames(med_ra), srt=70, cex=1, font=3, pos=2, xpd=NA)

dev.off()
