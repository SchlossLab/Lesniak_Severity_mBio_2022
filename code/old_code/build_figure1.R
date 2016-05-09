################################################################################
#
# build_figure1.R
#
# This script builds Figure 1, which is a barchart of the median relative
# abundance of each genus found in the mice D0
#
# Dependencies...
#   * data/process/humanGF.final.an.unique_list.0.03.subsample.shared
#   * data/process/humanGF.final.an.unique_list.0.03.cons.taxonomy
#   * data/raw/humanGF_metadata.txt
#
# Output...
#   * results/figures/figure1.tiff
#
################################################################################


# let's get the taxonomy data so that we have the string from the kingdom to
# the genus level name or whatever the next level up is that provided a robust
# classification.
taxonomy_file <- read.table(file="data/process/humanGF.final.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
names(taxonomy) <- rownames(taxonomy_file)
taxonomy <- gsub("\\(\\d*\\)", "", taxonomy) #gets rid of bootstrap values
taxonomy <- gsub(";unclassified", "", taxonomy) #gets rid of all unclassified
taxonomy <- gsub("/.*", "", taxonomy) 
taxonomy <- gsub(";$", "", taxonomy) #gets rid of extra ; at the end
taxonomy <- gsub("\"$", "", taxonomy) #gets rid of end parenthesis
taxonomy <- gsub(".*;", "", taxonomy) #eliminates all preceding text, except last classification
taxonomy <- gsub("^\"", "", taxonomy) #gets rid of beginning parenthesis
#taxonomy <- taxonomy[sig_otus]




# read in the metadata file
all_metadata <- read.table(file="data/raw/humanGF_metadata.txt", header=T, colClasses="character")
#top_dose <- counts_file[counts_file$experiment=="top_dose" | counts_file$abx=="control",]
cd431 <- all_metadata[all_metadata$cdiff_strain=="431",]

# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/humanGF.final.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% all_metadata[,1])]
rel_abund <- rel_abund[overlap,]
cd431 <- cd431[overlap,]

############################################################################
# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each antibiotic dose
control_rabund <- rel_abund[top_dose$abx == "control",]
control_metadata <- top_dose[top_dose$abx == "control",]
control_median <- apply(control_rabund, 2, median)
control_otus <- control_median > 1.0

otu_hyp_test <- function(drug){

    drug_rabund <- rel_abund[top_dose$abx == drug,]
    drug_metadata <- top_dose[top_dose$abx == drug,]
    drug_otus <- apply(drug_rabund, 2, median) > 3.0

    combined_otus <- control_otus | drug_otus
    combined_rabund <- rbind(drug_rabund, control_rabund)[combined_otus]

    drugged <- c(rep(TRUE, nrow(drug_rabund)), rep(FALSE, nrow(control_rabund)))
    n_otus <- ncol(combined_rabund)

    p_values <- rep(NA, n_otus)

    warn_orig <- options("warn")$warn
    options("warn"= -1)

    for(i in 1:n_otus){
        p_values[i] <- wilcox.test(combined_rabund[,i], g=drugged)$p.value
    }

    options("warn" = warn_orig)

    adj_p_values <- p.adjust(p_values, method="BH")
    colnames(combined_rabund)[adj_p_values<0.05]
}

amp_sig_otus <- otu_hyp_test("amp")
cef_sig_otus <- otu_hyp_test("cef")
cipro_sig_otus <- otu_hyp_test("cipro")
clinda_sig_otus <- otu_hyp_test("clinda")
metro_sig_otus <- otu_hyp_test("metro")
strep_sig_otus <- otu_hyp_test("strep")
vanc_sig_otus <- otu_hyp_test("vanc")

sig_otus <- sort(unique(c(amp_sig_otus, cef_sig_otus, cipro_sig_otus,
                            clinda_sig_otus, metro_sig_otus, strep_sig_otus,
                            vanc_sig_otus)))
x_max <- length(sig_otus) * 1.2

o <- order(control_median[sig_otus], decreasing=T)
rel_abund_sig <- rel_abund[,sig_otus[o]]


otu <- gsub("Otu0*", "OTU~", names(taxonomy))
names(otu) <- names(taxonomy)

tax_label <- paste0("italic('", taxonomy[sig_otus[o]], "')~(", otu[sig_otus[o]], ")")


single_drug_bars <- function(drug, drug_sig_otus, drug_label){

#    drug <- "control"
#    drug_sig_otus <- ""
#    drug_label <- "No antibiotics"
    drug_rabund <- rel_abund_sig[top_dose$abx == drug,]
    drug_metadata <- top_dose[top_dose$abx == drug,]

    n <- nrow(drug_rabund)

    drug_med <- apply(drug_rabund, 2, median)
    drug_uci <- apply(drug_rabund, 2, function(x){quantile(x, prob=0.75)})
    drug_lci <- apply(drug_rabund, 2, function(x){quantile(x, prob=0.25)})

    z <- barplot(drug_med, names.arg=rep("", length(drug_med)),
                    ylim=c(0,1+max(drug_uci)), xlim=c(0,x_max), axes=F,
                    col="white")

    warn_orig <- options("warn")$warn
    options("warn"= -1)

    arrows(x0=z, y0=drug_med, y1=drug_uci, angle=90, length=0.05)
    arrows(x0=z, y0=drug_med, y1=drug_lci, angle=90, length=0.05)

    options("warn"= warn_orig)

    text(x=z[sig_otus[o] %in% drug_sig_otus], y=-0.05*max(drug_uci),
                                            labels="*", cex=2, xpd=TRUE)

    axis(2, las=1)
    box()

    if(grepl(")", drug_label)){
        drug_label <- gsub(")", paste0("; N=", n, ")"), drug_label)
    } else {
        drug_label <- paste0(drug_label, "; N=", n, ")")
    }

    text(x=par("usr")[1], y=par("usr")[4]*1.175, label=drug_label,
                                adj=c(0,1), cex=1.2, font=2, xpd=TRUE)




#    summary_stats <- format(quantile(drug_metadata$CFU,
#                            prob=c(0.25, 0.50, 0.75)), scientific=T, digits=2)
#
#    summary_string <- paste0(summary_stats[2], "~(", summary_stats[1], "-", summary_stats[3], ")")
#
#    summary_string <- gsub("(\\d\\.\\d*)e\\+0", "plain('\\1x10')^", summary_string)
#    summary_string <- gsub("0e\\+00", "plain('<1x10')^2", summary_string)
#
#    summary_string <- paste0(summary_string, "~plain(' N=", n, "')")
#
#    text(x=par("usr")[2], y=1.05*par("usr")[4], labels=parse(text=summary_string),
#                                adj=c(1,0), pos=2, cex=0.8, xpd=TRUE)

    summary_stats <- quantile(drug_metadata$CFU, prob=c(0.25, 0.50, 0.75))
    z2 <- barplot(summary_stats[2]+1, width=0.3, xlim=c(-0.1,0.5), ylim=c(1,1e9), log="y", axes=FALSE, col="white", names.arg="")

    warn_orig <- options("warn")$warn
    options("warn"= -1)

    arrows(x0=z2, y0=summary_stats[2]+1, y1=summary_stats[3]+1, angle=90, length=0.05)
    arrows(x0=z2, y0=summary_stats[2]+1, y1=summary_stats[1]+1, angle=90, length=0.05)

    options("warn"= warn_orig)

    box()
    axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))

    if(summary_stats[2] == 0){
        text(x=0.2, y=10, label=expression(plain('<10')^2))
    }

    z

}



tiff(file="results/figures/figure1.tiff", width=4.5, height=10.0, unit="in", res=300)
    par(cex=1.2)

    layout_matrix <-matrix(c(
            1,  2,
            3,  4,
            5,  6,
            7,  8,
            9, 10,
           11, 12,
           13, 14,
           15, 16,
           17, 18
        ),nrow=9, byrow=T)
    layout_matrix <- cbind(c(rep(19,8),0), layout_matrix, c(rep(20,8),0))

    layout(layout_matrix, width=c(0.2, 1.25, 0.15, 0.2), height=c(rep(1,8),1.5))

#    par(mar=c(0.5,5,1.5,0.5))
    par(mar=c(0.75,0.25,1.5,0.5), oma=c(0,0,0,0))

    z <- single_drug_bars("control", "", "No antibiotics")
    z <- single_drug_bars("cipro", cipro_sig_otus, "Ciprofloxacin (10 mg/kg)")
    z <- single_drug_bars("vanc", vanc_sig_otus, "Vancomycin (0.625 mg/mL)")

    z <- single_drug_bars("amp", amp_sig_otus, "Ampicillin (0.5 mg/mL)")
    z <- single_drug_bars("clinda", clinda_sig_otus, "Clindamycin (10 mg/kg)")
    z <- single_drug_bars("strep", strep_sig_otus, "Streptomycin (5 mg/mL)")
    z <- single_drug_bars("metro", metro_sig_otus, "Metronidazole (1 mg/mL)")
    z <- single_drug_bars("cef", cef_sig_otus, "Cefoperazone (0.5 mg/mL)")

    plot(NA, type="n", axes=F, xlim=c(0, max(z)), ylim=c(0,1))
    text(x=z+0.5, y=1.2, xpd=NA, label=parse(text=tax_label), pos=2, srt=70, cex=1.0)


        plot.new()

    plot.new()
    par(mar=c(0.1,0.1,0.1,0.1))
    text(x=0.25, y=0.5, "Relative abundance (%)", srt=90, cex=1.5)

    plot.new()
    par(mar=c(0.1,0.1,0.1,0.1))
    text(x=0.75, y=0.5, expression(italic('C. difficile')~colonization~(CFU/g)), srt=-90, cex=1.5)

dev.off()
