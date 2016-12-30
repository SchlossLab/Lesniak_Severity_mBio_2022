# Mouse histopathology analysis code
#
#
# Dependencies:
#
#
#
# Output: 

# load files and trim
meta_file   <- 'data/process/human_CdGF_metadata.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, row.names = 'sample_id')
histo_file <- 'data/process/histology_scores_txt.txt'
histo_file <- read.table(histo_file, sep = '\t', header = T)

#plot of summary score by donor
ggplot(histo_file, aes(human_source, summary_score)) + geom_point(aes(fill=outcome, color=outcome), size=2) +ggtitle("histology summary score")

#plot of independent components by donor
ggplot(histo_file, aes(human_source, edema_tissue)) + geom_point(aes(fill=outcome, color=outcome), size=2) + ggtitle("edema score")
ggplot(histo_file, aes(human_source, inflammation_tissue)) + geom_point(aes(fill=outcome, color=outcome), size=2) + ggtitle("inflammation score")
ggplot(histo_file, aes(human_source, epithelial_damage)) + geom_point(aes(fill=outcome, color=outcome), size=2) + ggtitle("epithelial damage")

#wilcoxon test for summary scores and edema 
#store scores for mild and sevvere in separate variables
mild_summary <- histo_file[histo_file$outcome == 'Mild',]
severe_summary <- histo_file[histo_file$outcome == 'Severe', ]
wilcox.test(mild_summary$summary_score, severe_summary$summary_score)
# W = 15.5, p-value = 5.681e-09

wilcox.test(mild_summary$edema_tissue, severe_summary$edema_tissue)
# W = 12.5, p-value = 2.456e-09

wilcox.test(mild_summary$inflammation_tissue, severe_summary$inflammation_tissue)
# W = 114, p-value = 7.755e-06

wilcox.test(mild_summary$epithelial_damage, severe_summary$epithelial_damage)
#W = 61, p-value = 2.5e-07



