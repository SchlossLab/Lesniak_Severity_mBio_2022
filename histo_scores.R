# Mouse histopathology analysis code
#
#
# Dependencies:
#
#
#
# Output: 
library('tidyr')
library('ggplot2')
# load files and trim
meta_file   <- 'data/process/human_CdGF_metadata.txt'
meta_file   <- read.table(meta_file, sep = '\t', header = T, row.names = 'sample_id')
histo_file <- 'data/process/histology_scores_txt.txt'
histo_file <- read.table(histo_file, sep = '\t', header = T)

#make mouse_id column and merge
histo_file <- unite_(histo_file, "mouse_id", c("cage_ID", "mouse_ID"), sep="-", remove=FALSE)
histo_file<- unite_(histo_file, "mouse_id_day", c("mouse_id", "day"), sep="-D", remove=FALSE)
histo_full <- merge(meta_file, histo_file, by.x="row.names", by.y="mouse_id_day")

#plot weightloss and summary score to see if they're correlated
line <- lm(histo_full$percent_weightLoss~ histo_full$summary_score)
plot(histo_full$percent_weightLoss~ histo_full$summary_score, ylab= "% weight loss", xlab="histology summary score")
abline(line)
#Adjusted R-squared:  0.1448

cor(histo_full$percent_weightLoss, histo_full$summary_score)
# output: -0.4015095

#plot mild and severe separately 
severe <- subset(histo_full, outcome == "Severe")
mild <- subset(histo_full, outcome == "Mild")

severe_line <- lm(severe$percent_weightLoss~ severe$summary_score)
plot(severe$percent_weightLoss~ severe$summary_score, ylab= "% weight loss", xlab="histology summary score")
abline(severe_line)
#Adjusted R-squared:  0.009939 

mild_line <- lm(mild$percent_weightLoss~ mild$summary_score)
plot(mild$percent_weightLoss~ mild$summary_score, ylab= "% weight loss", xlab="histology summary score")
abline(mild_line)
#Adjusted R-squared:  0.05345 

#plot of summary score by donor
ggplot(histo_file, aes(human_source, summary_score)) + geom_point(aes(fill=outcome, color=outcome), size=2) +ggtitle("histology summary score")

#plot of independent components by donor
ggplot(histo_file, aes(human_source, edema_tissue)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ggtitle("edema score")
ggplot(histo_file, aes(human_source, inflammation_tissue)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ggtitle("inflammation score")
ggplot(histo_file, aes(human_source, epithelial_damage)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ggtitle("epithelial damage")

ggplot(histo_file, aes(human_source, epithelial_damage)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ggtitle("epithelial damage")

#results plotted by outcome 
sum_out_plot <- ggplot(histo_file, aes(outcome, summary_score)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ylab("Disease severity") + theme_bw() +theme(axis.title.x=element_blank(), legend.position="none")
edema_out_plot <-ggplot(histo_file, aes(outcome, edema_tissue)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ylab("Edema") + theme_bw() + theme(axis.title.x=element_blank(), legend.position="none")
inflam_out_plot <- ggplot(histo_file, aes(outcome, inflammation_tissue)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ylab("Inflammation") + theme_bw() + theme(axis.title.x=element_blank(), legend.position="none")
damage_out_plot <- ggplot(histo_file, aes(outcome, epithelial_damage)) + geom_jitter(aes(fill=outcome, color=outcome), size=2) + ylab("Epithelial damage") + theme_bw() + theme(axis.title.x=element_blank(), legend.position="none")
source('code/multiplot.R')
histo_file <- '~/Documents/Schloss_Lab/Schubert_humanCdGF_XXXX_2016/results/figures/figure_histo.pdf'
histo_multiplot <- multiplot(sum_out_plot, edema_out_plot, inflam_out_plot, damage_out_plot, cols=2)
pdf(file=histo_file, width=11, height=9)
multiplot(sum_out_plot, edema_out_plot, inflam_out_plot, damage_out_plot, cols=2)
dev.off()



plot(histo_file$summary_score ~ histo_file$outcome)

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

#figure code
plot_file <- '~/Documents/Schloss_Lab/Schubert_humanCdGF_XXXX_2016/results/figures/figure_histo.pdf'
pdf(file=plot_file, width=11, height=9)
layout(matrix(c(1,
                2,
                3,
                4), 
              nrow=2, byrow = TRUE))

sum_out_plot

edema_out_plot

inflam_out_plot

damage_out_plot

dev.off()


