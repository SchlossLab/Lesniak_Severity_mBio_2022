# Builds figure 2
# This script builds the mouse weight loss figure
#
#
# Dependencies: 
#     * data/process/human_CdGF_metadata.txt
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#
#
# Output: * results/figures/figure_2.pdf

#load dependencies 
pack_used <- c('tidyr','ggplot2','dplyr', 'knitr')
for (dep in pack_used){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), repos = 'http://cran.us.r-project.org', 
                     quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
# need to change these file paths so they match yours duhhh

meta_file   <- 'data/process/human_CdGF_metadata.txt'
shared_file <- 'data/process/human_CdGF.an.unique_list.0.03.subsample.shared'
tax_file <- 'data/process/gf_new.an.taxonomy'
tax_function <- 'code/tax_level.R'

# read in files
meta_file   <- read.table(meta_file, sep = '\t', header = T, row.names = 'sample_id')
shared_file <- read.table(shared_file, sep = '\t', header = T, row.names = 'sample_id')
tax_file <- read.table(tax_file, sep = '\t', header = T, row.names = 'OTU')
source(tax_function)

# select only mice with 431 cdiff
meta_file <- meta_file[meta_file$cdiff_strain == '431', ]

# subset shared
#Create df with relative abundances
stool_shared <- shared_file[rownames(meta_file[meta_file$day %in% 0, ]) ,]
# use %in% instead of == because to ignore NA
# rename day 0 community rows by mouse id - since meta rows have day in row names, use mouse id for both
rownames(stool_shared) <- meta_file$mouse_id[rownames(meta_file) %in% 
                                               rownames(stool_shared) ]
rel_abund_stool <- 100 * stool_shared / unique(apply(stool_shared, 1, sum))
#Create vector of OTUs with median abundances >1%
OTUs_stool_0_01 <- apply(rel_abund_stool, 2, max) > 1
#df of OTUs with abundances >1% - by cage and inoculum
rel_abund_stool <- rel_abund_stool[ ,OTUs_stool_0_01]
OTU_list <- colnames(rel_abund_stool)

# Weight loss  
wt_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(percent_weightLoss) ) %>%
  filter(!cage_id == 'NP1') %>% 
  select(mouse_id, day, Early_Euth, percent_weightLoss, cage_id) %>% 
  ggplot(aes(x = factor(day), y = percent_weightLoss, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.1) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Persistent", "Severe")) +
  stat_summary(fun.y=median, geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = '% weight change from day 0') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.1), 
        legend.title=element_blank()) +
  stat_summary(fun.data = 'median_hilow')

# C. Diff CFU
cfu_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  filter(!cage_id == 'NP1') %>% 
  select(mouse_id, day, Early_Euth, log_cfu, cage_id) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.2) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Persistent", "Severe")) +
  stat_summary(fun.y=median, geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = ' CFU (log10)') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.2), 
        legend.title=element_blank()) +
  stat_summary(fun.data = 'median_hilow')

# weight loss plots colored by histo score group/bins

#first need to make a table with summary scores in meta
bintest <- merge(meta_file, three_bins, by.x="mouse_id", by.y="mouse_id.x")

bins_wt_plot <- bintest %>% 
  filter(day %in% c(0:10), !is.na(percent_weightLoss) ) %>%
  filter(!cage_id == 'NP1') %>% 
  select(mouse_id, day, Early_Euth, percent_weightLoss, cage_id, severity) %>% 
  ggplot(aes(x = factor(day), y = percent_weightLoss, color = factor(severity))) + 
  geom_jitter(width = 0.5, alpha = 0.5)  +
  stat_summary(fun.y=median, geom="line", aes(group = severity)) +
  labs(x = 'Day', y = '% weight change from day 0') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.1), 
        legend.title=element_blank()) +
  stat_summary(fun.data = 'median_hilow') + theme_bw()

#just severe group, weight loss by mouse and cage 
severebin <- bintest %>% 
  filter(day %in% c(0:10), !is.na(percent_weightLoss) ) %>%
  filter(!cage_id == 'NP1') %>% 
  filter(severity == 'severe') %>%
  select(mouse_id, day, Early_Euth, percent_weightLoss, cage_id) %>%
  ggplot(aes(x=factor(day), y = percent_weightLoss, color=mouse_id)) +geom_point()

# just severe group, weight loss by cage
ggplot(severebin, aes(factor(day), percent_weightLoss))+ geom_point(aes(color=cage_id))

#cfu bins
cfu_plot <- bintest %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  filter(!cage_id == 'NP1') %>% 
  select(mouse_id, day, Early_Euth, log_cfu, cage_id, severity) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(severity))) + 
  geom_jitter(width = 0.5, alpha = 0.2) +
  stat_summary(fun.y=median, geom="line", aes(group = severity)) +
  labs(x = 'Day', y = ' CFU (log10)') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.2), 
        legend.title=element_blank()) +
  stat_summary(fun.data = 'median_hilow')





#build and save plot 

plot_file <- '~/Documents/Schloss_Lab/Schubert_humanCdGF_XXXX_2016/results/figures/figure_2.pdf'
pdf(file=plot_file, width=11, height=9)
layout(matrix(c(1,
                2), 
              nrow=2, byrow = TRUE))

wt_plot

cfu_plot 

dev.off()

