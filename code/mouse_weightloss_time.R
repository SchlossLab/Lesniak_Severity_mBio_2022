##########
#
#
# This script builds the mouse weight loss figure
#
#
# Dependencies: 
#     * data/process/human_CdGF_metadata.txt
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#     * data/mothur/gf_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#
#
# Output: * results/figures/weight_loss_time.pdf

install.packages("cowplot")
install.packages("dplyr")
library(cowplot)
library(dplyr)

meta_file <- read.table("data/process/human_CdGF_metadata.txt", header = TRUE, sep='\t', fill = TRUE, row.names=2)
shared_file <- read.table("data/process/human_CdGF.an.unique_list.0.03.subsample.shared", sep='\t', header = TRUE, row.names=1)
tax_file <- read.table('data/process/gf_new.an.taxonomy', sep='\t',header = T, row.names = 1)


#########All below is Nick's code that I am modifying ######

wt_ls_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(percent_weightLoss) ) %>%
  select(mouse_id, day, Early_Euth, percent_weightLoss, cage_id) %>% 
  ggplot(aes(x = factor(day), y = percent_weightLoss, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.1) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Persistent", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = ' Percent Relative to Day 0', title = 'Weight Loss') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.1), 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")


cfu_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  select(mouse_id, day, Early_Euth, log_cfu, cage_id) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.2) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Persistent", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = ' CFU (log10)', title = 'C. Difficile Colonization') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.2), 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

ggdraw() +
  draw_plot(wt_ls_plot, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(cfu_plot, x = 0, y = 0.5, width = 0.5, height = 0.5)

#half working solution for median error bars
stat_summary(fun.data = median_hilow, conf.int=.5)
