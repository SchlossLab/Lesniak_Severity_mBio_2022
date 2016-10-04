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

meta_file <- meta_file[meta_file$cdiff_strain == '431', ]
meta_file <- meta_file[meta_file$cage_id != 'NP1',]


wt_ls_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(percent_weightLoss) ) %>%
  select(mouse_id, day, Early_Euth, percent_weightLoss, cage_id) %>% 
  ggplot(aes(x = factor(day), y = percent_weightLoss, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.1) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Mild", "Severe")) +
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
                     labels=c("Mild", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = ' CFU (log10)', title = 'C. Difficile Colonization') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.2), 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

cfu_cage_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  select(mouse_id, day, Early_Euth, log_cfu, cage_id) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(cage_id))) + 
  geom_jitter(width = 0.5, alpha = 0.2) + 
  #scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     #labels=c("Mild", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = as.numeric(cage_id))) +
  labs(x = 'Day', y = ' CFU (log10)', title = 'C. Difficile Colonization') +
  theme(legend.position="right", 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

cfu_donor_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  select(mouse_id, day, Early_Euth, log_cfu, human_source) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(human_source))) + 
  geom_jitter(width = 0.5, alpha = 0.2) + 
  #scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
  #labels=c("Mild", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = as.numeric(human_source))) +
  labs(x = 'Day', y = ' CFU (log10)', title = 'C. Difficile Colonization') +
  theme(legend.position="right", 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

cfu_mouse_plot <- meta_file %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  select(mouse_id, day, Early_Euth, log_cfu, cage_id, ear_tag) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(cage_id))) + 
  geom_jitter(width = 0.5, alpha = 0.2) + 
  #scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
  #labels=c("Mild", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = as.numeric(ear_tag))) +
  labs(x = 'Day', y = ' CFU (log10)', title = 'C. Difficile Colonization') +
  theme(legend.position="right", 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

#plot SD for median, median for weight loss, cut colonization off at day 5


ggdraw() +
  draw_plot(wt_ls_plot, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(cfu_plot, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(cfu_cage_plot, x = 0, y=0, width=0.5, height=0.5) +
  draw_plot(cfu_donor_plot, x = 0.5, y=0, width=0.5, height=0.5)

ggdraw() +
  draw_plot(cfu_mouse_plot)


#5 day plots only
d5cfu_plot <- meta_file %>% 
  filter(day %in% c(0:5), !is.na(log_cfu) ) %>%
  select(mouse_id, day, Early_Euth, log_cfu, cage_id) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.2) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Mild", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = ' CFU (log10)', title = 'C. Difficile Colonization') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.2), 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

d5wt_ls_plot <- meta_file %>% 
  filter(day %in% c(0:5), !is.na(percent_weightLoss) ) %>%
  select(mouse_id, day, Early_Euth, percent_weightLoss, cage_id) %>% 
  ggplot(aes(x = factor(day), y = percent_weightLoss, color = factor(Early_Euth))) + 
  geom_jitter(width = 0.5, alpha = 0.1) + 
  scale_color_manual(values = c('black', 'red'), breaks=c("FALSE", "TRUE"),
                     labels=c("Mild", "Severe")) +
  stat_summary(fun.y="median", geom="line", aes(group = Early_Euth)) +
  labs(x = 'Day', y = ' Percent Relative to Day 0', title = 'Weight Loss') +
  theme(legend.justification=c(1,0), legend.position=c(1,0.1), 
        legend.title=element_blank()) +
  stat_summary(fun.y= "median", geom = "point")

ggdraw()+
  draw_plot(d5cfu_plot, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(d5wt_ls_plot, x = 0, y = 0, width = 0.5, height = 0.5)

#half working solution for median error bars
#stat_summary(fun.data = median_hilow, conf.int=.5)



####### Not 431 strains #######

#nick's code

non_431   <- 'data/process/human_CdGF_metadata.txt'

# read in files
non_431   <- read.table(non_431, sep = '\t', header = T, row.names = 'sample_id')

# select mice infected with two strains of c diff
cdiff_strains <- unique(non_431$cdiff_strain)
non_431_strains <- cdiff_strains[!unique(cdiff_strains) %in% c('431', NA)]
non_431_inocula <- unique(non_431[non_431$cdiff_strain %in% non_431_strains, 'human_source'])
non_431 <- non_431[non_431$human_source %in% non_431_inocula, ]

non_431_cfu <- non_431 %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  select(mouse_id, day, human_source, cdiff_strain, cage_id, ear_tag, log_cfu) %>% 
  ggplot(aes(x = factor(day), y = log_cfu, color = factor(cdiff_strain))) + geom_point() +
  stat_summary(fun.y=median, geom="line", aes(group = cage_id)) + 
  facet_grid(. ~ human_source) + theme_bw() + 
  labs(x = 'Day', y = 'C Diff CFU (log10)') +
  theme(legend.justification=c(1,0.5), legend.position=c(1,0.5), 
        legend.title=element_blank())


non_431_wt <- non_431 %>% 
  filter(day %in% c(0:10), !is.na(log_cfu) ) %>%
  select(mouse_id, day, human_source, Early_Euth, percent_weightLoss, cdiff_strain, ear_tag, cage_id, log_cfu) %>% 
  ggplot(aes(x = factor(day), y = percent_weightLoss, color = factor(cdiff_strain))) + geom_point() +
  stat_summary(fun.y=median, geom="line", aes(group = cage_id)) + 
  facet_grid(. ~ human_source) + theme_bw() + 
  labs(x = 'Day', y = 'Percent weight loss relative to day 0') +
  theme(legend.justification=c(1,0.5), legend.position=c(1,0.5), 
        legend.title=element_blank())

ggdraw() +
  draw_plot(non_431_wt,  x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(non_431_cfu,x = 0, y = 0.5, width = 0.5, height = 0.5)
