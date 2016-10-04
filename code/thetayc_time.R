###################
#
# thetayc_time.R
#
#Plot thetaYC over time distance from day 0 
#thetaYC distances cant have sample names that are just numbers, this script wont work 
#    
#
###################
install.packages("reshape")
library("reshape")
install.packages("gridExtra")
library("gridExtra")
library("ggplot2")
source('code/read.dist.R')

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/human_CdGF.an.unique_list.0.03.subsample.shared', sep='\t',header=T)
theta <- read.dist(file='data/process/shared_subset_2194.thetayc.0.03.lt.dist', input="lt", make.square=T, diag=0)



#then maybe write it as a function 

#################THE RIGHT ONE!!!#######
#compares day0 to every subsequent day. this asks, how different are the communities from the start? can this be linked to outcome?

mouseIDVector <- as.character(unique(meta_file$mouse_id))
mouseIDVector <- mouseIDVector[mouseIDVector != "inoculum_NA"]
stability <- as.data.frame(matrix(ncol=10, nrow=length(mouseIDVector)))
row.names(stability) <- mouseIDVector
names(stability) <- c(1,2,3,4,5,6,7,8,9,10)
for (k in 1:10){
  sampA <- rownames(meta_file)[meta_file$day %in% 0]
  sampB <- rownames(meta_file)[meta_file$day %in% k]
  dist <- theta[sampA, sampB]
  row.names(dist) <- meta_file$mouse_id[rownames(meta_file) %in% rownames(dist)]
  names(dist) <- meta_file$mouse_id[rownames(meta_file) %in% names(dist)]
  for (i in mouseIDVector){
    coord <- dist[i,i]
    try(
    stability[i,k] <- as.numeric(coord), silent = TRUE
    )
  }
}

#need to melt data to get it in plottable form. add mouse id and cage id column. WOO. 


melted <- melt(stability) %>% mutate(mouse_id = rep(rownames(stability), 10))
melted$cage_id <- sapply((strsplit(melted$mouse_id, "_")), '[', 1)


#ggplot wont plot NAs but if you remove them it plots the lines ok 
ggplot(na.omit(melted), aes(variable, value, group = mouse_id)) + geom_line(aes(color = cage_id)) + theme_bw() + facet_grid(cage_id~.) + ggtitle("ThetaYC vs day 0 over time") + xlab("day post infection") + ylab("Theta YC distance from day 0")


#not faceted or averaged
ggplot(na.omit(melted), aes(variable, value, group = mouse_id)) + geom_line(aes(color = cage_id)) + theme_bw() + ggtitle("ThetaYC vs day 0 over time") + xlab("day post infection") + ylab("Theta YC distance from day 0")

ggplot(na.omit(melted), aes(variable, value, group = mouse_id)) + geom_line(aes(color = mouse_id)) + theme_bw()

#trying multiple plots - put this in a loop and have plot all on one page? or use this to avg by cage 

ggplot(na.omit(subset(melted, melted$cage_id == '369')), aes(variable, value, group = mouse_id)) + geom_line() + ggtitle("Cage 369") + xlab("day post infection") +ylab("Theta YC distance from day 0") +theme_bw()

ggplot(na.omit(subset(melted, melted$cage_id == '430')), aes(variable, value, group = mouse_id)) + geom_line() + ggtitle("Cage 430") + xlab("day post infection") +ylab("Theta YC distance from day 0") + theme_bw()

ggplot(na.omit(subset(melted, melted$cage_id == 'OP')), aes(variable, value, group = mouse_id)) + geom_line(aes(color='red')) + ggtitle("Cage OP") + xlab("day post infection") +ylab("Theta YC distance from day 0") + theme_bw()

p <- list()
i <- 1
for(j in melted$cage_id){
  p[[i]] <- ggplot(na.omit(subset(melted, melted$cage_id == j)), aes(variable, value, group = mouse_id)) + geom_line(aes(color = mouse_id))
  i <- i +1
}
do.call(grid.arrange, p)

#now average by cage. woof. and add labels and shit 
#or at least split up plot by severity 
#facet by cage!!!

p <- list()
for(j in melted$cage_id){
  for(i in 1:length(unique(melted$cage_id))){
    p[[i]] <- ggplot(na.omit(subset(melted, melted$cage_id == j)), aes(variable, value, group = mouse_id)) + geom_line(aes(color = mouse_id))
  }
}
do.call(grid.arrange, p)
#build as function 


#now write version where you compare each day to previous day. this asks question: how stable is community over time?



