# rename_groups.R
# since there is - in the sample nanmes and this is the symbol mothur uses to split sames
# mothur attempts to split all the file names containing - 
# replace all - with another character

groups_file <- "data/mothur/humanGF_cdiff.contigs.good.groups"

groups <- read.table(groups_file)
groups$V2 <- gsub("[-_]", 'v', groups$V2)
file.rename(from = groups_file, to = paste0(groups_file, '.original'))
write.table(groups, groups_file, sep = '\t',
	qoute = F, row.names = F, col.names = F)