################################################################################
#
# make_metadata_file.R
#
# This script will build the metadata file that will indicate the density of C.
# difficile on the day after challenge and a break down of the groups. 
#
# Dependencies...
# * data/raw/toxin_data.txt
# * data/raw/humanGF_ids.txt
#
# Output...
# * data/process/humanGF_metadata.txt
#
################################################################################
setwd("~/Documents/Github/Schubert_humanCdGF_2015")


# Read in data
all_metadata <- read.table(file="data/raw/humanGF_ids.txt", header=T)
all_metadata <- all_metadata[, -which(colnames(all_metadata)=="file")]

# Merge toxin data with the other metadata
toxin_data <- read.table(file="data/raw/toxin_data.txt", header=T)
# Go through each sample and fill in toxin data where have it
for(i in 1:dim(all_metadata)[1]){
  if(all_metadata$group[i] %in% toxin_data$Cage_Mouse){
    all_metadata[i, "toxin_log"] <- toxin_data$Log_repiricoal_dilution[which(toxin_data$Cage_Mouse%in% all_metadata$group[i])]
  } 
}

#output metadata file
write.table(all_metadata, file="data/process/humanGF_metadata.txt", quote=FALSE, sep="\t")
