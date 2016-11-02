###################
#
# InvSimp_d0.R
#
#Plot Inverse Simpson diversity metric for day 0 communities +/or over time
#
#    
#
###################

#read in files
meta_file <- read.table('data/process/human_CdGF_metadata.txt', sep='\t',header = T, row.names = 2)
shared <- read.table('data/process/humanCdGF_full.shared', sep='\t',header=T)
simp <- read.table('data/process/invsimp.summary', sep='\t',header=T)

#subset simp to be day 0 only
day0list <- rownames(meta_file)[meta_file$day==0]
simp[2] <- gsub("-", "", simp$group)
day0list <- gsub("-", "", day0list)

#for the below to work, need to do some sort of gsubbing to get the dashes fixed 
day0_simp <- simp[simp$group %in% day0list,]
day0_meta <- subset(meta_file, day==0)

day0_simp$cage <- day0_meta$cage_id
day0_simp$donor <- day0_meta$human_source
day0donor <- day0_meta$human_source

fullsimp <- cbind(day0_simp[2:3], day0_simp[6:7])

#plot(fullsimp$cage, fullsimp$invsimpson, type="p", cex.lab=0.7, cex.axis=0.7)

stripchart(fullsimp$invsimpson ~ fullsimp$cage, vertical = TRUE, cex.axis=0.7, pch=16, main="Simpson diversity on day 0", xlab="cage", ylab="Inverse Simpson")

stripchart(fullsimp$invsimpson ~ fullsimp$donor, vertical = TRUE, method = "jitter", cex.axis=0.7, pch=16, main="Simpson diversity on day 0", xlab="donor", ylab="Inverse Simpson")

#would be good to label the points by severity.. dont have a column for this but can 
#force it by coloring donors that are severe 
severe <- fullsimp[grep('DA00884', fullsimp$donor), c(4,2)]
temp <- fullsimp[grep('DA00431', fullsimp$donor), c(4,2)]
severe <- rbind(severe, temp)
temp <- fullsimp[grep('DA01245', fullsimp$donor), c(4,2)]
severe <- rbind(severe, temp)
temp <- fullsimp[grep('DA10034', fullsimp$donor), c(4,2)]
severe <- rbind(severe, temp)
temp <- fullsimp[grep('DA10082', fullsimp$donor), c(4,2)]
severe <- rbind(severe, temp)

points(severe, pch=16, col = "red")
legend <- c("severe", "asymptomatic")
legend(x="topright", legend, col = c("red", "black"), pch=16)




#why the fuck are the labels not working 
text(cex=0.8, x=0+0.5, y=0, label=fullsimp$donor, xpd=TRUE, srt=45, pos=2)
text(6, 5, xpd=NA, label=parse(text=day0donor), srt=70, cex=0.8)

#x=0.5, y=-4,
#try making this but in a version of donors
#remove inoculum 



#then can use same principle for RA charts!! 


#plot day0 inv simpson with cage along bottom inv simpson on top
#need to merge w metadata file to get column of cages 

#but would be better to show over time like i did with the cdiff time paper plot 0-10 days, by donor? cage? but dont want to do average, too small n 
