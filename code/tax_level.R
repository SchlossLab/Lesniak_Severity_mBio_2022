#get taxonomy - df with columns for each level - 
# level - 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum), 6 (kingdom)
## - convert taxonomy text list to df, remove percentages
## - subset df by desired tax level, then replace any unclassified with next level up
get_tax <- function(tax_level=1,row_list=1:length(rownames(tax_file)),df=tax_file){
     if (tax_level %in% c(1:5)){
          row_list <- as.character(row_list)
          taxonomy <- df[row_list,]
          taxonomy <- strsplit(as.character(taxonomy$Taxonomy),';',fixed=TRUE)
          n <- max(sapply(taxonomy, length))
          taxonomy <- lapply(taxonomy, `[`, seq_len(n))
          taxonomy <- data.frame(do.call('rbind', taxonomy))
          taxonomy <- sapply(taxonomy,gsub,pattern="\\(\\d*\\)",replacement="")
          level <- 7-tax_level
          tax_out <- as.character(taxonomy[,level])
          for (i in level:2){
               next_level <- i-1
               tax_out[is.na(tax_out)] <- 
                 as.character(taxonomy[is.na(tax_out),next_level])
          }
          tax_out <- gsub('_unclassified', '', tax_out)
          label_out <- paste0(tax_out, ' (', 
                              gsub('tu0*', 'TU ', row_list),')')
          return(data.frame(tax= tax_out, tax_label = label_out, row.names=row_list))
     } else {print(
          'Error: Level must be 1 (genus), 2 (family), 3 (order), 4 (class), 5 (phylum)')
     }
}