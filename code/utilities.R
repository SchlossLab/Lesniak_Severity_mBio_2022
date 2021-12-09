###############################################################################
#
#	Common functions used across scripts
#
# need files:
# 	produced by code/tidy_taxonomy.R
#		data/process/final.taxonomy.tidy.tsv
#	produced through running mothur scripts
#		data/mothur/sample.final.0.03.subsample.shared
#		data/mothur/sample.final.thetayc.0.03.lt.ave.dist
#	produced by code/tidy_raw_data.R
#		data/process/toxin_tidy.tsv
#		data/process/histology_tidy.tsv
# 		data/process/metadata_tidy.tsv
#
# Nick Lesniak 2021-10-26
###############################################################################
# setup environment
###############################################################################
library(tidyverse)
library(cowplot)
library(here)

###############################################################################
# functions to load data
###############################################################################
read_taxonomy <- function(){
	read_tsv(here('data/process/final.taxonomy.tidy.tsv'),
		col_types = cols(.default = col_character()))
}

read_metadata <- function(){
	read_tsv(here('data/process/metadata_tidy.tsv'),
		col_type = 'cccccDddddDcdddddcdl') 
}

read_shared <- function(){
	read_tsv(here('data/mothur/sample.final.0.03.subsample.shared'),
		col_types = cols(.default = col_double(),
			Group = col_character())) %>% 
	mutate(Group = ifelse(grepl('_D|DA', Group), Group,
		gsub('D', '_D', Group)))
}

read_beta <- function(){
	source(here('code/read_dist.R'))
	read_dist(here('data/mothur/sample.final.thetayc.0.03.lt.ave.dist')) %>% 
		mutate(rows = ifelse(grepl('_D|DA', rows), rows,
				gsub('D', '_D', rows)),
			columns = ifelse(grepl('_D|DA', columns), columns,
				gsub('D', '_D', columns)))
}

read_toxin <- function(){
	read_tsv(here('data/process/toxin_tidy.tsv'),
		col_type = 'cccccdc')
}

read_histology <- function(){
	read_tsv(here('data/process/histology_tidy.tsv'),
		col_type = 'ccccdddcdcc')
}


donor_df <- tibble(	
	donor_colors = c("#41bbc5", "#7e2640", "#9ad859", "#682dbd", "#728e24",
		"#e26df8", "#23980d", "#fd048f", "#46e33e", "#fa1bfc",
		"#155126", "#f9b2d1", "#154975", "#e2c627", "#6294ce"),
	human_source = c('DA01324', 'DA00953', 'DA00581', 'DA01134', 'DA10148',
		'DA00369', 'DA00430', 'DA10093', 'DA10027', 'DA10034', 
		'DA10082', 'DA00884', 'DA00578', 'DA00431', 'DA01245'),
	donor_labels = factor(c(paste0('N', 1:9), paste0('M', 1:6)), 
			levels = c(paste0('N', 1:9), paste0('M', 1:6))))