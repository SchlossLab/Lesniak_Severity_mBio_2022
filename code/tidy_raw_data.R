################################################################################
#
# tidy_raw_data.R
#
# This script clean up the metadata, toxin, and histology excel files 
#
# Dependencies...
# * data/raw/humanGF_ids.xlsx
# * data/raw/Alyx_Humice_toxinassay_results.xlsx
# * data/raw/histopathology_scores_raw.xlsx
#
# Output...
# * data/process/metadata_tidy.tsv
# * data/process/toxin_tidy.tsv
# * data/process/histology_tidy.tsv
#
################################################################################

# load packages
library(tidyverse)
library(readxl)

# file inputs/outputs
input_metadata <- 'data/raw/humanGF_ids.xlsx'
input_toxin <- 'data/raw/Alyx_Humice_toxinassay_results.xlsx'
input_histology <- 'data/raw/histopathology_scores_raw.xlsx'
input_human_source <- 'data/raw/MIMARKS_cdclinical.xlsx'
output_metadata <- 'data/process/metadata_tidy.tsv'
output_toxin <- 'data/process/toxin_tidy.tsv'
output_histology <- 'data/process/histology_tidy.tsv'
output_human_source <- 'data/process/human_source_tidy.tsv'

# Read in metadata
metadata <- read_xlsx(input_metadata, 
		sheet = 'complete metadata', range = 'A1:R565',
		col_types = c(rep('text', 5), 'date', rep('numeric', 4), 'date', 'text', 
			rep('numeric', 5), 'text'), 
		na = c('NA', 'none')) %>% 
	select(-file) %>% 
	mutate(mouse_DOB = as.Date(mouse_DOB), date = as.Date(date))%>% 
	filter(day != -7) %>% # only a few mice had day -7 samples (only available for INA and OUTA)
	filter(!group %in% c('581-inoculum', 'DA581')) %>% # remove duplicate rows for inoculum DA00581
	mutate(ear_tag = mouse_id,
		   mouse_id = paste0(cage_id, '_', ear_tag), # create unique mouse id
		   group = gsub('-', '_', group)) %>% # adjust sample labels to be mothur friendly
	group_by(mouse_id) %>% 
	mutate(mouse_endpoint = max(day),
		   early_euth = ifelse(mouse_endpoint < 10, T, F)) # create column for mice which were euthanized early

# Read in toxin data
toxin_data <- read_xlsx(input_toxin, sheet = 'Sheet1',
		col_types = c('text', 'text', 'text', 'numeric')) %>% 
		rename(group = Cage_Mouse,
			   toxin_sample_day = Experiment_day,
			   toxin_sample_type = Sample_source) %>% 
		separate(group, c('cage_id', 'ear_tag', 'day'), remove = F, fill = 'left') %>% 
		mutate(mouse_id = paste0(cage_id, '_', ear_tag), # create unique mouse id
			   group = gsub('-', '_', group), # adjust sample labels to be mothur friendly
			   mouse_id = ifelse(mouse_id == 'NA_NA', group, mouse_id)) %>% # set mouse_id for inoculum
		select(-day) # eliminate column that is the same as toxin_day

# Read in histology data
histology_data <- read_xlsx(input_histology, 
		sheet = '16B037SchlossCecum', range = 'A13:L71',
		col_types = c('numeric', rep('text', 5), rep('numeric', 3), 'text', 'numeric', 'text'),
		na = c('na', '')) %>% 
	select(-outcome, -'No.') %>% # remove outcome column, previously created but not based on current analysis
	rename(cage_id = cage_ID, ear_tag = mouse_ID) %>% 
	mutate(mouse_id = paste0(cage_id, '_', ear_tag)) # create unique mouse id

# Read in human source data
human_source_data <- read_xlsx(input_human_source,
		sheet = 'Sheet1', range = 'A1:AW353') %>% 
	filter(sample_id %in% c(unique(metadata$human_source), 
							'DA00299', 'DA00395', 'DA00458')) %>% 
	select(sample_id, biome, samp_mat_process, age, race, gender, 
		   recent_antibiotic_use = "antibiotics >3mo", protonpump, h2receptor, 
		   antacid, Healthworker, historyCdiff, Surgery6mos, Vegetarian, weight, 
		   disease_stat)

#output tidied data
write_tsv(metadata, path = output_metadata)
write_tsv(toxin_data, path = output_toxin)
write_tsv(histology_data, path = output_histology)
write_tsv(human_source_data, path = output_human_source)