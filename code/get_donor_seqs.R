####
#
# read in metadata from Microbiome Data Distinguish Patients with Clostridium difficile Infection and Non-C. difficile-Associated Diarrhea from Healthy Controls. AM Schubert, et al. DOI: 10.1128/mBio.01021-14
# extraxt sequence file information 
#
####

library(tidyverse)
library(readxl)

donor_data <- read_xlsx('data/raw/MIMARKS_cdclinical.xlsx', sheet = 'Sheet1')
exp_data <- read_xlsx('data/raw/humanGF_ids.xlsx', sheet = 'complete metadata')

donor_source <- exp_data %>% 
	filter(human_source != 'AMS') %>% 
	pull(human_source) %>% 
	unique 
	
strain_source <- exp_data %>% 
	filter(!cdiff_strain %in% c('NA', 'none')) %>% 
	pull(cdiff_strain) %>% 
	unique %>% 
	paste0('DA00',.)

sample_id_list <- unique(c(donor_source, strain_source))

donor_data <- donor_data %>% 
	filter(sample_id %in% sample_id_list)

seq_file_list <- donor_data %>% 
	mutate(sfffile_id = gsub('.sff', '', sfffile_id)) %>% 
	pull(sfffile_id) %>% 
	unique

mock_list <- paste0(seq_file_list, '.MOCK') 

write_delim(donor_data, 'data/process/donor_metadata.txt', delim = '\t')
write_lines(c(sample_id_list, mock_list), 'data/process/donor_sample_ids.accnos')
write_lines(seq_file_list, 'data/process/donor_seq_files.txt')
write_lines(mock_list, 'data/process/donor_mock_ids.accnos')