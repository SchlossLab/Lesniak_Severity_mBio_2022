################################################################################
#
# setup_sra_submission.R
#
# This script creates the files for submitting to the sra 
#
# Dependencies...
# * data/raw/*fastq
# * data/mothur/humanGF_cdiff.files
# * data/process/metadata_tidy.tsv
# * data/mothur/cdiff_severity.project 
#		see https://mothur.org/wiki/creating_a_new_submission/ for project file
#
# Output...
# * data/mothur/humanGF_cdiff.tsv
# * data/mothur/submission.xml
# * data/mothur/submit.ready
#
################################################################################

# Prior to running this, 
# Email the SRA at NCBI (sra\@ncbi.nlm.nih.gov) to notify your will be submitting mothur created files via ftp
# Create a .project file with necessary information

# load packages
library(tidyverse)

# Setup the MIMarks file
system('code/mothur/mothur "#get.mimarkspackage(inputdir=data/raw, file=data/mothur/humanGF_cdiff.files, package=host_associated, requiredonly=t)"')

# Read in MIMarks file
mimarks <- read_tsv('data/mothur/humanGF_cdiff.tsv', skip = 9) %>% 
	select(sample_name = `*sample_name`) %>% 
	mutate(group = sample_name,
		metadata_label = gsub('(\\d)D', '\\1_D', sample_name))
# Required columns
# "sample_name"     "description"     "sample_title"    "seq_methods"    
# "organism"        "collection_date" "env_biome"       "env_feature"    
# "env_material"    "geo_loc_name"    "host"            "lat_lon"

# Add metadata and required information
sample_metadata <- read_tsv('data/process/metadata_tidy.tsv') %>% 
	rename(metadata_label = group,
		human_fecal_donor = human_source,
		collection_date = date) %>% 
	select(-nextDayCFU_colon, -nextDayCFU_any, -percent_weightLoss) %>% 
	inner_join(mimarks, by = 'metadata_label') %>% 
	mutate(description = paste('Day', day, ',relative to C. difficile challenge, fecal sample from mouse', 
			ear_tag, 'cage', cage_id, ', which received its gut community from human donor', human_fecal_donor, 
			'and was challenged with C. difficile strain', cdiff_strain),
		sample_title = paste('Mouse', mouse_id, ', Day', day),
		seq_methods = 'http://www.mothur.org/w/images/0/0c/Wet-lab_MiSeq_SOP.pdf ',
		organism = 'mouse gut metagenome',
		env_biome = 'mouse gut',
		env_feature = 'feces',
		env_material = 'feces',
		geo_loc_name = 'USA: Michigan, Ann Arbor',
		host = 'Mus musculus',
		lat_lon = '42.2820973 N 83.7338904 W')

mock_metadata <- mimarks %>% 
	filter(grepl('ock', metadata_label)) %>% 
	mutate(description = 'mock community',
		sample_title = 'mock community',
		seq_methods = 'http://www.mothur.org/w/images/0/0c/Wet-lab_MiSeq_SOP.pdf ',
		organism = 'synthetic metagenome',
		collection_date = as.Date("2014-01-01"),
		env_biome = 'mock',
		env_feature = 'mock',
		env_material = 'mock',
		geo_loc_name = 'USA: Michigan, Ann Arbor',
		host = 'mock',
		lat_lon = '42.2820973 N 83.7338904 W')
# save mimarks data
write_tsv(relocate(bind_rows(sample_metadata, mock_metadata), group),
	'data/mothur/TEMP_humanGF_cdiff.tsv')
# add mimarks header back
system('head -9 data/mothur/humanGF_cdiff.tsv | cat - data/mothur/TEMP_humanGF_cdiff.tsv > data/mothur/TEMP.tsv')
system('mv data/mothur/TEMP.tsv data/mothur/humanGF_cdiff.tsv')
system('rm data/mothur/TEMP_humanGF_cdiff.tsv')
# Create xml file for submitting fastqs
system('code/mothur/mothur "#make.sra(inputdir=data/raw, file=data/mothur/humanGF_cdiff.files, project=data/mothur/cdiff_severity.project, mimark=data/mothur/humanGF_cdiff.tsv)"')
system('touch data/mothur/submit.ready')
# On an interactive node,
# Follow NCBI instructions to submit
#	data files
#	submission.xml file
#	an empty file called "submit.ready"
# via ftp @ https://www.ncbi.nlm.nih.gov/sra/docs/submitfiles/