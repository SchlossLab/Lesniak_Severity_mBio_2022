Otu000033	3061	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu000152	8701	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(98);
Otu000153	29231	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu000628	2805	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu001100	228	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu001148	13	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(93);
Otu001149	108	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu001521	9	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu001824	209	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Peptostreptococcus(100);
Otu001845	7	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu002284	3	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu002316	4	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu003821	9	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);
Otu003946	4	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu004148	2	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu004423	2	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu004549	1	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu005188	3	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu005234	2	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);
Otu005391	1	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);
Otu005403	1	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);Clostridium_XI(100);
Otu005577	3	Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Peptostreptococcaceae(100);

Peptostreptococcaceae_OTUs <- c('Otu000033', 'Otu000152', 'Otu000153', 'Otu000628')
	# Following OTUs have a single sample with a count of 1,1,0,0,0,1,3,1 repsectively
	# 'Otu001100', 'Otu001148', 'Otu001149', 'Otu001521', 'Otu001824', 'Otu001845', 'Otu002284', 'Otu002316', 'Otu003821', 'Otu003946', 'Otu004148', 'Otu004423', 'Otu004549', 'Otu005188', 'Otu005234', 'Otu005391', 'Otu005403', 'Otu005577')

shared_file %>% 
	mutate(sample = row.names(shared_file)) %>% 
	gather(OTU, counts, contains('Otu')) %>% 
#	group_by(sample) %>% 
#	summarize(sum(counts)) # check subsample level (2000)
	filter(OTU %in% c('Otu000033', 'Otu000152', 'Otu000153', 'Otu000628', 'Otu001100', 'Otu001149', 'Otu001521', 'Otu001845', 'Otu002284', 'Otu002316', 'Otu003946', 'Otu004148', 'Otu004423', 'Otu004549', 'Otu005188', 'Otu005403')) %>% 
	group_by(sample) %>% 
	summarize(counts = sum(counts)) %>% 
	filter(counts > 1, grepl('-D0', sample)) %>% 
	data.frame %>% 
	left_join(select(mutate(meta_file, sample = row.names(meta_file), human_source = as.character(human_source)), sample, human_source)) %>% 
	left_join(donor_counts)

donor_counts <- shared_file %>% 
	mutate(human_source = row.names(shared_file)) %>% 
	filter(human_source %in% as.character(unique(meta_file$human_source))) %>% 
	gather(OTU, counts, contains('Otu')) %>% 
	filter(OTU %in% c('Otu000033', 'Otu000152', 'Otu000153', 'Otu000628', 'Otu001100', 'Otu001149', 'Otu001521', 'Otu001845', 'Otu002284', 'Otu002316', 'Otu003946', 'Otu004148', 'Otu004423', 'Otu004549', 'Otu005188', 'Otu005403')) %>% 
	group_by(human_source) %>% 
	summarize(cdiff_donor_counts = sum(counts))


Does having C. difficile 16s on day 0 mean they are already infected? or is it possible these donors were asypmtomatic carriers?
Do we have plating info of day 0?

select(meta_file, human_source, cage_id) %>% 
	unique

histo_file <- 'data/process/histology_scores_txt.txt'
histo_file <- read.table(histo_file, sep = '\t', header = T)
histo_file %>% 
	mutate(severity = case_when(day < 5 ~ 'severe',
			summary_score > 5 ~ 'mild',
			summary_score < 5 ~ 'persistent',
			T ~ 'NA'),
			timepoint = case_when(severity == 'severe' ~ 'Early Deceased (Prior to Day 4)',
				T ~ 'End Point')) %>% select(cage_ID, severity) %>% 
				unique
	filter(severity != 'NA') %>% 
	mutate(severity = factor(severity, levels = c('persistent', 'mild', 'severe')),
		timepoint = factor(timepoint, levels = c('End Point', 'Early Deceased (Prior to Day 4)'))) %>% 
	ggplot(aes(x = severity, y = summary_score, color = severity)) +
		geom_jitter() + 
		facet_grid(.~timepoint, scales = 'free_x', space = 'free_x')		


What are C diff counts on days prior to day 0?
What is the distribution Peptostreptococcaceae OTUs over the experiment?
Is there one Peptostreptococcaceae OTU associated with the Cdiff used for colonization?
what portion of day 0 samples have 16S counts?

shared_file %>% 
	mutate(sample = row.names(shared_file)) %>% 
	gather(Cdiff_OTU, counts, contains('Otu')) %>% 
	filter(Cdiff_OTU %in% Peptostreptococcaceae_OTUs) %>% 
	left_join(mutate(meta_file, sample = row.names(meta_file)), by = 'sample') %>% 
	mutate(day = case_when(grepl('DA', sample) ~ as.numeric(-2.0),
		is.na(day) ~ as.numeric(gsub('.*D', '', sample)),
		T ~ as.numeric(day)),
		mouse_id = case_when(grepl('DA', sample) ~ as.character(sample),
			is.na(mouse_id) ~ gsub('-D\\d*', '', sample),
			T ~ as.character(mouse_id))) %>% 
	ggplot(aes(x = day, y = counts, group = mouse_id)) + 
#	ggplot(aes(x = day, y = log10(cdiff_cfu + 60), color = OTU, group = mouse_id)) + 
		geom_line(alpha = 0.4) + geom_point() + 
		facet_wrap(.~Cdiff_OTU, scales = 'free_y') + 
		theme_bw() +
		labs(title = 'Clostridium XI counts \nSubsampled at 2000, Day -2 Samples are Inocula')

16 donors, 4 with 16S (not one 431 was isolated from)
66 mice, 

donors
grep('DA', row.names(shared_file), value = T)
"DA10027" "DA10034" "DA10082" "DA10093" "DA10148" "DA01134" "DA01146"
"DA01245" "DA01324" "DA00369" "DA00430" "DA00431" "DA00578" "DA00581"
"DA00884" "DA00953"

strains
"ribotype","severe","tcdA","tcdB","cdtA","cdtB","tcdC","sporulation","toxin","CFE","growth","germ_tc_only","germ_tc_and_gly"
"DA00299","UM11","Y","+","+","-","-","ND","33900","0.173","52.98","0.1576",0.28,94.93
"DA00395","027","Y","+","+","+","+","ND","16433","2.780","42.79","0.0645",0.08,46.71
"DA00431","027","Y","+","+","+","+","ND","64000","3.053","39.56","0.0637",0.04,80.64
"DA00458","027","Y","+","+","+","+","ND","27850","2.993","55.19","0.0777",0.07,53.43

Download the Cdiff sequence of the one used for the infections
	DA00431, assemble genome and extract 16S V4
Compare 16S V4 of Cdifficile
	In all C difficile genomes available (85) on refseq, there is very little variation in the 16S V4 region. There are a total of 819 V4 regions ranging from 1 to 16 per genome with the median/mean of 12.  792 are identical and the remaining are 98.8% similar or greater except for one sequence that is 95.6% similar. This outlier sequence is from C. difficile M68 ribotype 017 and has 15 other 16S V4 genes that are the consensus sequence.
C difficile sequence stored in data/process/cdiff_16S/consensus_cdifficile_16S_V4.txt
	TACGTAGGGGGCTAGCGTTATCCGGATTTACTGGGCGTAAAGGGTGCGTAGGCGGTCTTTCAAGTCAGGAGTGAAAGGCTACGGCTCAACCGTAGTAAGCTCTTGAAACTGGGAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTTGCGAAGGCGGCTCTCTGGACTGTAACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG
Mapped it to sequence >M00967_67_000000000-A4RV6_1_1104_20443_12527
	in file humanGF_cdiff.trim.contigs.good.unique.fasta
*** Reran mothur so OTU numbers have changed ***
#clust_list <- read.table('data/mothur/humanGF_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list')
#clust_list <- data.frame(t(clust_list)) # convert to columns
#colnames(clust_list) <- c('otu', 'sequences')
#cdiff_otu <- clust_list[grepl('M00967_67_000000000-A4RV6_1_1104_20443_12527', clust_list$sequences),'otu']
# Otu0013
#otu13_sequences <- clust_list[grepl('M00967_67_000000000-A4RV6_1_1104_20443_12527', clust_list$sequences),'sequences']
#otu13_sequences <- gsub('^ *', '', otu13_sequences) # remove leading spaces
#otu13_sequences <- gsub(',', '\n', otu13_sequences) # replace commas wth a newline
#write.table(otu13_sequences, 'data/process/cdiff_16S/OTU00013_sequences.txt',
#	quote = F, col.names = F, row.names = F)
This sequence ends up in OTU 13 with 427 other sequences
	saved as data/process/cdiff_16S/OTU00013_sequences.txt
Pull in data/mothur/humanGF_cdiff.trim.contigs.good.unique.fasta
Select sequences and compare to C difficile
#awk '{printf "%s%s",$0,(NR%2?FS:RS)}' data/mothur/humanGF_cdiff.trim.contigs.good.unique.fasta > data/process/cdiff_16S/humanGF_cdiff.wide.fasta
#library(tidyverse);library(msa)
#all_fasta <- read.table('data/process/cdiff_16S/humanGF_cdiff.wide.fasta')
#colnames(all_fasta) <- c('name', 'ee', 'sequence')
#cdiff_concensus_seq <- read.table("data/process/cdiff_16S/consensus_cdifficile_16S_V4.txt", stringsAsFactors = F)
#otu13_sequence_names <- read.table("data/process/cdiff_16S/OTU00013_sequences.txt", stringsAsFactors = F)
#otu13_sequence_df <- otu13_sequence_names %>% 
#	mutate(name = paste0('>', V1)) %>% 
#	select(name) %>% 
#	left_join(all_fasta) %>% 
#	mutate(sequence = as.character(sequence))

# determine count distribution of sequences in OTU 13
counts_by_seqeunce <- read.table('data/mothur/humanGF_cdiff.trim.contigs.good.good.count_table', 
	stringsAsFactors = F, skip = 2, fill = T)[-1,1:2]
colnames(counts_by_seqeunce) <- c('Representative_Sequence', 'total')
counts_by_seqeunce <- counts_by_seqeunce %>% 
	filter(Representative_Sequence %in% unlist(otu13_sequence_names)) %>% 
	mutate(total = as.numeric(total))
n_cdiff <- counts_by_seqeunce %>% 
	filter(Representative_Sequence == 'M00967_67_000000000-A4RV6_1_1104_20443_12527') %>% 
	pull(total)
n_wo_cdiff <- counts_by_seqeunce %>% 
	filter(Representative_Sequence != 'M00967_67_000000000-A4RV6_1_1104_20443_12527') %>% 
	pull(total) %>% 
	sum
counts_by_seqeunce
# determine the how similar the other clustered sequences are to the consensus C. difficile sequence
otu13_sequence_identity_df <- 
	map_dfr(otu13_sequence_df$sequence, function(i){ 
		p_aln <- pairwiseAlignment(cdiff_concensus_seq$V1, i)
		percent_alignment <- pid(p_aln)	
		return(tibble(sequence = i, percent_identity = percent_alignment))
	})
otu13_sequence_identity_df <- otu13_sequence_identity_df %>% 
	full_join(otu13_sequence_df) %>% 
	mutate(name = gsub('^>', '', name)) %>% 
	full_join(counts_by_seqeunce, by = c('name' = 'Representative_Sequence'))
library(cowplot)
counts_plot <- otu13_sequence_identity_df %>% 
	ggplot(aes(percent_identity)) + geom_histogram(binwidth = .1) + 
		labs(x = 'Percent Identity to C. difficile Consensus Sequence', y = 'Count of unique sequences', 
			title = 'OTU 13 - Peptostreptococcaceae_unclassified sequences',
			subtitle = '428 Sequences that cluster with the C. difficile consensus sequence (1)') +
		coord_cartesian(xlim = c(96,100.2)) + theme_bw() + coord_flip() + scale_y_reverse()

sequences_plot <- otu13_sequence_identity_df %>% 
	group_by(percent_identity) %>% 
	summarise(`16S` = sum(as.numeric(total))) %>% 
	filter(!is.na(percent_identity)) %>% 
	ggplot(aes(x = percent_identity, y = `16S`)) + 
		geom_bar(stat='identity', width = .125, alpha = .5) + 
		labs(x = 'Percent Identity to C. difficile Consensus Sequence', y = 'OTU 13 Sequence Counts', 
			title = 'OTU 13 - Peptostreptococcaceae_unclassified sequences',
			subtitle = '428 Sequences that cluster with the C. difficile consensus sequence (1)') +
				coord_cartesian(xlim = c(96,100.2)) + theme_bw() + coord_flip() + scale_y_log10(
   			breaks = scales::trans_breaks("log10", function(x) 10^x),
   			labels = scales::trans_format("log10", scales::math_format(10^.x)))# + 
   		#annotation_logticks(sides = 'b')
# using plot of cdiff v4 sequences from code/get_cdiff_variation.R
ggsave('~/Desktop/cdiff_otu.jpg',
	plot_grid(counts_plot, sequences_plot, cdiff_v4_seq_dist_plot, nrow = 1),
	width = 12, height = 8)


otu13_aligned <- msa(otu13_sequence_df$sequence, type = 'dna')
# for alignment and logo plot use following code
msaPrettyPrint(otu13_aligned, output = 'pdf')	
# compare sequence similarity to cdiff

