################################################################################
#
#	Part 1: Get the reference files
#
#	Here we give instructions on how to get the necessary reference files that
#	are used throughout the rest of the analysis. These are used to calculate
#	error rates, generate alignments, and classify sequences. We will pull down
#	the mock community reference (HMP_MOCK.fasta), the silva reference alignment
#	(silva.bacteria.align), and the RDP training set data (trainset18_062020).
#	Finally, we use the HMP_MOCK.align to get the alignment coordinates for the
#	V4 data. These data will be stored in the data/references/ folder.
#
#	The targets in this part are all used as dependencies in other rules
#
################################################################################

#Location of files
REFS = data/references/
RAW = data/raw/
MOTHUR = data/mothur/
PROC = data/process/

# create all data and reference files
#get v4 region of the silva reference alignment, rdp training set data and hmp mock
$(REFS)silva.v4.align $(REFS)trainset18_062020.v4.fasta $(REFS)trainset18_062020.v4.tax $(REFS)HMP_MOCK.v4.align : code/get_references.batch
	bash code/get_references.batch

references : $(REFS)silva.v4.align $(REFS)trainset18_062020.v4.fasta $(REFS)trainset18_062020.v4.tax $(REFS)HMP_MOCK.v4.align

################################################################################
#
#	Part 2: Get fastq files
#
################################################################################

$(RAW)MIMARKS_cdclinical.xlsx : 
	wget -P $(RAW) https://mothur.s3.us-east-2.amazonaws.com/data/CDI_MicrobiomeModeling/MIMARKS_cdclinical.xlsx 

# Clean metadata for this project
$(PROC)metadata_tidy.tsv $(PROC)toxin_tidy.tsv $(PROC)histology_tidy.tsv $(PROC)human_source_tidy.tsv : $(RAW)humanGF_ids.xlsx\
										$(RAW)Alyx_Humice_toxinassay_results.xlsx\
										$(RAW)histopathology_scores_raw.xlsx\
										$(RAW)MIMARKS_cdclinical.xlsx\
										code/tidy_raw_data.R
	Rscript code/tidy_raw_data.R

# get the fastq files
$(RAW)get_data : code/get_fastqs.sh $(MOTHUR)SRR_Acc_List.txt
	bash code/get_fastqs.sh\
	touch $(RAW)get_data

# build the files file. 
$(MOTHUR)humanGF_cdiff.files : code/make_files_file.R $(RAW)human_source_tidy.tsv
	Rscript code/make_files_file.R


################################################################################
#
#	Part 3: Run data through mothur 
#
################################################################################

SUB = 2107

# install mothur in the project code directory
code/mothur/mothur : code/install_mothur.sh
	sh code/code/install_mothur.sh

# here we go from the raw fastq files and the files file to generate a fasta,
# taxonomy, and count_table file that has had the chimeras removed as well as
# any non bacterial sequences
# then we get the sequencing error as seen in the mock community samples
# then we go from the good sequences and generate a shared file and a cons.taxonomy 
# file based on OTU data. Lastly, we rarefy the number of reads to $(SUB) sequences 
# per library for the alpha and beta diversity analyses and modeling and 
# for the shared file
$(MOTHUR)complete.sample.final.shared $(MOTHUR)sample.final.shared $(PROC)mock.sample.final.shared $(MOTHUR)final.taxonomy $(MOTHUR)sample.error.count : code/get_good_seqs.batch\
										code/get_error.batch\
										code/get_shared_otus.batch\
										$(REFS)silva.v4.align\
										$(REFS)trainset18_062020.v4.fasta\
										$(REFS)trainset18_062020.v4.tax\
										$(REFS)HMP_MOCK.v4.fasta
	mothur code/get_good_seqs.batch
	mothur code/get_error.batch
	mothur code/get_shared_otus.batch
	mv data/mothur/humanGF_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared data/mothur/complete.sample.final.shared
	code/mothur/mothur "#remove.groups(shared=data/mothur/complete.sample.final.shared, groups=mock_S282-mock_S343-mock_S344-mock_S102-mock_S96)"
	mv data/mothur/complete.sample.final.0.03.pick.shared data/mothur/sample.final.shared
	code/mothur/mothur "#set.dir(input=data/mothur, output=data/mothur);
	summary.single(shared=sample.final.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=$(SUB));
	dist.shared(shared=sample.final.shared, calc=thetayc-jclass, subsample=$(SUB));
	nmds(phylip=current);
	sub.sample(shared=sample.final.shared, size=$(SUB));
	get.groups(shared=complete.sample.final.shared, groups=mock_S282-mock_S343-mock_S344-mock_S102-mock_S96)"
	mv data/mothur/complete.sample.final.0.03.pick.shared data/process/mock.sample.final.shared
	mv data/mothur/humanGF_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy data/mothur/final.taxonomy
	mv data/mothur/humanGF_cdiff.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.error.count data/mothur/sample.error.count

################################################################################
#
#	Part 4: Write the paper
#
################################################################################


################################################################################
# Process raw data for analysis
################################################################################

# Convert taxonomy file into a dataframe with OTU labels
$(PROC)final.taxonomy.tidy.tsv : $(MOTHUR)final.taxonomy\
									code/tidy_taxonomy.R
	Rscript code/tidy_taxonomy.R


################################################################################
#
# Run LEfSE and RF
#
################################################################################

# Run LEfSe analysis
$(PROC)lefse/day_0_severity_Genus.0.03.lefse_summary $(PROC)lefse/day_0_severity_Genus.design $(PROC)lefse/day_0_toxin_Genus.0.03.lefse_summary $(PROC)lefse/day_0_toxin_Genus.design : $(PROC)metadata_tidy.tsv\
										$(MOTHUR)sample.final.0.03.subsample.shared\ 
										$(PROC)toxin_tidy.tsv\
										$(PROC)histology_tidy.tsv\
										$(PROC)final.taxonomy.tidy.tsv\
										code/run_lefse_analysis.R
	Rscript code/run_lefse_analysis.R

# Run LEfSe on temporal trends
$(PROC)lefse/temporal_trend_hilo.0.03.lefse_summary $(PROC)lefse/temporal_trend_hilo.design : $(PROC)metadata_tidy.tsv\
										$(MOTHUR)sample.final.0.03.subsample.shared\ 
										$(PROC)toxin_tidy.tsv\
										$(PROC)histology_tidy.tsv\
										$(PROC)final.taxonomy.tidy.tsv\
										code/run_lefse_analysis.R
	Rscript code/run_lefse_temporal_trend.R

# Run random forest to predict C. difficile challenge outcome from day 0 community
$(PROC)ml/lr_performance.tsv $(PROC)ml/lr_feature_imp.tsv $(PROC)ml/lr_hp_performance.tsv : $(PROC)metadata_tidy.tsv\
										$(PROC)toxin_tidy.tsv\
										$(PROC)histology_tidy.tsv\
										$(MOTHUR)sample.final.0.03.subsample.shared\
										code/ml/setup_ml_data.R\
										$(PROC)ml/day_0*\
										code/ml/tune_ml.R\
										code/ml/run_ml.R\
										$(PROC)ml/temp/*\
										code/ml/concat_ml_results.sh\
										code/ml/concat_tune_ml_results.sh
	for seed in {1..100}
	do
		for tax_rank in {Phylum,Class,Order,Family,Genus,OTU}
		do
			echo 'Rscript code/ml/run_ml.R' $seed $tax_rank
		done
	done
	# instead of for loops, run SBATCH code/ml/run_ml_{tax_rank}.sbat and SBATCH code/ml/tune_ml_{tune_frac}.sbat
	# combine the output from all the seeds and runs
	bash code/ml/concat_ml_results.sh
	bash code/ml/concat_tune_ml_results.sh


################################################################################
# Create figures
################################################################################

# Figure 1
submission/Figure_1.tiff : code/plot_figure_1.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)final.taxonomy.tidy.tsv\
							$(PROC)sample.final.thetayc.0.03.lt.ave.dist\
							code/read.dist.R
	Rscript code/plot_figure_1.R

# Figure 2
submission/Figure_2.tiff : code/plot_figure_2.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv
	Rscript code/plot_figure_2.R

# Figure 3 and S1
submission/Figure_3.tiff submission/Figure_S1.tiff : code/plot_figure_3.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(PROC)toxin_tidy.tsv\
							$(PROC)histology_tidy.tsv
	Rscript code/plot_figure_3.R

# Figure 4
submission/Figure_4.tiff : code/plot_figure_4.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)final.taxonomy.tidy.tsv\
							$(PROC)histology_tidy.tsv\
							$(PROC)lefse/day_0_severity_Genus.0.03.lefse_summary\
							$(PROC)lefse/day_0_severity_Genus.design\
							$(PROC)lefse/day_0_toxin_Genus.0.03.lefse_summary\
							$(PROC)lefse/day_0_toxin_Genus.design
	Rscript code/plot_figure_4.R

# Figure 5
submission/Figure_5.tiff : code/plot_figure_5.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)final.taxonomy.tidy.tsv\
							$(PROC)histology_tidy.tsv\
							$(PROC)toxin_tidy.tsv\
							$(PROC)ml/ml_feature_imp.tsv
	Rscript code/plot_figure_5.R

# Figure S2
submission/Figure_S2.tiff : code/plot_figure_S2.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
	Rscript code/plot_figure_S2.R

# Figure S3
submission/Figure_S3.tiff : code/plot_figure_S3.R\
							$(PROC)ml/lr_performance.tsv\
							$(PROC)ml/lr_feature_imp.tsv
	Rscript code/plot_figure_S3.R

# Figure S4
submission/Figure_S4.tiff : code/plot_figure_S4.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)final.taxonomy.tidy.tsv\
							$(PROC)histology_tidy.tsv\
							$(PROC)/lefse/temporal_trend_hilo.0.03.lefse_summary\
							$(PROC)/lefse/temporal_trend_hilo.design
	Rscript code/plot_figure_S3.R

# Table S1
submission/Table_S1.pdf : code/render_table_S1.sh\
							results/table/Table_S1.Rmd\
							results/table/header.tex
	bash code/render_table_S1.sh

################################################################################
# Create manuscript
################################################################################


submission/manuscript.pdf submission/manuscript.docx : submission/manuscript.Rmd\
		submission/figure_1.tiff\
		submission/figure_2.tiff\
		submission/figure_3.tiff\
		submission/figure_4.tiff\
		submission/figure_5.tiff\
		submission/figure_S1.tiff\
		submission/figure_S2.tiff\
		submission/figure_S3.tiff\
		submission/figure_S4.tiff\
		submission/mbio.csl\
		submission/references.bib
	R -e 'library(rmarkdown);render("submission/manuscript.Rmd", output_format="all")'

