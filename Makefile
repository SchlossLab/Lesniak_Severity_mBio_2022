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

# Clean metadata for this project
$(PROC)metadata_tidy.tsv $(PROC)toxin_tidy.tsv $(PROC)histology_tidy.tsv $(PROC)human_source_tidy.tsv : $(RAW)humanGF_ids.xlsx\
										$(RAW)Alyx_Humice_toxinassay_results.xlsx\
										$(RAW)histopathology_scores_raw.xlsx\
										$(RAW)MIMARKS_cdclinical.xlsx\
										code/tidy_raw_data.R
	Rscript code/tidy_raw_data.R

# get the fastq files
#### UPDATE USING SRA to download fastqs into data/mothur ####
#$(RAW)get_data : code/get_fastqs.sh $(MOTHUR)humanGF_cdiff.files
#	bash code/get_fastqs.sh $(MOTHUR)humanGF_cdiff.files;\
#	touch $(RAW)get_data

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
#$(PROC)lefse/temporal_trend_hilo.0.03.lefse_summary $(PROC)lefse/temporal_trend_hilo.design : $(PROC)metadata_tidy.tsv\
#										$(MOTHUR)sample.final.0.03.subsample.shared\ 
#										$(PROC)toxin_tidy.tsv\
#										$(PROC)histology_tidy.tsv\
#										$(PROC)final.taxonomy.tidy.tsv\
#										code/run_lefse_analysis.R
#	Rscript code/run_lefse_temporal_trend.R

####SEARCH_DIR=data/temp/l2_otu
####FINAL_DIR=data/process/l2_otu
##### Create dataframe of subset samples classification column (cleared), and features (OTUs)
##### and create correlation matrix of features
####data/process/otu_input_data.csv data/process/otu_sample_names.txt data/process/sig_flat_corr_matrix_otu.csv : code/R/setup_model_data.R\
####																	code/R/compute_correlation_matrix.R\
####																	data/process/abx_cdiff_metadata_clean.txt\
####																	data/process/abx_cdiff_taxonomy_clean.tsv\
####																	data/mothur/sample.final.0.03.subsample.shared\
####
####	Rscript code/R/setup_model_data.R l2_otu
####
##### Run pipeline array
####$SEARCH_DIR/walltime_L2_Logistic_Regression_1.csv : code/R/main.R\
####						code/R/run_model.R\
####						code/R/model_pipeline.R\
####						code/R/tuning_grid.R\
####						code/R/permutation_importance.R\
####						code/run_logit.sbat\
####						code/R/auprc.R\
####						code/R/functions.R
####	for seed in {1..100}
####	do
####		Rscript code/R/main.R --seed $seed --model L2_Logistic_Regression --level l2_otu --data  data/process/l2_otu_input_data.csv --hyperparams data/default_hyperparameters.csv --outcome clearance --permutation
####	done
####	# or run SBATCH code/run_logit.sbat on the Great Lakes cluster
####
##### concatenate results and determine feature importance
####$FINAL_DIR/combined_all_imp_features_cor_results_L2_Logistic_Regression.csv : code/bash/process_l2_output.sh\
####																				code/R/get_feature_rankings.R
####	bash code/bash/process_l2_output.sh


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

# Figure 3
submission/Figure_3.tiff : code/plot_figure_3.R\
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

# Figure 6
submission/Figure_6.tiff : code/plot_figure_6.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(PROC)toxin_tidy.tsv\
							$(PROC)histology_tidy.tsv
	Rscript code/plot_figure_6.R

# Figure S1
submission/Figure_S1.tiff : code/plot_figure_S1.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
	Rscript code/plot_figure_S1.R

# Figure S2
submission/Figure_S2.tiff : code/plot_figure_S2.R\
							$(PROC)ml/ml_performance.tsv\
							$(PROC)ml/ml_feature_imp.tsv
	Rscript code/plot_figure_S2.R

# Figure 5
submission/Figure_S3.tiff : code/plot_figure_S3.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)final.taxonomy.tidy.tsv\
							$(PROC)ml/ml_feature_imp.tsv
	Rscript code/plot_figure_S3.R

# Figure S4
submission/Figure_S4.tiff : code/plot_figure_S4.R\
							code/utlities.R\
							$(PROC)metadata_tidy.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)final.taxonomy.tidy.tsv\
							$(PROC)histology_tidy.tsv\
							$(PROC)lefse/temporal_trend_hilo.0.03.lefse_summary\
							$(PROC)lefse/temporal_trend_hilo.design
	Rscript code/plot_figure_S4.R

################################################################################
# Create manuscript
################################################################################


submission/manuscript.pdf submission/manuscript.docx : submission/manuscript.Rmd\
		submission/figure_1.tiff\
		submission/figure_2.tiff\
		submission/figure_3.tiff\
		submission/figure_4.tiff\
		submission/figure_5.tiff\
		submission/figure_6.tiff\
		submission/figure_S1.tiff\
		submission/figure_S2.tiff\
		submission/figure_S3.tiff\
		submission/mbio.csl\
		submission/references.bib
	R -e 'library(rmarkdown);render("submission/manuscript.Rmd", output_format="all")'

