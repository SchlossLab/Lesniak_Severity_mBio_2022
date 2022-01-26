#!/bin/bash

# Author: Nick Lesniak
# Date: 2021-09-16
#
#######################################################################################
# This script will:
#	Combine them together to have the results for 100 datasplits in one file
#   
########################################################################################

SEARCH_DIR=data/process/ml/temp
FINAL_DIR=data/process/ml/

# Keep the first line of File1 and remove the first line of all the others and combine
# combine binary models
for taxa_level in "OTU" "Genus" "Family" "Order" "Class" "Phylum"
	do
	head -1 $SEARCH_DIR/tune_frac_50_ml_performance_"$taxa_level"_1.tsv  > $SEARCH_DIR/combined_tune_ml_performance_"$taxa_level".tsv; tail -n +2 -q $SEARCH_DIR/tune_frac_+([0-9])_ml_performance_"$taxa_level"_+([0-9]).tsv >> $SEARCH_DIR/combined_tune_ml_performance_"$taxa_level".tsv
	head -1 $SEARCH_DIR/tune_frac_50_ml_hp_performance_"$taxa_level"_1.tsv  > $SEARCH_DIR/combined_tune_ml_hp_performance_"$taxa_level".tsv; tail -n +2 -q $SEARCH_DIR/tune_frac_+([0-9])_ml_hp_performance_"$taxa_level"_+([0-9]).tsv >> $SEARCH_DIR/combined_tune_ml_hp_performance_"$taxa_level".tsv

	mv $SEARCH_DIR/combined_tune_ml_performance_"$taxa_level".tsv $FINAL_DIR/combined_tune_ml_performance_"$taxa_level".tsv
	mv $SEARCH_DIR/combined_tune_ml_hp_performance_"$taxa_level".tsv $FINAL_DIR/combined_tune_ml_hp_performance_"$taxa_level".tsv
done

head -1 $FINAL_DIR/combined_tune_ml_performance_Phylum.tsv  > $FINAL_DIR/ml_tune_performance.tsv; tail -n +2 -q $FINAL_DIR/combined_tune_ml_performance_*.tsv >> $FINAL_DIR/ml_tune_performance.tsv
head -1 $FINAL_DIR/combined_tune_ml_hp_performance_Phylum.tsv  > $FINAL_DIR/ml_tune_hp_performance.tsv; tail -n +2 -q $FINAL_DIR/combined_tune_ml_hp_performance_*.tsv >> $FINAL_DIR/ml_tune_hp_performance.tsv