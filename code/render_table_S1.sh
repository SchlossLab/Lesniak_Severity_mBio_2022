#!/usr/bin/env bash

# author: Nick Lesniak
# input: 
#	results/tables/Table_S1.Rmd
#	results/tables/header.tex
# outputs: submission/Table_S1.pdf

Rscript -e "rmarkdown::render('results/tables/Table_S1.Rmd')"
mv results/tables/Table_S1.pdf submission/Table_S1.pdf