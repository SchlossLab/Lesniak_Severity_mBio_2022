#!/bin/bash

################################################################################
#
# get_fastqs.sh
#
# Dependencies...
#   SRA Toolkit installed
#   Accession List file downloaded
#       data/mothur/SRR_Acc_List.txt
#
# Output...
#   data/mothur/*.fastq.gz
#
################################################################################

################################################################################
# Setup enviroment for download
################################################################################

# Following NIH SRA download guide - https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
#   ensure SRA Toolkit is installed by running:
which prefetch

# Download Accession List file
#   Goto https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP350240&o=acc_s%3Aa
#   Download the Accession List file, which is a text file of all the SRR run ids
#   Move Accession List file into /data/mothur

################################################################################
# Prepare directory for download
################################################################################

# set location for sequence files
DATA=data/mothur

# remove fastq files from raw dir
rm $DATA/*.fastq

# check for Accession List file
if test -f "$DATA/SRR_Acc_List.txt"; then
    echo "$DATA/SRR_Acc_List.txt exists"
else
    echo "ERROR: $DATA/SRR_Acc_List.txt is missing. Check dependencies in code/get_fastqs.sh"
    exit
fi

################################################################################
# Download sequence data
################################################################################

# For each SRR run file,
# Prefetch, download, split and gzip each fastq 
for sample in `cat $DATA/SRR_Acc_List.txt`
do
  echo $sample
    prefetch $sample
  fasterq-dump --split-files $sample -O $DATA

    gzip -f $DATA/${sample}_1.fastq
    gzip -f $DATA/${sample}_2.fastq
done