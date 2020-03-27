# get files if not already obtained
FILE=data/raw/MIMARKS_cdclinical.xlsx
if [ -f $FILE ]; then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist, downloading now."
   curl -o $FILE https://mothur.s3.us-east-2.amazonaws.com/data/CDI_MicrobiomeModeling/MIMARKS_cdclinical.xlsx
fi

# run rscript to create donor files
Rscript code/get_donor_seqs.R

# download sequence files
for temp_file in $(cat data/process/donor_seq_files.txt)
do
	echo Processing $temp_file
	curl -o data/raw/$temp_file.sff.bz2 https://mothur.s3.us-east-2.amazonaws.com/data/CDI_MicrobiomeModeling/$temp_file.sff.bz2
	curl -o data/raw/$temp_file.oligos https://mothur.s3.us-east-2.amazonaws.com/data/CDI_MicrobiomeModeling/$temp_file.oligos
	bzip2 -d data/raw/$temp_file.sff.bz2
done
