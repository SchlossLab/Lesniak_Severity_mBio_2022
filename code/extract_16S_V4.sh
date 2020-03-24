# download all available complete C difficle genomes
rsync --copy-links --times --verbose rsync://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Clostridioides_difficile/assembly_summary.txt data/process/cdiff_16S
Rscript get_cdiff_genome_paths.R
mkdir data/process/cdiff_16S/cdiff_genomes
for file in $(cat data/process/cdiff_16S/cdiff_genome_paths.txt)
do
	rsync --copy-links --times --verbose $file data/process/cdiff_16S/cdiff_genomes
done
gunzip data/process/cdiff_16S/cdiff_genomes/*.fna.gz

# get file list of all the C difficile refseq genomes 
# includes the genome used as the reference for assembling the 431 genome
# GCF_000009205.2_ASM920v2 >NC_009089.1 Clostridioides difficile 630, complete genome
# and the DA00431 genome
genome_list=(data/process/cdiff_16S/consensus_cdifficile_431.fasta \
 data/process/cdiff_16S/cdiff_genomes/*)

# read in fasta file, convert to one line, pull out all matches to primer patterns
# extract 16S sequences from each genome
# 16S V4 primers https://github.com/SchlossLab/MiSeq_WetLab_SOP/blob/master/MiSeq_WetLab_SOP.md
# 16Sf V4: GTGCCAGCMGCCGCGGTAA (reverse complement TTACCGCGGC.GCTGGCAC)
# 16Sr V4: GGACTACHVGGGTWTCTAAT (reverse complement ATTAGA.ACCC..GTAGTCC)
echo -e 'sequence_16S\tcdiff_genome' > data/process/cdiff_16S/cdifficile_16S_V4.txt
for temp_file in ${genome_list[@]}
do
	echo Extracting 16S V4 from $temp_file
	# get genome name from header, ( replace any / with _ )
	genome=$(head -n 1 $temp_file | tr "/" "_")
	# remove new lines splitting sequences
	cat $temp_file  | tr -d '\n' | 
		grep -E -o '(GTGCCAGC.GCCGCGGTAA.{200,255}ATTAGA.ACCC..GTAGTCC)|(GGACTAC..GGGT.TCTAAT.{200,255}TTACCGCGGC.GCTGGCAC)' | # extract 16S gene
		sed "s/$/$(printf '\t')$genome/" >> data/process/cdiff_16S/all_cdifficile_16S_V4.txt # add genome name to each line
done

rm -r data/process/cdiff_16S/cdiff_genomes
