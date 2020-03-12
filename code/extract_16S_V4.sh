# download all available complete C difficle genomes



gcsplit -k data/process/cdiff_16S/cdifficile_genomes.fasta '/>/' '{*}' \
	-f data/process/cdiff_16S/temp/cdifficile_genome \
	--elide-empty-files

# read in fasta file, convert to one line, pull out all matches to primer patterns

# use header of fasta to label 16S regions pulled from genome
head -n 1 data/reference/cdifficile630_reference_genome.fasta > data/process/cdiff_16S/cdifficile_16S_V4.txt
# remove new lines splitting sequences
# extract 16S sequences from each genome
# !6S V4 primers https://github.com/SchlossLab/MiSeq_WetLab_SOP/blob/master/MiSeq_WetLab_SOP.md
# 16Sf V4: GTGCCAGCMGCCGCGGTAA (reverse complement TTACCGCGGC.GCTGGCAC)
# 16Sr V4: GGACTACHVGGGTWTCTAAT (reverse complement ATTAGA.ACCC..GTAGTCC)

cat data/reference/cdifficile630_reference_genome.fasta | tr -d '\n' | 
	grep -E -o '(GTGCCAGC.GCCGCGGTAA.{200,255}ATTAGA.ACCC..GTAGTCC)|(GGACTAC..GGGT.TCTAAT.{200,255}TTACCGCGGC.GCTGGCAC)' >> data/process/cdiff_16S/cdifficile_16S_V4.txt

echo '> C. difficile Genome from DA00431' >> data/process/cdiff_16S/cdifficile_16S_V4.txt
cat data/process/cdiff_16S/consensus_cdifficile_431.fasta | tr -d '\n' | 
	grep -E -o '(GTGCCAGC.GCCGCGGTAA.{200,255}ATTAGA.ACCC..GTAGTCC)|(GGACTAC..GGGT.TCTAAT.{200,255}TTACCGCGGC.GCTGGCAC)' >> data/process/cdiff_16S/cdifficile_16S_V4.txt

genome_files="data/process/cdiff_16S/temp/*"
for temp_file in $genome_files
do
	# add genome label
	head -n 1 $temp_file >> data/process/cdiff_16S/cdifficile_16S_V4.txt
	# remove new lines splitting sequences
	cat $temp_file  | tr -d '\n' | 
		grep -E -o '(GTGCCAGC.GCCGCGGTAA.{200,255}ATTAGA.ACCC..GTAGTCC)|(GGACTAC..GGGT.TCTAAT.{200,255}TTACCGCGGC.GCTGGCAC)' >> data/process/cdiff_16S/cdifficile_16S_V4.txt

done