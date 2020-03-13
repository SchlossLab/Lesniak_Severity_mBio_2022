# Download DA00431 Cdiff sequence
fastq-dump --split-files --origfmt --gzip SRR804616 -O data/raw
# Download Cdiff reference genome
curl -L -o data/reference/cdifficile630_reference_genome.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/205/GCF_000009205.1_ASM920v1/GCF_000009205.1_ASM920v1_genomic.fna.gz
gunzip data/reference/cdifficile630_reference_genome.fasta.gz

# Get read quality
fastqc data/raw/SRR5804616_*.fastq.gz -o data/process/cdiff_16S/

# Clean reads
#cp  ~/bin/Trimmomatic-0.38/adapters/NexteraPE-PE.fa data/reference
trimmomatic PE data/raw/SRR5804616_1.fastq.gz data/raw/SRR5804616_2.fastq.gz \
	data/process/cdiff_16S/SRR5804616_1.trim.fastq.gz data/process/cdiff_16S/SRR5804616_1un.trim.fastq.gz \
	data/process/cdiff_16S/SRR5804616_2.trim.fastq.gz data/process/cdiff_16S/SRR5804616_2un.trim.fastq.gz \
	ILLUMINACLIP:data/reference/NexteraPE-PE.fa:2:40:15 \
	SLIDINGWINDOW:4:20 MINLEN:25

# align to reference genome
mkdir -p data/process/cdiff_16S/cdiff_genome_alignment
bwa index data/reference/cdifficile630_reference_genome.fasta
bwa mem data/reference/cdifficile630_reference_genome.fasta data/process/cdiff_16S/SRR5804616_1.trim.fastq.gz data/process/cdiff_16S/SRR5804616_2.trim.fastq.gz > data/process/cdiff_16S/cdiff_genome_alignment/SRR804616.aligned.sam
samtools view -S -b data/process/cdiff_16S/cdiff_genome_alignment/SRR804616.aligned.sam > data/process/cdiff_16S/cdiff_genome_alignment/SRR804616.aligned.bam

# get consensus sequence
bcftools mpileup -f data/reference/cdifficile630_reference_genome.fasta data/process/cdiff_16S/cdiff_genome_alignment/SRR804616.aligned.bam | bcftools call -mv -Oz -o data/process/cdiff_16S/cdiff_genome_alignment/calls.vcf.gz
tabix data/process/cdiff_16S/cdiff_genome_alignment/calls.vcf.gz
cat data/reference/cdifficile630_reference_genome.fasta | bcftools consensus data/process/cdiff_16S/cdiff_genome_alignment/calls.vcf.gz > data/process/cdiff_16S/consensus_cdifficile_431.fasta
