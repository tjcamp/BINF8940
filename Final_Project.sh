#!/bin/bash
#SBATCH --job-name=CYP2D6_final_project		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/tc88074/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=tc88074@uga.edu                   # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/tc88074/final_project"             # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi



cd $OUTDIR

#extract paired end illumina reads from 1000 Genomes Project phase 3: 30X coverage whole genome sequencing (Accession: PRJEB31736)
module load SRA-Toolkit/2.11.1-centos_linux64
prefetch -O $OUTDIR ERR4048409
fastq-dump --split-files --gzip $OUTDIR/ERR4048409 -O $OUTDIR

#extract paired end illumina reads from 1000 Genomes Project phase 3: 30X coverage whole genome sequencing of 698 samples to complete trios (Accession: PRJEB36890)
#prefetch -O $OUTDIR ERR3989458
#fastq-dump --split-files --gzip $OUTDIR/ERR3989458 -O $OUTDIR

#download reference genome GCF_000001405.40_GRCh38.p14
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz | gunzip -c > GCF_000001405.fna

#construct BWA index and map the reads to reference genome and generate samtools index
module load BWA/0.7.17-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0

bwa index GCF_000001405.fna
bwa mem -t 6 GCF_000001405.fna ERR4048409_1.fastq.gz ERR4048409_2.fastq.gz | samtools view - -O BAM | samtools sort --threads 6 > ERR4048409.sorted.bam
samtools index -@ 6 ERR4048409.sorted.bam

#variant mapping with bcftools
module load BCFtools/1.10.2-GCC-8.3.0

bcftools mpileup -Oz --threads 6 --min-MQ 60 -f GCF_000001405.fna ERR4048409.sorted.bam > ERR4048409.sorted.mpileup.vcf.gz
bcftools call -Oz -m -v --threads 6 --ploidy 2 ERR4048409.sorted.mpileup.vcf.gz > ERR4048409.sorted.mpileup.call.vcf.gz
bcftools filter -Oz -e 'QUAL<40 || DP<10' ERR4048409.sorted.mpileup.call.vcf.gz > ERR4048409.sorted.mpileup.call.filter.vcf.gz
bcftools view -H -v snps ERR4048409.sorted.mpileup.call.filter.vcf.gz | wc -l > results.snps.txt
bcftools view -H -v indels ERR4048409.sorted.mpileup.call.filter.vcf.gz | wc -l > results.indels.txt

#create IGV readable index file
bcftools index ERR4048409.sorted.mpileup.call.filter.vcf.gz
