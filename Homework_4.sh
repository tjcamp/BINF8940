#!/bin/bash
#SBATCH --job-name=variant_calling		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/tc88074/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=tc88074@uga.edu                   # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/tc88074/homework_4"             # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi



cd $OUTDIR

#extract paired end illumina reads
module load SRA-Toolkit/2.11.1-centos_linux64
vdb-config --interactive
prefetch -O $OUTDIR SRR8082143
fastq-dump --split-files --gzip $OUTDIR/SRR8082143 -O $OUTDIR

#download reference genome
curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > GCF_000005845.fna

#construct BWA index and map the reads to reference genome and generate samtools index
module load BWA/0.7.17-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0

bwa index GCF_000005845.fna
bwa mem -t 6 GCF_000005845.fna SRR8082143_1.fastq.gz SRR8082143_2.fastq.gz | samtools view - -O BAM | samtools sort --threads 6 > SRR8082143.sorted.bam
samtools index -@ 6 SRR8082143.sorted.bam

#variant mapping
module load BCFtools/1.10.2-GCC-8.3.0

bcftools mpileup -Oz --threads 6 --min-MQ 60 -r NC_000913.3:2173000-2174000 -f GCF_000005845.fna SRR8082143.sorted.bam > SRR8082143.sorted.mpileup.vcf.gz
bcftools call -Oz -m -v --threads 6 --ploidy 1 SRR8082143.sorted.mpileup.vcf.gz > SRR8082143.sorted.mpileup.call.vcf.gz
bcftools filter -Oz -e 'QUAL<40 || DP<10' SRR8082143.sorted.mpileup.call.vcf.gz > SRR8082143.sorted.mpileup.call.filter.vcf.gz

#create IGV readable index file
bcftools index SRR8082143.sorted.mpileup.call.filter.vcf.gz
