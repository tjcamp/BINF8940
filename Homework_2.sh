#!/bin/bash
#SBATCH --job-name=GC-content		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=2		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=2:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/tc88074/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=tc88074@uga.edu                   # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/tc88074/homework_2"             # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi


#transfer e.coli str using curl
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz | gunzip -c > $OUTDIR/ecoli_MG1655.gff
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ecoli_MG1655.fna

#load modules
module load BEDOPS/2.4.39-foss-2019b
module load BEDTools/2.29.2-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load ucsc/359

#convert gff file to bed format
convert2bed --input=gff < $OUTDIR/ecoli_MG1655.gff > $OUTDIR/ecoli_MG1655.bed
grep "CDS" $OUTDIR/ecoli_MG1655.bed > $OUTDIR/ecoli_MG1655_cds.bed
samtools faidx $OUTDIR/ecoli_MG1655.fna
cut -f1,2 $OUTDIR/ecoli_MG1655.fna.fai > $OUTDIR/ecoli_MG1655.genome.txt
bedtools complement -i $OUTDIR/ecoli_MG1655_cds.bed -g $OUTDIR/ecoli_MG1655.genome.txt > $OUTDIR/ecoli_MG1655_intergenic.bed
bedtools getfasta -fi $OUTDIR/ecoli_MG1655.fna -bed $OUTDIR/ecoli_MG1655_cds.bed -fo $OUTDIR/ecoli_MG1655_cds.fna
bedtools getfasta -fi $OUTDIR/ecoli_MG1655.fna -bed $OUTDIR/ecoli_MG1655_intergenic.bed -fo $OUTDIR/ecoli_MG1655_noncds.fna

faCount $OUTDIR/ecoli_MG1655.cds.fna -summary > $OUTDIR/results_cds.txt
faCount $OUTDIR/ecoli_MG1655.noncds.fna -summary > $OUTDIR/results_noncds.txt
