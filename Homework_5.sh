#!/bin/bash
#SBATCH --job-name=RNA_seq		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/tc88074/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=tc88074@uga.edu                   # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/tc88074/homework_5"             # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#set output kallisto directory variable
KALLISTODIR="/work/gene8940/tc88074/homework_5/kallisto"

#if kallisto directory doesn't exist, create it
if [ ! -d $KALLISTODIR ]
then
  mkdir -p $KALLISTODIR
fi

cd $OUTDIR

#download reference genome refseq CDS fasta file
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2 | gunzip -c > $OUTDIR/ecoli_MG1655_refseq_cds.fa

#load kallisto module
module load kallisto/0.46.1-foss-2019b

#make kallisto index for CDS refseq fasta file
kallisto index -i $OUTDIR/ecoli_MG1655_refseq_cds.fa.idx $OUTDIR/ecoli_MG1655_refseq_cds.fa

#run kallisto quant on all 4 samples
for i in SRR5344681 SRR5344682 SRR5344683 SRR5344684
do
  kallisto quant -t 6 $THREADS -b 100 -i $OUTDIR/ecoli_MG1655_refseq_cds.fa.idx -o $KALLISTODIR/$i /work/gene8940/instructor_data/${i}_1.fastq.gz /work/gene8940/instructor_data/${i}_2.fastq.gz
done

source activate R
R --no-save < $Tasneem/Desktop/GitHub-BINF8940/BINF8940/homework5.r
