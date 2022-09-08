#!/bin/bash
#SBATCH --job-name=BLAST-test		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=2		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=2:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/cbergman/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=cbergman@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/cbergman/blast-test"             # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load blast module
module load BLAST+/2.9.0-gompi-2019b

# run blast against local copy of NCBI nucleotide database
blastn -num_threads 2 -query /usr/local/training/sample.fasta -db /db/ncbiblast/nt/06042020/nt -out $OUTDIR/sample.fa.blastn.${SLURM_JOB_ID}.tsv -outfmt 6 -max_target_seqs 2

# same command as above, but split over multiple lines for improved readability
# blastn -num_threads 2 -query /usr/local/training/sample.fasta \
#        -db /db/ncbiblast/nt/06042020/nt \
#        -out $OUTDIR/sample.fa.blastn.${SLURM_JOB_ID}.tsv \
#        -outfmt 6 \
#        -max_target_seqs 2
