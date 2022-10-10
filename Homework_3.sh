#!/bin/bash
#SBATCH --job-name=denovo-assembly		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/tc88074/log.%j			    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=tc88074@uga.edu                   # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory variable
OUTDIR="/work/gene8940/tc88074/homework_3"             # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi


#transfer e.coli str using curl
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ecoli_MG1655.fna

#load modules
module load canu/1.9-GCCcore-8.3.0-Java-11
module load SPAdes/3.14.1-GCC-8.3.0-Python-3.7.4
module load QUAST/5.0.2-foss-2019b-Python-3.7.4
module load MUMmer/3.23_conda
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

#run jobs
canu -p ecoli -d $OUTDIR/canu genomeSize=4.8m useGrid=false -pacbio-raw /work/gene8940/instructor_data/ecoli_p6_25x.filtered.fastq.gz

spades.py -t 6 -k 21,33,55,77 --isolate --memory 24 --pe1-1 /work/gene8940/instructor_data/s_6_1.fastq.gz --pe1-2 /work/gene8940/instructor_data/s_6_2.fastq.gz -o $OUTDIR/spades

quast.py -o quast -t 6 -r $OUTDIR/ecoli_MG1655.fna $OUTDIR/canu/ecoli.contigs.fasta $OUTDIR/spades/scaffolds.fasta

mkdir mummer

nucmer $OUTDIR/ecoli_MG1655.fna $OUTDIR/canu/ecoli.contigs.fasta -p mummer/ecoli.canu
delta-filter -1 mummer/ecoli.canu.delta > mummer/ecoli.canu.1delta
show-coords mummer/ecoli.canu.1delta
mummerplot --size large -layout --color -f --png mummer/ecoli.canu.1delta -p mummer/ecoli.canu


nucmer $OUTDIR/ecoli_MG1655.fna $OUTDIR/spades/scaffolds.fasta -p mummer/ecoli.spades
delta-filter -1 mummer/ecoli.spades.delta > mummer/ecoli.spades.1delta
show-coords mummer/ecoli.spades.1delta
mummerplot --size large -layout --color -f --png mummer/ecoli.spades.1delta -p mummer/ecoli.spades
