#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --error=/data/users/lfernandez/rnaseq_course/quantification/errors/error_transcriptome_fasta_%j.e
#SBATCH --partition=pall

#Create a directory if it does not exist
mkdir -p /data/users/lfernandez/rnaseq_course/quantification/output_files
mkdir -p /data/users/lfernandez/rnaseq_course/quantification/errors

GENOME_FILE=/data/courses/rnaseq_course/lncRNAs/Project2/references/GRCh38.genome.fa
INPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/quantification/output_files


module load UHTS/Assembler/cufflinks/2.2.1
module load UHTS/Analysis/kallisto/0.46.0


#Create a transcriptome fasta file
#  -w     write a fasta file with spliced exons for each GTF transcript
gffread -w $OUTPUT_DIR/transcriptome.fasta -g $GENOME_FILE $INPUT_DIR/meta_assembly.gtf

#Create a kallisto index
kallisto index -i $OUTPUT_DIR/kallisto_transcriptome_index $OUTPUT_DIR/transcriptome.fasta

