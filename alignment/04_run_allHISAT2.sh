#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=04_run_HISAT2
#SBATCH --error=/data/users/lfernandez/rnaseq_course/alignment/errors/error_HISAT_%j.e
#SBATCH --partition=pall
#SBATCH --array=0-5

INPUT_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
GENOME_DIR=/data/users/lfernandez/rnaseq_course/alignment/genome_index
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/alignment/output_files


module load UHTS/Aligner/hisat/2.2.1
cd $INPUT_DIR

FILES_mate1=( $( ls | grep -E '^[3P].+R1' ) )
FILES_mate2=( $( ls | grep -E '^[3P].+R2' ) )
current_file=${FILES_mate2[$SLURM_ARRAY_TASK_ID]}
name=${current_file:0:3}

hisat2 -x $GENOME_DIR/GRCh38_index -1 $FILES_mate1 -2 $FILES_mate2  -S $OUTPUT_DIR/$name.sam --rna-strandness RF