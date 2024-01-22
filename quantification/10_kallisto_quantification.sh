#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --error=/data/users/lfernandez/rnaseq_course/quantification/errors/error_kallisto_quantification_%j.e
#SBATCH --partition=pall
#SBATCH --array=0-5

module load UHTS/Analysis/kallisto/0.46.0

FASTQ_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
TRANSCRIPTOME_DIR=/data/users/lfernandez/rnaseq_course/quantification/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/quantification/output_files

cd $FASTQ_DIR

FILES_mate1=( $( ls | grep -E '^[3P].+R1' ) )
FILES_mate2=( $( ls | grep -E '^[3P].+R2' ) )

current_file=${FILES_mate2[$SLURM_ARRAY_TASK_ID]}
name=${current_file:0:3}

#create one directory to store results per each sample
mkdir -p /data/users/lfernandez/rnaseq_course/quantification/output_files/$name


kallisto quant -i $TRANSCRIPTOME_DIR/kallisto_transcriptome_index -o $OUTPUT_DIR/$name  --rf-stranded -b 100  ${FILES_mate1[$SLURM_ARRAY_TASK_ID]} ${FILES_mate2[$SLURM_ARRAY_TASK_ID]}