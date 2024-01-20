#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=200MB
#SBATCH --time=00:05:00
#SBATCH --job-name=02_run_fastqc
#SBATCH --error=/data/users/lfernandez/rnaseq_course/quality/errors/error_%j.e
#SBATCH --output=/data/users/lfernandez/rnaseq_course/quality/jobs/slurm_%j.o
#SBATCH --partition=pall
#SBATCH --array=0-11


READS_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/quality/output_files

#load package fastqc
module load UHTS/Quality_control/fastqc/0.11.9

#change to raw data directory
cd $READS_DIR

INPUT_FILES=( $(ls | grep ^[3P]) )
current_file=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}


fastqc -o $OUTPUT_DIR $current_file





