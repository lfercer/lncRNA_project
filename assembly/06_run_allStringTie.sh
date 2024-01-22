#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=06_run_StringTie
#SBATCH --error=/data/users/lfernandez/rnaseq_course/assembly/errors/error_StringTie_%j.e
#SBATCH --output=/data/users/lfernandez/rnaseq_course/assembly/output_StringTie_%j.out
#SBATCH --partition=pall
#SBATCH --array=0-5


INPUT_DIR=/data/users/lfernandez/rnaseq_course/alignment/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
ANNOTATION_DIR=/data/courses/rnaseq_course/lncRNAs/Project2/references

cd $INPUT_DIR

SAM_FILES=( $( ls | grep '^sorted.*bam$' ))
current_file=${SAM_FILES[$SLURM_ARRAY_TASK_ID]}
name=${current_file:7:3}


module load UHTS/Aligner/stringtie

stringtie $current_file -G $ANNOTATION_DIR/gencode.v44.annotation.gtf --rf -o $OUTPUT_DIR/assembly_$name.gtf

