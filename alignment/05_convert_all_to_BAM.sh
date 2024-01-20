s#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --job-name=04_run_HISAT2
#SBATCH --error=/data/users/lfernandez/rnaseq_course/alignment/errors/error_convertBAM_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/alignment/job_convertBAM_%j.out
#SBATCH --partition=pall
#SBATCH --array=0-5

INPUT_DIR=/data/users/lfernandez/rnaseq_course/alignment/output_files


cd $INPUT_DIR

module load UHTS/Analysis/samtools/1.10

SAM_files=( $( ls | grep -E '\.sam$') )
current_file=${SAM_files[$SLURM_ARRAY_TASK_ID]}
name=${current_file:0:3}

samtools view -Sb $current_file -o $name.bam
samtools sort -O bam $name.bam -o sorted_$name.bam
samtools index sorted_$name.bam