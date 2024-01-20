#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=200MB
#SBATCH --time=00:02:00
#SBATCH --job-name=01_count_reads
#SBATCH --error=/data/users/lfernandez/rnaseq_course/quality/errors/error_%j.e
#SBATCH --output=/data/users/lfernandez/rnaseq_course/quality/jobs/output_%j.o
#SBATCH --partition=pall
#SBATCH --array=0-11

READS_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/quality/output_files

#change to raw data directory
cd $READS_DIR

#create the array ("list") of all my samples (paraclonal and parental)
INPUT_FILES=( $(ls | grep -E '^[3P].+\.fastq\.gz') )


# Sort the array in ascending order
IFS=$'\n' sorted_array=($(printf "%s\n" "${INPUT_FILES[@]}" | sort))
unset IFS

# Access the sorted array elements
current_file=${sorted_array[$SLURM_ARRAY_TASK_ID]}


count=$(zcat $current_file | awk 'BEGIN {n=0} /^@/ {++n} END {print n}')
echo "$current_file $count" >> "$OUTPUT_DIR/count_reads.txt" 

