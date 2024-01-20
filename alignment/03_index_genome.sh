#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=03_index_genome
#SBATCH --error=/data/users/lfernandez/rnaseq_course/alignment/errors/error_%j.e
#SBATCH --partition=pall

INPUT_DIR=/data/courses/rnaseq_course/lncRNAs/Project2/references
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/alignment/output_files

module load UHTS/Aligner/hisat/2.2.1

cd $INPUT_DIR

hisat2-build GRCh38.genome.fa $OUTPUT_DIR/GRCh38_index