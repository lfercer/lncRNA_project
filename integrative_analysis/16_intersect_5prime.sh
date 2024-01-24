#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --job-name=16_run_5prime
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_intersect5_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_intersect5_%j.out
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
REF_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/references

module load UHTS/Analysis/BEDTools/2.29.2;

cd $INPUT_DIR

#same strand (+) (overlaps be on the same strand)
bedtools intersect -wa -s -a novel5window.bed -b $REF_DIR/refTSS_v4.1_human_coordinate.hg38.bed > $OUTPUT_DIR/overlap5prime.bed

