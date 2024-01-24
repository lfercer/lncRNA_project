#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --job-name=17_intersect_3prime
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_intersect3_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_intersect3_%j.out
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
REF_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/references

module load UHTS/Analysis/BEDTools/2.29.2;

cd $INPUT_DIR

#same strand (+) (overlaps be on the same strand)
bedtools intersect -wa -s -a novel3window_corrected.bed -b $REF_DIR/atlas.clusters.2.0.GRCh38.96.bed > $OUTPUT_DIR/overlap3prime.bed

