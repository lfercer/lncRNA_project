#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --job-name=14_run_intersect
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_intersect_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_intersect_%j.out
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files

module load UHTS/Analysis/BEDTools/2.29.2;

cd $INPUT_DIR

#Strand and windows are not needed in this step
bedtools intersect -v -a novel_transcripts.bed -b annotated_transcripts.bed > $OUTPUT_DIR/novel_intergenic.bed

#Count of novel intergenic transcripts identified
#awk '{++n} END{print n}' $OUTPUT_DIR/novel_intergenic.bed > $OUTPUT_DIR/step6_answers.txt