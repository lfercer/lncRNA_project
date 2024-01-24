#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=13_create_bed
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_BEDfile_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_BEDfile_%j.out
#SBATCH --partition=pall

MERGE_GFT=/data/users/lfernandez/rnaseq_course/assembly/output_files/meta_assembly.gtf
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files

# Exclude transcripts that are in chromosome patches. Only include the one that are in the chromosomes.

#Create bed file with novel transcripts

awk '$3 == "transcript" && $1 ~ /ch/ && $12 !~ /ENST/ {gsub(/;/,"",$12)  ; print $1"\t"$4"\t"$5"\t"$12"\t"$6"\t"$7}' $MERGE_GFT > $OUTPUT_DIR/novel_transcripts.bed


#Create bed file with annotated transcripts
awk '$3 == "transcript" && $1 ~ /ch/ && $12 ~ /ENST/ {gsub(/;/,"",$12)  ; print $1"\t"$4"\t"$5"\t"$12"\t"$6"\t"$7}' $MERGE_GFT > $OUTPUT_DIR/annotated_transcripts.bed

#Create bed file with all transcripts 
#awk '$3 == "transcript" && $1 ~ /ch/ {gsub(/;/,"",$12)  ; print $1"\t"$4"\t"$5"\t"$12"\t"$6"\t"$7}' $MERGE_GFT > $OUTPUT_DIR/all_transcripts.bed