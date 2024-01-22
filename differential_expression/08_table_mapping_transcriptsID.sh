#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:05:00
#SBATCH --error=/data/users/lfernandez/rnaseq_course/assembly/errors/error_table_transcriptsID_%j.e
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
INPUT_FILE=/data/users/lfernandez/rnaseq_course/assembly/output_files/meta_assembly.gtf
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
GENCODE_DIR=/data/courses/rnaseq_course/lncRNAs/Project2/references

#Create a file with a table with 2 columns: gene_id transcript_id  from meta_assembly.gtf
#Transcripts and genes in "GL00.." chromosome patches are excluded.
#It only includes transcripts in "ch".


cd $INPUT_DIR

echo 'transcript_id gene_id gene_name' > $OUTPUT_DIR/transcripts_mapped.txt

#create table from metaassembly.gtf
#gsub is used to remove the ";"
awk '$1 ~ /ch/ && $3 == "transcript" {gsub(/;/, "", $12); gsub(/;/, "", $10); gsub(/;/, "", $14) ; print $12, $10, $14}' $INPUT_FILE | sort  >> $OUTPUT_DIR/transcripts_mapped.txt



#create table from gencode.gtf
#transcript_id gene_type
echo 'transcript_id biotype' > $OUTPUT_DIR/biotype_table.txt
awk '$3 == "transcript" {gsub(/;/, "", $12) ; gsub(/;/, "", $14) ; print $12,$14}' $GENCODE_DIR/gencode.v44.annotation.gtf | sort >> $OUTPUT_DIR/biotype_table.txt


#join both tables (add gene_type as 4rd column from reference genome GTF)
join --nocheck-order -a 1 -a 2 -1 1 -2 1 transcripts_mapped.txt biotype_table.txt > my_annotation_table.txt









