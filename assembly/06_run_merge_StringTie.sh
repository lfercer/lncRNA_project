#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBACTH --error=/data/users/lfernandez/rnaseq_course/assembly/errors/error_merge_StringTie_%j.e
#SBATCH --output=/data/users/lfernandez/rnaseq_course/assembly/output_merge_StringTie_%j.out
#SBATCH --partition=pall



INPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
ANNOTATION_DIR=/data/courses/rnaseq_course/lncRNAs/Project2/references

cd $INPUT_DIR

GTF_FILES=$(ls | grep '\.gtf')




module load UHTS/Aligner/stringtie
stringtie --merge $GTF_FILES -G $ANNOTATION_DIR/gencode.v44.annotation.gtf --rf -o $OUTPUT_DIR/meta_assembly.gtf