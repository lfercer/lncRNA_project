#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:05:00
#SBATCH --error=/data/users/lfernandez/rnaseq_course/quantification/errors/error_quantification_%j.e
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/quantification/output_files

cd $INPUT_DIR
SAMPLES_DIR=$(ls | grep -E "^[3P]")

for sample in $SAMPLES_DIR
do
    cd "$INPUT_DIR/$sample"

    # Total of estimated counts
    est_c=$(awk -F'\t' 'BEGIN{n=0} NR>1 {n=n+$4} END{print n}' abundance.tsv)
    echo "total of estimated counts for $sample: $est_c" >> $INPUT_DIR/step4_questions.txt

    # Number of transcripts detected
    num_transcripts=$(awk -F'\t' 'NR>1 {if($4 > 0) ++n} END{print n}' abundance.tsv)
    echo "number of transcripts detected for $sample: $num_transcripts" >> $INPUT_DIR/step4_questions.txt

    # Number of genes detected
    num_genes=$(awk -F'\t' '$1 ~ /\.1$/ {if($4 > 0) ++n } END{print n}' abundance.tsv)
    echo "number of genes detected for $sample: $num_genes" >> $INPUT_DIR/step4_questions.txt

    #Number of novel transcripts detected
    novel_transcripts=$(awk -F'\t' 'NR>1 && $1 !~ /'ENS'/ {if($4 > 0) ++n} END{print n}' abundance.tsv) 
    echo "number of novel transcripts detected for $sample: $novel_transcripts" >> $INPUT_DIR/step4_questions.txt

    #Number of novel genes detected
    novel_genes=$(awk -F'\t' '$1 ~ /\.1$/ && $1 !~ /'ENS'/ {if($4 > 0) ++n } END{ print n}' abundance.tsv)
    echo "number of novel genes detected for $sample: $novel_genes" >> $INPUT_DIR/step4_questions.txt
done




