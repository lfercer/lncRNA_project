#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:05:00
#SBACTH --error=/data/users/lfernandez/rnaseq_course/assembly/errors/error_count_transcripts_%j.e
#SBATCH --output=/data/users/lfernandez/rnaseq_course/assembly/output_count_transcripts_%j.out
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/assembly/output_files


cd $INPUT_DIR

#Count the number of transcripts
n_transcripts=$( awk '$3 == "transcript" {count++} END {print count}' meta_assembly.gtf )
echo "Number of transcripts: $n_transcripts" > $OUTPUT_DIR/count_assembly.txt

#Count the number of exons
n_exons=$( awk '$3 == "exon" {count++} END {print count}' meta_assembly.gtf ) 
echo "Number of exons: $n_exons" >> $OUTPUT_DIR/count_assembly.txt

#Count the number of genes
n_genes=$( awk 'NR>1 {print $10}' meta_assembly.gtf |  sort | uniq | wc -l )
echo "Number of genes: $n_genes" >> $OUTPUT_DIR/count_assembly.txt

#Count novel transcripts
n_novel_transcripts=$(awk ' $3 == "transcript" {if($12 !~ "ENST") print $12}' meta_assembly.gtf | sort | uniq |  wc -l )
echo "Number of novel transcripts: $n_novel_transcripts" >> $OUTPUT_DIR/count_assembly.txt

#Count novel genes
n_novel_genes=$( awk '{if($10 !~ "ENSG") print $10}' meta_assembly.gtf |  sort | uniq | wc -l )
echo "Number of novel genes: $n_novel_genes " >> $OUTPUT_DIR/count_assembly.txt

#Count single exon transcripts
n_single_exon_transcripts=$( awk '$3 == "exon" {print $12}' meta_assembly.gtf | sort | uniq -c | awk '$1 == 1' | wc -l )
echo "Number of single exon transcripts: $n_single_exon_transcripts" >> $OUTPUT_DIR/count_assembly.txt

#Count single exon genes
n_single_exon_genes=$( awk '$3 == "exon" {print $10}' meta_assembly.gtf | sort | uniq -c | awk '$1 == 1' | wc -l )
echo "Number of singe exon genes: $n_single_exon_genes" >> $OUTPUT_DIR/count_assembly.txt