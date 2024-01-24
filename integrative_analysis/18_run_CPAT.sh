#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --job-name=19_CPAT
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_cpat_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_cpat_%j.out
#SBATCH --partition=pall

REFERENCES_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/references
NOVEL_BED=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files/novel_transcripts.bed
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files



module load UHTS/Analysis/BEDTools/2.29.2;
module load SequenceAnalysis/GenePrediction/cpat/1.2.4;


#get fasta file for novel transcripts
# -s force strandedness
bedtools getfasta -s -name -fi $REFERENCES_DIR/GRCh38.genome.fa -bed $NOVEL_BED -fo $OUTPUT_DIR/novel_transcripts.fa

#run cpat
cpat.py -x $REFERENCES_DIR/Human_Hexamer.tsv -d $REFERENCES_DIR/Human_logitModel.RData -g $OUTPUT_DIR/novel_transcripts.fa -o $OUTPUT_DIR/novel_coding_potential
