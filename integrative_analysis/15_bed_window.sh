#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=14_run_intersect
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_bedwindow_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_bedwindow_%j.out
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files


cd $INPUT_DIR

#create a 5´ bed file of novel transcripts
#create a window of 100 nt
#2nd field for + and 3rd field for -
awk '{if($6=="+") print $1"\t"($2-50)"\t"($2+50)"\t"$4"\t"$5"\t"$6 ; else print $1"\t"($3-50)"\t"($3+50)"\t"$4"\t"$5"\t"$6}' novel_transcripts.bed > novel5window.bed

#create a 3´ bed file of novel transcripts
#create a window of 100 nt
#3rd field for + and 2nd field for -
awk '{if($6=="+") print $1"\t"($3-50)"\t"($3+50)"\t"$4"\t"$5"\t"$6; else print $1"\t"($2-50)"\t"($2+50)"\t"$4"\t"$5"\t"$6}' novel_transcripts.bed > novel3window.bed
#correction of negative number
awk '{if($2 < 0 ) print $1"\t"0"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' novel3window.bed > novel3window_corrected.bed



#check after if there are negative numbers after creating the window
#grep -e '-[123456789]' novel5window.bed --> no negative numbers
#grep -e '-[123456789]' novel3window.bed --> only one
