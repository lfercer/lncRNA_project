#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=03:00:00
#SBATCH --job-name=18_count_5_3_overlaps
#SBATCH --error=/data/users/lfernandez/rnaseq_course/integrative_analysis/errors/error_count_%j.e
#SBACTH --output=/data/users/lfernandez/rnaseq_course/integrative_analysis/job_count_%j.out
#SBATCH --partition=pall

INPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files
OUTPUT_DIR=/data/users/lfernandez/rnaseq_course/integrative_analysis/output_files


cd $INPUT_DIR

#Question 1
#take transcripts_id from each bed file and sort them
awk '{print $4}' overlap5prime.bed | sort > $OUTPUT_DIR/transcripts5prime_sorted.txt
awk '{print $4}' overlap3prime.bed | sort > $OUTPUT_DIR/transcripts3prime_sorted.txt


# join overlaps5prime and overlaps3prime
join -1 1 -2 1 transcripts5prime_sorted.txt transcripts3prime_sorted.txt > $OUTPUT_DIR/all_overlaps5and3.txt


#count your total of transcripts --> lines in novel_transcripts.bed
total=$(awk '{++n} END{print n}' novel_transcripts.bed)


#count lines of all_overlaps5and3.txt  
transcripts5and3=$(awk '{++n} END{print n}' all_overlaps5and3.txt)
transcripts5=$(awk '{++n} END{print n}' overlap5prime.bed)
transcripts3=$(awk '{++n} END{print n}' overlap3prime.bed)



#Question 3
#Count of novel intergenic transcripts identified
novel_intergenic=$(awk '{++n} END{print n}' novel_intergenic.bed)
result5=$(echo "scale=3 ; ($novel_intergenic * 100 / $total ) " | bc )
echo "Total of novel transcripts: $total" > $OUTPUT_DIR/step6_answers.txt
echo "Number of novel intergenic transcripts: $novel_intergenic" >> $OUTPUT_DIR/step6_answers.txt
echo "Percentage of novel transcripts considered as intergenic: $result5%" >> $OUTPUT_DIR/step6_answers.txt


#Question 1 answers:
#calculate percentage
result=$(echo "scale=3 ; ($transcripts5and3 * 100 / $total ) " | bc )
echo "Percentage of novel transcripts with 5' and 3' good annotations: $result%" >> $OUTPUT_DIR/step6_answers.txt

result2=$(echo "scale=3 ; ($transcripts5 * 100 / $total ) " | bc )
echo "Percentage of novel transcripts with 5' good annotations: $result2%" >> $OUTPUT_DIR/step6_answers.txt

result3=$(echo "scale=3 ; ($transcripts3 * 100 / $total ) " | bc )
echo "Percentage of novel transcripts with  3' good annotations: $result3%" >> $OUTPUT_DIR/step6_answers.txt

#Question 2
# Cutoff chosen for human: 0.364 (optimum determined from TG-ROC)
#Transcripts with protein coding probability > 0.364 are considerated protein coding
prot_coding=$(awk '{if($5 > 0.364) ++n} END{print n}' novel_coding_potential )
result4=$(echo "scale=3 ; ($prot_coding * 100 / $total ) " | bc )
echo "Using a threshold of 0.364 to consider as protein coding transcript:"
echo "Percentage of novel transcripts considered as protein coding: $result4%" >> $OUTPUT_DIR/step6_answers.txt
