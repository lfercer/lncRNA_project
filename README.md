# lncRNA_project
# RNA sequencing course
# Name: Laura Fernández Cerro

Each directory refers to one step of the project. The description of the contents of each directory is indicated afterwards.
Bash scripts are enumerated by the order they were submitted.

STEP 1: quality
- Bash scripts
- count_reads.txt : count of the reads for each sample. Output file from 01_count_allreads_2.sh
- /output_fastqc : it contains htmls for each sample and mate to fastqc output.
  
STEP 2: alignment
- Bash scripts
- /output_files : summary of statistics and alignment rates of each sample
- SAM and BAM files can be found in the following directory in the cluster: /data/users/lfernandez/rnaseq_course/alignment/output_files
