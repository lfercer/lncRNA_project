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

STEP 3: assembly
- Bash scripts: 06_run_allStringTie.sh, 06_run_merge_StringTie.sh and 07_count_transcripts.sh to answer the questions of step 3.
- /output_files: count_assembly.txt -> answers to the questions of step 3.
- meta-assembly GTF file can be found in the cluster: /data/users/lfernandez/rnaseq_course/assembly/output_files/meta_assembly.gtf
        
STEP 4: quantification
- Bash scripts: 
  09_transcriptome_fasta.sh : create a transcriptome fasta file from meta-assembly.gtf and the kallisto indexed file
  10_kallisto_quantification.sh : run Kallisto per each sample
  11_summary_quantification.sh : to answer the questions of step 4.
- step4_questions.txt :  answers to questions of step 4. Total of estimated counts, and number of transcripts, genes, novel transcripts and novel genes detected per sample.
- transcript_expression_table.txt : deliverable file with transcript_id as row names and samples as column names. The expression levels are expressed in estimated counts.
- In the directory (/data/users/lfernandez/rnaseq_course/quantification/output_files) of the cluster you can find the transcriptome fasta files as well as the kallisto indexed file, and one directory per sample with the results from kallisto quantification.

STEP 5: differential_expression
- 08_table_mapping_transcriptsID.sh : bash script to create my_annotation_table.txt from the transcriptome assembly "meta_assembly.gtf".
- 13_sleuth.R : R script to perform the differential expression analysis and subsequents plots. It requires the files: samples.txt and my_annotation_table.txt
  
- /required_files: 
-   samples.txt : file that specifies the sample, condition and path to kallisto quantification output. It is required for the preparation step ('sleuth_prep') in Sleuth
-   my_annotation_table.txt : it compiles transcript_id, gene_id, gene_name and biotype information for each transcript present in the transcriptome assembly (meta_assembly.gtf). It is necessary to map each transcript to one gene due to the differential expression analysis is performed at transcript level.
  
- differential_expression_table.txt : deliverable file. It contains transcript_id, gene_id, gene_name, biotype, log2FC and q-value.
