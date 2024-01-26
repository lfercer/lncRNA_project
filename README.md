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
  - 09_transcriptome_fasta.sh : create a transcriptome fasta file from meta-assembly.gtf and the kallisto indexed file
  - 10_kallisto_quantification.sh : run Kallisto per each sample
  - 11_summary_quantification.sh : to answer the questions of step 4.
- step4_questions.txt :  answers to questions of step 4. Total of estimated counts, and number of transcripts, genes, novel transcripts and novel genes detected per sample.
- transcript_expression_table.txt : deliverable file with transcript_id as row names and samples as column names. The expression levels are expressed in estimated counts.
- In the directory (/data/users/lfernandez/rnaseq_course/quantification/output_files) of the cluster you can find the transcriptome fasta files as well as the kallisto indexed file, and one directory per sample with the results from kallisto quantification.

STEP 5: differential_expression
- 08_table_mapping_transcriptsID.sh : bash script to create my_annotation_table.txt from the transcriptome assembly "meta_assembly.gtf".
- 12_sleuth.R : R script to perform the differential expression analysis and subsequents plots. It requires the files: samples.txt and my_annotation_table.txt
  
- /required_files: 
  - samples.txt : file that specifies the sample, condition and path to kallisto quantification output. It is required for the preparation step ('sleuth_prep') in Sleuth
  - my_annotation_table.txt : it compiles transcript_id, gene_id, gene_name and biotype information for each transcript present in the transcriptome assembly (meta_assembly.gtf). It is necessary to map each transcript to one gene due to the differential expression analysis is performed at transcript level.
  
- differential_expression_table.txt : deliverable file. It contains transcript_id, gene_id, gene_name, biotype, log2FC and q-value.

  STEP 6: integrative_analysis
- Bash scripts:
  - 13_create_bed.sh : to create two BED files (novel_transcripts.bed and annotated_transcripts.bed, they can be found in the cluster: /data/users/lfernandez/rnaseq_course/integrative_analysis/output_files ) from the meta_assembly.gtf
  - 14_intersect_intergenic.sh : to find the novel intergenic transcripts.
  - 15_bed_window.sh :  to create a BED files with a window of 100 nucleotides around the transcription start sites (TSS) and transcription end sites (TES) of the novel transcripts. Required for the next step.
  - 16_intersect_5prime.sh : to find the overlaps of the TSS of the novel transcripts and the FANTOM CAGE of TSS.
  -  17_intersect_3prime.sh :to find the overlaps of the TES of the novel transcripts and the FANTOM CAGE of TES.
  -  18_run_CPAT.sh : to estimate the protein coding potential of the novel transcripts.
  -  19_step6_questions.sh :  to answer the questions of the step 6.
- 20_prioritization.R : R script to join information from DE and integrative analysis and make a sorted table of the transcripts. It also includes the code to generate the plots of the integrative analysis.

- /output_files :
  - step6_answers.txt : deliverable. Statistics and porcentages answering the questions of the step
  - novel_intergenic_transcripts_plot.pdf : deliverable. It shows the percentage of intergenic, 5’ well-annotated, 3’ well-annotated, 5’ and 3’ well-annotated and protein coding transcripts above the total of novel transcripts.
  - novel_transcripts_plot.pdf : deliverable. It shows the percentage of intergenic, 5’ well-annotated, 3’ well-annotated, 5’ and 3’ well-annotated and protein coding transcripts above the total of novel intergenic transcripts.
  - prioritization_novel_table.tsv : list of novel transcripts sorted regarding the prioritization criteria.

 - All the BED files required can be found in the cluster: /data/users/lfernandez/rnaseq_course/integrative_analysis/output_files

 - All the reference files used can be found: /data/courses/rnaseq_course/lncRNAs/Project1/references

- /Supplementary : contains Supplementary materials mentionated in the scientific report.



    
