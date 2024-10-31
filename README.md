# Annotation and characterization of lncRNAs in Lung Cancer
## RNA sequencing course - University of Bern

## **Project Description**

This project is focused on the annotation and characterization of long non-coding RNAs (lncRNAs) in lung cancer, specifically non-small cell lung carcinoma (NSCLC), using RNA sequencing data from two phenotypically distinct A549 cell populations.Through differential expression and integrative analysis, prioritized lncRNA targets were identified based on intergenic nature, transcription site annotation quality, and protein-coding potential.

Full report can be found: [report](docs/Report_lncRNA.pdf)

## **Project Structure**

The project is divided into several steps, each in its own directory.
Each directory refers to one step of the project. The description of the contents of each directory is indicated afterwards.
Bash and R scripts are enumerated regarding its execution in the workflow.

STEP 1 - [quality](/quality) :  Initial data quality checks
- Bash scripts
- [count_reads.txt](/quality/count_reads.txt) : Read counts per sample . Output file from [01_count_allreads_2.sh](/quality/01_count_allreads_2.sh)
- [output_fastqc](/quality/output_fastqc) : FastQC reports (HTMLs) for each sample and mate
  
STEP 2 - [alignment](/alignment) : Alignment to the reference genome
- Bash scripts
- [output_files](alignment/output_files) : summary of statistics and alignment rates of each sample
- SAM and BAM files can be found in the following directory in the cluster: /data/users/lfernandez/rnaseq_course/alignment/output_files

STEP 3 - [assembly](/assembly) :Transcriptome assembly using StringTie
- Bash scripts: [06_run_allStringTie.sh](/assembly/06_run_allStringTie.sh), [06_run_merge_StringTie.sh](/assembly/06_run_merge_StringTie.sh) and [07_count_transcripts.sh](/assembly/07_count_transcripts.sh)
- /output_files: [count_assembly.txt](assembly/count_assembly.txt) -> summary of the assembly
- meta-assembly GTF file can be found in the cluster: /data/users/lfernandez/rnaseq_course/assembly/output_files/meta_assembly.gtf
        
STEP 4 - [quantification](/quantification)
- Bash scripts: 
  - [09_transcriptome_fasta.sh](/quantification/09_transcriptome_fasta.sh) : create a transcriptome fasta file from meta-assembly.gtf and the kallisto indexed file
  - [10_kallisto_quantification.sh](/quantification/10_kallisto_quantification.sh) : run Kallisto per each sample
  - [11_summary_quantification.sh](/quantification/11_summary_quantification.sh)
- [step4_questions.txt](/quantification/step4_questions.txt) : summary with the total of estimated counts, and the number of transcripts, genes, novel transcripts and novel genes detected per sample.
- [transcript_expression_table.txt](/quantification/transcript_expression_table.txt) : file with transcript_id as row names and samples as column names. The expression levels are expressed in estimated counts.
- In the directory (/data/users/lfernandez/rnaseq_course/quantification/output_files) of the cluster you can find the transcriptome fasta files as well as the kallisto indexed file, and one directory per sample with the results from kallisto quantification.

STEP 5 - [differential_expression](/differential_expression)
- [08_table_mapping_transcriptsID.sh](/differential_expression/08_table_mapping_transcriptsID.sh) : bash script to create my_annotation_table.txt from the transcriptome assembly "meta_assembly.gtf".
- [12_sleuth.R](/differential_expression/12_sleuth.R) : R script to perform the differential expression analysis and subsequents plots. It requires the files: samples.txt and my_annotation_table.txt
  
- [/required_files](/required_files): 
  - [samples.txt](/required_files/samples.txt) : file that specifies the sample, condition and path to kallisto quantification output. It is required for the preparation step ('sleuth_prep') in Sleuth
  - [my_annotation_table.txt](/required_files/my_annotation_table.txt) : it compiles transcript_id, gene_id, gene_name and biotype information for each transcript present in the transcriptome assembly (meta_assembly.gtf). It is necessary to map each transcript to one gene due to the differential expression analysis is performed at transcript level.
  
- [differential_expression_table.txt](/required_files/differential_expression_table.txt) : Summary with transcript_id, gene_id, gene_name, biotype, log2FC and q-value.

  STEP 6: [integrative_analysis](/integrative_analysis)
- Bash scripts:
  - [13_create_bed.sh](/integrative_analysis/13_create_bed.sh) : to create two BED files (novel_transcripts.bed and annotated_transcripts.bed, they can be found in the cluster: /data/users/lfernandez/rnaseq_course/integrative_analysis/output_files ) from the meta_assembly.gtf
  - [14_intersect_intergenic.sh](/integrative_analysis/14_intersect_intergenic.sh) : to find the novel intergenic transcripts.
  - [15_bed_window.sh](/integrative_analysis/15_bed_window.sh) :  to create a BED files with a window of 100 nucleotides around the transcription start sites (TSS) and transcription end sites (TES) of the novel transcripts. Required for the next step.
  - [16_intersect_5prime.sh](/integrative_analysis/16_intersect_5prime.sh) : to find the overlaps of the TSS of the novel transcripts and the FANTOM CAGE of TSS.
  -  [17_intersect_3prime.sh](/integrative_analysis/16_intersect_5prime.sh) :to find the overlaps of the TES of the novel transcripts and the FANTOM CAGE of TES.
  -  [18_run_CPAT.sh](/integrative_analysis/18_run_CPAT.sh) : to estimate the protein coding potential of the novel transcripts.
  -  [19_step6_questions.sh](/integrative_analysis/19_step6_questions.sh) :  to answer the questions of the step 6.
- [20_prioritization.R](/integrative_analysis/20_prioritization.R) : R script to join information from DE and integrative analysis and make a sorted table of the transcripts. It also includes the code to generate the plots of the integrative analysis.

- [/output_files](/output_files) :
  - [step6_answers.txt](/output_files/step6_answers.txt) :  Statistics and porcentages of 5' and 3' well-annotated transcripts, protein coding potential transcritpts and intergentic transcripts.
  - [novel_intergenic_transcripts_plot.pdf](/output_files/novel_intergenic_transcripts_plot.pdf) : It shows the percentage of intergenic, 5’ well-annotated, 3’ well-annotated, 5’ and 3’ well-annotated and protein coding transcripts above the total of novel transcripts.
  - [novel_transcripts_plot.pdf](/output_files/novel_transcripts_plot.pdf) : deliverable. It shows the percentage of intergenic, 5’ well-annotated, 3’ well-annotated, 5’ and 3’ well-annotated and protein coding transcripts above the total of novel intergenic transcripts.
  - [prioritization_novel_table.tsv](/output_files/prioritization_novel_table.tsv) : list of novel transcripts sorted regarding the prioritization criteria.

 - All the BED files required can be found in the cluster: /data/users/lfernandez/rnaseq_course/integrative_analysis/output_files

 - All the reference files used can be found in the cluster: /data/courses/rnaseq_course/lncRNAs/Project1/references

 - For any questions or material, please feel free to reach me [lfercer.2014@gmail.com](mailto:lfercer.2014@gmail.com)




    
