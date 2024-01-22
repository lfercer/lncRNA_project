#DE SLEUTH

#install Sleuth
install.packages("BiocManager")
BiocManager::install("pachterlab/sleuth")
BiocManager::install('EnhancedVolcano')
devtools::install_github('kevinblighe/EnhancedVolcano')
# Install ggplot2
install.packages("ggplot2")
# Install ggrepel
install.packages("ggrepel")

library('ggrepel')
library(BiocManager)
library(sleuth)
library(EnhancedVolcano)

#Directory with all the input files: samples.txt, 13_sleuth.R, my_annotation_table.txt
working_dir <- "//wsl.localhost/Ubuntu/home/laurafernandez/rnaseq_course/quantification"

#Change to the working directory
setwd(working_dir)

# Load metadata about samples
samples <- read.table("samples.txt", header = TRUE, stringsAsFactors = FALSE)
samples
colnames(samples)

# Load transcript_id and gene_id table
transcripts_mapped <- read.table("my_annotation_table.txt", header = TRUE, sep = " ", fill = TRUE, 
                                 stringsAsFactors = FALSE)
colnames(transcripts_mapped)
#change name of the column to target_id
names(transcripts_mapped)[1] <- "target_id"
colnames(transcripts_mapped)

# Create Sleuth object

so <- sleuth_prep(samples, target_mapping = transcripts_mapped, formula = ~condition,
                  read_bootstrap_tpm = TRUE,
                  transform_fun_counts = function(x) log2(x + 0.5))
                  
# I think IM DOING PARENTAL VS PARACLONE
# Fit the model
so <- sleuth_fit(so, ~condition)
so <- sleuth_fit(so,~1,'reduced')

#Perform Wald test
so <- sleuth_wt(so, 'conditionparental',which_model = "full")


models(so)

#Extract Wald test result
sleuth_results <- sleuth_results(so, 'conditionparental',
                                test_type = "wt",
                                which_model = "full")

sleuth_results


# Filter results based on FDR threshold (adjust as needed)
#significant_results <- dplyr::filter(sleuth_results, qval <= 0.05)
#significant_results



#----------------------------------------------------
#VOLCANO PLOT SLEUTH
p.volcano <- plot_volcano(so, 'conditionparental',
                          test_type = "wt",
                          which_model = "full",
                          sig_level = 0.05)

# Extract gene IDs and coordinates from significant_results
gene_labels <- data.frame(x = highligthed_results$b,
                          y = -log10(highligthed_results$qval),
                          label = highligthed_results$gene_name)

p.volcano

# Add gene labels to the plot
p.volcano_labeled <- p.volcano + geom_text(data = gene_labels, aes(x = x, y = y, label = label), vjust = 1.5, hjust = 0)
p.volcano_labeled

#-------------------------------------------------------------------------------------------

# For table step 4
#sleuth function to get counts per gene per replicate:
counts_per_replicate <- sleuth_to_matrix(so, "obs_norm", "est_counts")
counts_per_replicate
# do you only have the transcripts?

#------------------------
#ENHANCED VOLCANO PLOT

#General volcano (no filters)
parental_paraclone_vol <- EnhancedVolcano(toptable = step5_table, 
                lab = step5_table$gene_name,
                x = "log2FC",
                y = "FDR",
                FCcutoff = 2,
                pCutoff = 0.05,
                title = 'Parental vs Paraclonal',
                labSize = 4.0,
                legendLabels=c('NS','Log2FC','p-value',
                               'p-value & Log2FC'),
                legendPosition = 'top',
                legendLabSize = 10,
                legendIconSize = 3.0)
parental_paraclone_vol




#----------------------------------

#Sample heat map
sample_heatmap <- plot_sample_heatmap(so)

#Transcripts heatmap
transcripts_heatmap<- plot_transcript_heatmap(so, sleuth_results$target_id)
transcripts_heatmap

#---------------------------------------------------------------------------------------
#Add novel in biotype to novel transcripts
transcripts_mapped$biotype[transcripts_mapped$biotype == ""] <- "novel"
sleuth_results$biotype[sleuth_results$biotype == ""] <- "novel"

#Deliverable table step 5:
step5_table <- data.frame(transcript_id = sleuth_results$target_id, gene_id = sleuth_results$gene_id, gene_name = sleuth_results$gene_name, biotype = sleuth_results$biotype,log2FC = sleuth_results$b, qval = sleuth_results$qval)
step5_table

#export file
write.table(step5_table,"differential_expression_table.txt", sep = " ", row.names=FALSE, col.names=TRUE)

#Novel transcripts volcano 






