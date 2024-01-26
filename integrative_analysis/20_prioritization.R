install.packages("dplyr")
library("dplyr")

#Set working directory where the input files are
setwd("//wsl.localhost/Ubuntu/home/laurafernandez/rnaseq_course/quantification")

#Import table from step5 (DE analysis from Sleuth)
table <- read.table("step5_table.txt", header = TRUE, stringsAsFactors = FALSE)

#Import results from step 6
annotation_5 <- read.table("overlap5prime.bed", header = FALSE, stringsAsFactors = FALSE)
annotation_3 <- read.table("overlap3prime.bed", header = FALSE, stringsAsFactors = FALSE)
intergenic <- read.table("novel_intergenic.bed", header = FALSE, stringsAsFactors = FALSE)
coding_potential <- read.table("//wsl.localhost/Ubuntu/home/laurafernandez/rnaseq_course/integrative_analysis/novel_coding_potential", header = TRUE, stringsAsFactors = FALSE,
                               sep = '\t')

# transcript_id : Extract characters before the first ":"
coding_potential$'transcript_id' <- sub("^(.*?):.*$", "\\1", rownames(coding_potential))


#subset for novel genes
novel_table <- table[grepl("novel", table$biotype),] %>%
  dplyr::select(transcript_id, gene_id, biotype, log2FC, qval)
  

#add columns: five_prime, three_prime, intergenic, prot_coding_potential
novel_table_all_info <- novel_table %>%
  mutate(five_prime = transcript_id %in% annotation_5$V4) %>%
  mutate(three_prime = transcript_id %in% annotation_3$V4) %>%
  mutate(intergenic = transcript_id %in% intergenic$V4) %>%
  left_join(coding_potential, by = "transcript_id") %>%
  mutate(mRNA_size = NULL, ORF_size = NULL, Fickett_score = NULL, Hexamer_score = NULL) 

write.table(novel_table_all_info, file = "//wsl.localhost/Ubuntu/home/laurafernandez/rnaseq_course/integrative_analysis/output_files/novel_complete_table.txt", col.names = TRUE, row.names = FALSE, sep = " ", quote = FALSE )


#Sort the table
#Filter out non-significant qval
novel_table_filtered <- novel_table_all_info %>%
  filter(qval < 0.05) %>%
  filter(!is.na(log2FC) & !is.na(qval)) %>%
  arrange(factor(biotype, levels = c("novel", levels(biotype))),
          desc(five_prime ), 
          desc(three_prime ), 
          desc(intergenic))

#Sort the table without filtering
novel_table_sorted <- novel_table_all_info %>%
  filter(!is.na(log2FC) & !is.na(qval)) %>%
  arrange(factor(biotype, levels = c("novel", levels(biotype))),
          desc(five_prime ), 
          desc(three_prime ), 
          desc(intergenic),
          qval)

write.table(novel_table_sorted, file = "//wsl.localhost/Ubuntu/home/laurafernandez/rnaseq_course/integrative_analysis/output_files/prioritization_novel_table.tsv", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE )


#-----------------------------------------------------------------
#Plots Integrative Analysis (step 6)

#Novel transcripts plot ---------------------------------------------
dev.off()
par(las = 2, mgp = c(2.5, 0.5, 0), mar = c(7, 4, 2, 2))

#data from step6_answers.txt
intergenic <- 1.809
prime_5 <- 82.230
prime_3 <- 81.351
protein_coding <- 7.740
prime_5_and_prime_3 <- 70.393


categories <- c("intergenic","5' annot.", "3' annot.", "5' and 3' annot.", "protein coding")
percentages <- c(intergenic, prime_5, prime_3, prime_5_and_prime_3, protein_coding)



barcolors = c("#800080", "skyblue", "skyblue", "skyblue", "skyblue")
barplot(percentages, names.arg = categories, col = barcolors, 
        ylim = c(0,100),
        main = 'Novel transcripts',
        ylab = "Percentages")



# Adding percentage labels on top of bars
text(x = barplot(percentages, col = barcolors, ylim = c(0, 100), axes = TRUE,
                 main = 'Novel transcripts',
                 ylab = "Percentage",
                 names.arg = categories),
     y = percentages + 0.02,
     labels = paste0(percentages, "%"),
     pos = 3, cex = 0.8, col = "black")

#Intergenic plot ---------------------------------------
sorted_table <- novel_table_sorted

total_intergenic <- sorted_table$transcript_id[sorted_table$intergenic == 'TRUE']
intergenic_and_5 <- sorted_table$transcript_id[(sorted_table$intergenic == 'TRUE') & (sorted_table$five_prime == 'TRUE')]
intergenic_and_3 <- sorted_table$transcript_id[(sorted_table$intergenic == 'TRUE') & (sorted_table$three_prime == 'TRUE')]
intergenic_and_prot <- sorted_table$transcript_id[(sorted_table$intergenic == 'TRUE') & (sorted_table$protein_coding > 0.364)]
intergenic_and_5_and_3 <- sorted_table$transcript_id[(sorted_table$intergenic == 'TRUE') & (sorted_table$five_prime == 'TRUE') & (sorted_table$three_prime == 'TRUE')]

#calculate percentages
total <- length(total_intergenic)
prime_5 <- length(intergenic_and_5)/total*100
prime_3 <- length(intergenic_and_3)/total*100
prime_5_3 <- length(intergenic_and_5_and_3)/total*100
protein_coding <- length(intergenic_and_prot)/total*100

categories <- c("5' annot.", "3' annot.", "5' and 3' annot.", "protein coding")
percentages <- c(round(prime_5, digits = 3), round(prime_3, digits = 3),round(prime_5_3, digits = 3) , protein_coding)


barplot(percentages, names.arg = categories, col = "#800080", 
        ylim = c(0,100),
        main = 'Novel intergenic transcripts',
        ylab = "Percentages")

par(las = 2, mgp = c(2.5, 0.5, 0), mar = c(7, 4, 2, 2))

# Adding percentage labels on top of bars
text(x = barplot(percentages, col = "#800080", ylim = c(0, 100), axes = TRUE,
                 main = 'Novel intergenic transcripts',
                 ylab = "Percentage",
                 names.arg = categories),
     y = percentages + 0.02,
     labels = paste0(percentages, "%"),
     pos = 3, cex = 0.8, col = "black")
