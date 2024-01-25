install.packages("dplyr")
library("dplyr")

#Set working directory
setwd("//wsl.localhost/Ubuntu/home/laurafernandez/rnaseq_course/quantification")

#Import table from step5 (DE analysis from Sleuth)
table <- read.table("step5_table.txt", header = TRUE, stringsAsFactors = FALSE)

#Import results from step 6
annotation_5 <- read.table("overlap5prime.bed", header = FALSE, stringsAsFactors = FALSE)
annotation_3 <- read.table("overlap3prime.bed", header = FALSE, stringsAsFactors = FALSE)
intergenic <- read.table("novel_intergenic.bed", header = FALSE, stringsAsFactors = FALSE)
coding_potential <- read.table("novel_coding_potential.dat", header = FALSE, stringsAsFactors = FALSE)


#add empty columns to table 
table$'five_prime' <- c(rep(NA, length.out=length(table[,1])))
table$'three_prime' <- c(rep(NA, length.out=length(table[,1])))
table$'intergenic' <- c(rep(NA, length.out=length(table[,1])))
table$'protein_coding' <- c(rep(NA, length.out=length(table[,1])))


#create function that checks if transcript_id from table is in annotation_5
check5 <- function(row){
  if(row[1] %in% annotation_5$V4){
    return('x')
  } else {
    return('')
  }
}


#create function that checks if transcript_id from table is in annotation_3
check3 <- function(row){
  if(row[1] %in% annotation_3$V4){
    return('x')
  } else {
    return('')
  }
}

#create function that checks if transcript_id from table is in intergenic
check_intergenic <- function(row){
  if(row[1] %in% intergenic$V4){
    return('x')
  } else {
    return('')
  }
}

#create function that add the probability of coding region for a transcript_id
check_coding <- function(row){
  if(row[1] %in% coding_potential$V1){
    return(coding_potential$V6[coding_potential$V1 == row[1]])
  } else {
    return('')
  }
}

#Apply functions for all transcripts
table$five_prime <- apply(table,1,check5)
table$three_prime <- apply(table,1,check3)
table$intergenic <- apply(table,1,check_intergenic)
table$protein_coding <- apply(table,1,check_coding)

write.table(table, file = "complete_table.txt", col.names = TRUE, row.names = FALSE, sep = " " )

#------------
# Sorting the table
#Exclude transcripts with "NA" in log2FC and qval

sorted_table <- table %>%
  arrange(factor(biotype, levels = c("novel", levels(biotype))),
          desc(five_prime == 'x'), 
          desc(three_prime == 'x'), 
          desc(intergenic == 'x'),
          desc(abs(log2FC)),
          qval) %>%
  filter(!is.na(log2FC) & !is.na(qval))

write.table(sorted_table, file = "sorted_table.tsv", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE )



#filter out transcripts with non significant qval
sorted_table2 <- table %>%
  arrange(factor(biotype, levels = c("novel", levels(biotype))),
          desc(five_prime == 'x'), 
          desc(three_prime == 'x'), 
          desc(intergenic == 'x'),
          desc(abs(log2FC)),
          qval) %>%
  filter(!is.na(log2FC) & !is.na(qval)) %>%
  filter(qval < 0.05)

sorted_table2





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



#Intergenic plot ----------------------------

total_intergenic <- sorted_table$transcript_id[sorted_table$intergenic == 'x']
intergenic_and_5 <- sorted_table$transcript_id[(sorted_table$intergenic == 'x') & (sorted_table$five_prime == 'x')]
intergenic_and_3 <- sorted_table$transcript_id[(sorted_table$intergenic == 'x') & (sorted_table$three_prime == 'x')]
intergenic_and_prot <- sorted_table$transcript_id[(sorted_table$intergenic == 'x') & (sorted_table$protein_coding > 0.364)]
intergenic_and_5_and_3 <- sorted_table$transcript_id[(sorted_table$intergenic == 'x') & (sorted_table$five_prime == 'x') & (sorted_table$three_prime == 'x')]

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

