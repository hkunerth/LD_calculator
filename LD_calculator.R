# tool to calculate LD between pairs of SNPS

library(tidyverse)
library(seqinr)
library(plyr)

setwd("Path to working directory here")


### Part 1: filtering loci 


locus_X <- read.alignment("locus_X.aln", format = "FASTA") # .aln (Clustal) formatted alignment file with each sample and the locus of interest

count <- 0 # set counter to 0
locus_X_matrix <- as.matrix.alignment(locus_X) # reformat as matrix to work with each character
locus_X_filtered_matrix <- matrix(, nrow = nrow(locus_X_matrix), ncol = 0) # create empty matrix for storage
for (i in 1:(ncol(locus_X_matrix)))
{
  freq_table <- count(locus_X_matrix[, i])
  min_freq <- min(freq_table$freq)
   if(min_freq > 8){     if(nrow(freq_table) > 1){       if(freq_table[1, 1] != '-'){ # removing SNPs with frequency below 8, invariant sites, and indels
    count <- count + 1
    locus_X_filtered_matrix <- cbind(locus_X_filtered_matrix, locus_X_matrix[, i])}}} # bind resulting matrix
}

write.table(locus_X_filtered_matrix, file = "locus_X_filtered_matrix.csv", sep = '	', quote = FALSE) # write out filtered table 

# repeat for each locus

### filter by sample

# filter to samples of interest, in this case by population
Dover_females <- c("DOV_10_1", "DOV_11_1", "DOV_113_1", "DOV_2_1", "DOV_22_1", "DOV_24_1", "DOV_26_1", "DOV_29_1", "DOV_34_1", "DOV_35_1", "DOV_36_1", "DOV_37_1", "DOV_58_1", "DOV_60_1", "DOV_61_1", "DOV_62_1", "DOV_63_1", "DOV_67_1", "DOV_68_1", "DOV_69_1", "DOV_70_1", "DOV_71_1")


# LD between two loci, X and Y. substitute loci of interest or iterate over entire dataset using this as a function

locus_X_Dover_f <- subset(locus_X_filtered_matrix, rownames(locus_X_filtered_matrix) %in% Dover_females) # grab SNPs
locus_Y_Dover_f <- subset(locus_Y_filtered_matrix, rownames(locus_Y_filtered_matrix) %in% Dover_females)
locus_X_Y_Dover_f <- merge(locus_X_Dover_f, locus_Y_Dover_f, by = 0) # combine into one dataframe

tt_1 <- table(locus_X_Y_Dover_f[, 2]) # make a test table
A <- names(which.max(tt_1)) # major allele
freq_A <- tt_1[[A]] / sum(tt_1) # calculate frequency
freq_a <- 1 - freq_A # minor allele frequency

tt_2 <- table(locus_X_Y_Dover_f[, 3]) 
B <- names(which.max(tt_2))
freq_B <- tt_2[[B]] / sum(tt_2)
freq_b <- 1 - freq_B

hap_table <- locus_X_Y_Dover_f[, c(2,3)] # create a haplotype table

ht <- table(hap_table) # format as R table
ht[[A, B]] 
freq_AB <- ht[[A, B]] / sum(ht) # get frequency of AB

D_AB <- freq_AB - (freq_A*freq_B) # calculate D

r_2 <- (D_AB*D_AB)/(freq_A*freq_a*freq_B*freq_b) # calculate R_2

output <- c(Locus_X, Locus_Y, r_2) # bind results - add loop if iterating. 
