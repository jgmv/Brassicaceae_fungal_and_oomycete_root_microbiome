### script version: 2019-02-21


### load packages and functions
library(ape)
library(vegan)

# load custom functions
source("custom_functions.R")

# view package versions
print(sessionInfo())


### data input
cdm_raw <- read.table("cdm_fwd.csv", h = T, row.names = 1)
tax_raw <- read.table("taxonomy.csv", h = T, row.names = 1, sep = ";")
seq_raw <- read.dna("ASV_fwd_seqs.fasta", format = "fasta")


### modify taxonomy and remove unwanted ASVs
colnames(tax_raw)[1] <- "domain"
colnames(tax_raw)[2] <- "supergroup"
tax <- droplevels(tax_raw[tax_raw$supergroup == "Stramenopiles", ])
nrow(tax_raw) - nrow(tax)
cdm <- cdm_raw[, rownames(tax)]
ncol(cdm_raw) - ncol(cdm)
seq <- seq_raw[rownames(tax)]
length(seq_raw) - length(seq)


### rename ASVs
rownames(tax) <- gsub("ASV", "OO", rownames(tax))
colnames(cdm) <- gsub("ASV", "OO", colnames(cdm))
names(seq) <- gsub("ASV", "OO", names(seq))


### modify taxonomy file
tax_edited <- modify_taxonomy(tax)


### export data
if (!dir.exists("output")) dir.create("output")
write.table(cdm, file = "output/cdm_oomycota.csv", sep = ";", col.names = NA)
write.table(tax_edited, file = "output/taxonomy_oomycota.csv", sep = ";",
  col.names = NA)
write.dna(seq, file = "output/ASVs_oomycota.fasta", format = "fasta",
  colsep = "")


### end
