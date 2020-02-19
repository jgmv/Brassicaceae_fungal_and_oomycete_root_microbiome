### load packages and functions
require(ape)
require(vegan)

# load custom functions
source("R_code/functions.R")


### create output folder
if (!dir.exists("pre-processed_data")) dir.create("pre-processed_data")


### data input
cdm_raw <- read.table("raw_data/cdm_fwd_cox.csv", h = T, row.names = 1)
tax_raw <- read.table("raw_data/taxonomy_cox.csv", h = T, row.names = 1)
seq_raw <- read.dna("raw_data/ASV_fwd_seqs_cox.fasta", format = "fasta")


### modify taxonomy and remove non-oomycete ASVs
colnames(tax_raw)[1] <- "domain"
colnames(tax_raw)[2] <- "supergroup"
tax <- droplevels(tax_raw[tax_raw$supergroup == "Stramenopiles", ])
cdm <- cdm_raw[, rownames(tax)]
seq <- seq_raw[rownames(tax)]


### rename ASVs
rownames(tax) <- gsub("ASV", "OO", rownames(tax))
colnames(cdm) <- gsub("ASV", "OO", colnames(cdm))
names(seq) <- gsub("ASV", "OO", names(seq))


### modify taxonomy file
tax_edited <- modify_taxonomy(tax)


### export data
write.table(cdm, file = "pre-processed_data/cdm_oomycota.csv", sep = ";",
  col.names = NA)
write.table(tax_edited, file = "pre-processed_data/taxonomy_oomycota.csv",
  sep = ";", col.names = NA)
write.dna(seq, file = "pre-processed_data/ASVs_oomycota.fasta",
  format = "fasta", colsep = "")


### end
