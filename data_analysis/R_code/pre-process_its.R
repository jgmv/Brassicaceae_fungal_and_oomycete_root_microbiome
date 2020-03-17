### load packages and functions
require(ape)
require(vegan)

# load custom functions
source("R_code/functions.R")


### create output folder
if (!dir.exists("pre-processed_data")) dir.create("pre-processed_data")


### data input
cdm_raw <- read.table("raw_data/cdm_its.csv", h = T, row.names = 1)
tax_raw <- read.table("raw_data/taxonomy_its.csv", h = T, row.names = 1)
seq_raw <- read.dna("raw_data/ASV_seqs_its.fasta", format = "fasta")
samples <- read.table("raw_data/libraries_its.csv", h = T)


### rename samples
cdm_raw <- cdm_raw[as.character(samples$sequencing_library), ]
rownames(cdm_raw) <- samples$sample


### summarize number of reads per kingdom
n_reads <- t(apply(cdm_raw, 1, function(x) tapply(x, tax_raw$kingdom, sum)))
n_reads <- as.data.frame(colSums(n_reads))
n_reads$perc <- n_reads[, 1] * 100 / sum(n_reads)
colnames(n_reads) <- c("reads", "perc")
write.table(n_reads, "pre-processed_data/summary_reads_its.csv", sep = ";",
  col.names = NA)


### check and remove negative controls
nc <- cdm_raw[grep("H2O", rownames(cdm_raw)), ]
nc <- nc[, colSums(nc) > 0]

# remove control treatments, and ASVs therein
cdm_raw <- cdm_raw[grep("H2O", rownames(cdm_raw), invert = T), ]
tax_raw <- droplevels(tax_raw[!(rownames(tax_raw) %in% colnames(nc)), ])


### remove non-fungal ASVs
tax <- droplevels(tax_raw[tax_raw$kingdom == "Fungi", ])
cdm <- cdm_raw[, rownames(tax)]
seq <- seq_raw[rownames(tax)]


### rename ASVs
rownames(tax) <- gsub("ASV", "FU", rownames(tax))
colnames(cdm) <- gsub("ASV", "FU", colnames(cdm))
names(seq) <- gsub("ASV", "FU", names(seq))


### modify taxonomy file
tax_edited <- modify_taxonomy(tax)


### export data
write.table(cdm, file = "pre-processed_data/cdm_fungi.csv", sep = ";",
  col.names = NA)
write.table(tax_edited, file = "pre-processed_data/taxonomy_fungi.csv",
  sep = ";", col.names = NA)
write.dna(seq, file = "pre-processed_data/ASVs_fungi.fasta", format = "fasta",
  colsep = "")


### end
