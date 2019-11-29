### script version: 2019-02-21


### load packages and functions
library(ape)
library(vegan)

# load custom functions
source("custom_functions.R")

# view package versions
print(sessionInfo())


### data input
cdm_raw <- read.table("cdm.csv", h = T, row.names = 1)
tax_raw <- read.table("taxonomy.csv", h = T, row.names = 1, sep = ";")
seq_raw <- read.dna("ASV_seqs.fasta", format = "fasta")
samples <- read.table("samples.txt", h = T, sep = "\t")


### rename samples
samples$sequencing_library %in% rownames(cdm_raw)
cdm_raw <- cdm_raw[as.character(samples$sequencing_library), ]
rownames(cdm_raw) <- samples$sample


### check and remove negative controls
nc <- cdm_raw[grep("H2O", rownames(cdm_raw)), ]
nc <- nc[, colSums(nc) > 0]
tax_raw[colnames(nc), ]

# remove control treatments, and ASVs therein
cdm_raw <- cdm_raw[grep("H2O", rownames(cdm_raw), invert = T), ]
tax_raw <- droplevels(tax_raw[!(rownames(tax_raw) %in% colnames(nc)), ])


### remove non fungal ASVs
tax <- droplevels(tax_raw[tax_raw$kingdom == "Fungi", ])
cdm <- cdm_raw[, rownames(tax)]
seq <- seq_raw[rownames(tax)]
nrow(tax)
ncol(cdm)
length(seq)


### rename ASVs
rownames(tax) <- gsub("ASV", "FU", rownames(tax))
colnames(cdm) <- gsub("ASV", "FU", colnames(cdm))
names(seq) <- gsub("ASV", "FU", names(seq))


### modify taxonomy file
tax_edited <- modify_taxonomy(tax)


### export data
if (!dir.exists("output")) dir.create("output")
write.table(cdm, file = "output/cdm_fungi.csv", sep = ";", col.names = NA)
write.table(tax_edited, file = "output/taxonomy_fungi.csv", sep = ";",
  col.names = NA)
write.dna(seq, file = "output/ASVs_fungi.fasta", format = "fasta", colsep = "")


### end
