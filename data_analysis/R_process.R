### load custom functions
source("R_code/functions.R")


### create output folder
if (!dir.exists("output")) dir.create("output")


### pre-process data
source("R_code/pre-process_its.R")
source("R_code/pre-process_cox.R")


### data input
sample <- read.csv("raw_data/samples.csv", h = T, sep = ";", row.names = 1)
cdm_fu <- read.csv("pre-processed_data/cdm_fungi.csv", h = T, sep = ";",
  row.names = 1)
tax_fu <- read.csv("pre-processed_data/taxonomy_fungi.csv", h = T, sep = ";",
  row.names = 1)
cdm_oo <- read.csv("pre-processed_data/cdm_oomycota.csv", h = T, sep = ";",
  row.names = 1)
tax_oo <- read.csv("pre-processed_data/taxonomy_oomycota.csv", h = T, sep = ";",
  row.names = 1)

# reorder treatments
rownames(cdm_fu) %in% rownames(sample)
cdm_fu <- cdm_fu[rownames(sample), ]

rownames(cdm_oo) %in% rownames(sample)
cdm_oo <- cdm_oo[rownames(sample), ]


### remove rare ASVs
cdm_fu <- remove_rare_ASVs(cdm_fu)
tax_fu <- droplevels(tax_fu[colnames(cdm_fu), ])
sample_fu <- droplevels(sample[rownames(cdm_fu), ])
cdm_oo <- remove_rare_ASVs(cdm_oo)
tax_oo <- droplevels(tax_oo[colnames(cdm_oo), ])
sample_oo <- droplevels(sample[rownames(cdm_oo), ])


### rarefaction curves
plot_rarefaction_curves(cdm_fu, sample_fu, "rarefaction_fungi.pdf")
plot_rarefaction_curves(cdm_oo, sample_oo, "rarefaction_oomycota.pdf")


### sampling depth
sample_fu$reads <- rowSums(cdm_fu)
sample_oo$reads <- rowSums(cdm_oo)

# plot sample depth
plot_sample_depth()

# test for differences in sample depth
kruskal.test(sample_fu$reads ~ sample_fu$compartment)
kruskal.test(sample_fu$reads ~ sample_fu$species)
wilcox.test(sample_fu$reads ~ sample_fu$date)
kruskal.test(sample_oo$reads ~ sample_oo$compartment)
kruskal.test(sample_oo$reads ~ sample_oo$species)
wilcox.test(sample_oo$reads ~ sample_oo$date)


### sample richness and diversity
# calculate richness and Shannon diversity
sample_fu <- calculate_diversity(cdm_fu, sample_fu)
sample_oo <- calculate_diversity(cdm_oo, sample_oo)

# explore richness vs reads
plot(sample_fu$reads ~ sample_fu$Sobs)
cor.test(sample_fu$reads, sample_fu$Sobs)
plot(sample_fu$reads ~ sample_fu$Sha)
cor.test(sample_fu$reads, sample_fu$Sha)

plot(sample_oo$reads ~ sample_oo$Sobs)
cor.test(sample_oo$reads, sample_oo$Sobs)
plot(sample_oo$reads ~ sample_oo$Sha)
cor.test(sample_oo$reads, sample_oo$Sha)

# test sample depth differences vs richness
kruskal.test(sample_fu$Sobs ~ sample_fu$compartment)
kruskal.test(sample_fu$Sobs ~ sample_fu$species)
wilcox.test(sample_fu$Sobs ~ sample_fu$date)
kruskal.test(sample_oo$Sobs ~ sample_oo$compartment)
kruskal.test(sample_oo$Sobs ~ sample_oo$species)
wilcox.test(sample_oo$Sobs ~ sample_oo$date)

# test sample depth differences vs richness
kruskal.test(sample_fu$Shan ~ sample_fu$compartment)
kruskal.test(sample_fu$Shan ~ sample_fu$species)
wilcox.test(sample_fu$Shan ~ sample_fu$date, paired = T)
kruskal.test(sample_oo$Shan ~ sample_oo$compartment)
kruskal.test(sample_oo$Shan ~ sample_oo$species)
wilcox.test(sample_oo$Shan ~ sample_oo$date)

# plot richness and diversity per sample
diversity_boxplots()
diversity_boxplots(div = T, "shannon_x_treatment.pdf")

# create tables summarizing reads and richness per sample
summarize_diversity(cdm_fu, sample_fu, "diversity_fungi.csv")
summarize_diversity(cdm_oo, sample_oo, "diversity_oomycetes.csv")

# compare richness across dates, accounting for missing samples in oomycota
wilcox.test(sample_fu$Sobs ~ sample_fu$date, paired = T)
wilcox.test(sample_fu$reads ~ sample_fu$date, paired = T)
wilcox.test(sample_fu$Shan ~ sample_fu$date, paired = T)

sample_oo_p <- pair_oo_samples()
wilcox.test(sample_oo_p$Sobs ~ sample_oo_p$date, paired = T)
wilcox.test(sample_oo_p$reads ~ sample_oo_p$date, paired = T)
wilcox.test(sample_oo_p$Shan ~ sample_oo_p$date, paired = T)

# diversity comparison plots
plot_diversity()
plot_diversity(div = T, "shannon_comparison.pdf")


### normalization of read abundance data
# normalize reads
cdm_vs_fu <- DSeq_normalization(cdm_fu, sample_fu, tax_fu)
cdm_vs_oo <- DSeq_normalization(cdm_oo, sample_oo, tax_oo)

# estimate mean-variance relationships
mean_var_relationship(cdm_fu)
mean_var_relationship(cdm_vs_fu)
mean_var_relationship(cdm_oo)
mean_var_relationship(cdm_vs_oo)

# prepare new sample and taxonomic data
tax_vs_fu <- droplevels(tax_fu[colnames(cdm_vs_fu), ])
sample_vs_fu <- droplevels(sample_fu[rownames(cdm_vs_fu), ])
tax_vs_oo <- droplevels(tax_oo[colnames(cdm_vs_oo), ])
sample_vs_oo <- droplevels(sample_oo[rownames(cdm_vs_oo), ])


### effects of factors on community structure
# calculate ecological distances
cdm_dist_fu <- vegdist(cdm_vs_fu, method = "bray")
cdm_dist_oo <- vegdist(cdm_vs_oo, method = "bray")

# calculate ecological distances per compartment
comp_dists_fu <- distance_per_compartment(cdm_vs_fu, sample_vs_fu)
comp_dists_oo <- distance_per_compartment(cdm_vs_oo, sample_vs_oo)

# calculate and plot PCoAs
pcoa_fu <- calculate_PCoA(cdm_dist_fu, "pcoa_fu.pdf")
pcoa_oo <- calculate_PCoA(cdm_dist_oo, "pcoa_oo.pdf")

# calculate and plot PCoAs per compartment
calculate_PCoAs_per_compartment(comp_dists_fu, "pcoa_fu")
calculate_PCoAs_per_compartment(comp_dists_oo, "pcoa_oo")

# calculate dbRDAs with explanatory factors
dbrda_fu <- calculate_dbRDA(cdm_dist_fu, cdm_vs_fu, sample_vs_fu,
  "dbrda_fu.pdf")
dbrda_oo <- calculate_dbRDA(cdm_dist_oo, cdm_vs_oo, sample_vs_oo,
  "dbrda_oo.pdf")

# calculate dbRDAs per compartment
calculate_dbRDAs_per_compartment(comp_dists_fu, cdm_vs_fu, sample_vs_fu,
  "dbrda_fu")
calculate_dbRDAs_per_compartment(comp_dists_oo, cdm_vs_oo, sample_vs_oo,
  "dbrda_oo")

# plot variance explained by dbRDAs
plot_dbRDA_var(dbrda_fu, "dbrda_var_fu.pdf")
plot_dbRDA_var_per_compartment("dbrda_fu_", "dbrda_var_compartment_fu.pdf")
plot_dbRDA_var(dbrda_oo, "dbrda_var_oo.pdf")
plot_dbRDA_var_per_compartment("dbrda_oo_", "dbrda_var_compartment_oo.pdf")


### taxonomy analysis
# calculate number of orders (no incerta sedis or unidentified)
number_of_taxa(tax_vs_fu)
number_of_taxa(tax_vs_oo)

# barplots with proportion of taxa
cdm_genus_fu <- taxon_barplot(cdm_vs_fu, sample_vs_fu, taxlevel = "genus",
  tax_vs_fu, n = 60, outfile = "genus_sample_fu.pdf")
cdm_order_fu <- taxon_barplot(cdm_vs_fu, sample_vs_fu, tax_vs_fu, n = 15, 
  outfile = "orders_sample_fu.pdf")
cdm_phylum_fu <- taxon_barplot(cdm_vs_fu, sample_vs_fu, tax_vs_fu,
  taxlevel = "phylum", n = 10, outfile = "phylum_sample_fu.pdf")
cdm_genus_oo <- taxon_barplot(cdm_vs_oo, sample_vs_oo, taxlevel = "genus",
  tax_vs_oo, n = 16, outfile = "genus_sample_oo.pdf")
cdm_order_oo <- taxon_barplot(cdm_vs_oo, sample_vs_oo, tax_vs_oo,
  outfile = "orders_sample_oo.pdf")

# distribution of reads across taxa
reads_per_taxa(cdm_vs_fu, tax_vs_fu, outfile = "orders_overall_fu.pdf")
reads_per_taxa(cdm_vs_fu, tax_vs_fu, "phylum",
  outfile = "phylum_overall_fu.pdf")
reads_per_taxa(cdm_vs_oo, tax_vs_oo, outfile = "orders_overall_oo.pdf")

# significance of taxa distribution across compartments
genus_sig_fu <- taxa_across_compartments(cdm_vs_fu, sample_vs_fu, tax_vs_fu,
  taxlevel = "genus", outfile = "genus_per_compartment_fu.csv")
order_sig_fu <- taxa_across_compartments(cdm_vs_fu, sample_vs_fu, tax_vs_fu,
  outfile = "orders_per_compartment_fu.csv")
genus_sig_oo <- taxa_across_compartments(cdm_vs_oo, sample_vs_oo, tax_vs_oo,
  taxlevel = "genus", outfile = "genus_per_compartment_oo.csv")
order_sig_oo <- taxa_across_compartments(cdm_vs_oo, sample_vs_oo, tax_vs_oo,
  outfile = "orders_per_compartment_oo.csv")

shared_ASVs <- function(cdm, sam, filename = "venn") {
  require(gplots)
  asv <- sapply(as.character(unique(sam$compartment)), function(x) NULL)
  for(i in names(asv)) {
    cdm_temp <- cdm[sam$compartment == i, ]
    cdm_temp <- cdm_temp[, colSums(cdm_temp) > 0]
    asv[[i]] <- colnames(cdm_temp)
  }
  pdf(paste0("output/", filename, ".pdf"))
  plot(venn(asv), names = c("Root", "Endosphere", "Rhizosphere", "Bulk soil",
    "Root zone soil"))
  dev.off()
}
shared_ASVs(cdm_vs_fu, sample_vs_fu)

### enrichment of ASVs across compartments
enrichment_analysis(cdm_fu, sample_fu, tax_fu, cdm_order_fu,
  outfile = "dif_abund_fu")
enrichment_analysis(cdm_oo, sample_oo, tax_oo, cdm_order_oo,
  outfile = "dif_abund_oo")


### end
