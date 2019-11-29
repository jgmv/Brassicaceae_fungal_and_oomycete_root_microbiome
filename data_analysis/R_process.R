### version: 2019-02-22


### load packages and functions
library(DESeq)
library(gplots)
library(Hmisc)
library(igraph)
library(indicspecies)
#library(iNEXT)
#library(metacoder)
library(mvabund)
#library(phyloseq)
library(RColorBrewer)
library(statmod)
library(vegan)
library(viridis)

# load custom functions
source("input/custom_functions.R")

# view package versions
sessionInfo()


### data input
sample <- read.csv("input/samples.csv", h = T, sep = ";", row.names = 1)
cdm_fu <- read.csv("input/cdm_fungi.csv", h = T, sep = ";", row.names = 1)
tax_fu <- read.csv("input/taxonomy_fungi.csv", h = T, sep = ";", row.names = 1)
cdm_oo <- read.csv("input/cdm_oomycota.csv", h = T, sep = ";", row.names = 1)
tax_oo <- read.csv("input/taxonomy_oomycota.csv", h = T, sep = ";",
  row.names = 1)

# reorder treatments
rownames(cdm_fu) %in% rownames(sample)
cdm_fu <- cdm_fu[rownames(sample), ]
rownames(cdm_oo) %in% rownames(sample)
cdm_oo <- cdm_oo[rownames(sample), ]


### create output folder
if (!dir.exists("output")) dir.create("output")


### create color and symbol codes for treatments
# color key for species
col_spec <- brewer.pal(length(levels(sample$sp_code)), "Set1")
names(col_spec) <- levels(sample$sp_code)
sample$col_spec <- col_spec[sample$sp_code]

# symbol key for species
pch_spec <- unique(as.numeric(sample$sp_code)) + 20
names(pch_spec) <- levels(sample$sp_code)
sample$pch_spec <- pch_spec[sample$sp_code]

# color key for compartment
col_comp <- brewer.pal(length(levels(sample$compartment)), "Set2")
names(col_comp) <- levels(sample$compartment)
sample$col_comp <- col_comp[sample$compartment]


### remove rare ASVs
# Fungi
n_min_fu <- 5
cdm_fu <- cdm_fu[colSums(cdm_fu) >= n_min_fu]
tax_fu <- droplevels(tax_fu[colnames(cdm_fu), ])
cdm_fu <- cdm_fu[rowSums(cdm_fu) > 0, ]
sample_fu <- droplevels(sample[rownames(cdm_fu), ])
ncol(cdm_fu)
nrow(cdm_fu)

# Oomycota
n_min_oo <- 5
cdm_oo <- cdm_oo[colSums(cdm_oo) >= n_min_oo]
tax_oo <- droplevels(tax_oo[colnames(cdm_oo), ])
cdm_oo <- cdm_oo[rowSums(cdm_oo) > 0, ]
sample_oo <- droplevels(sample[rownames(cdm_oo), ])
ncol(cdm_oo)
nrow(cdm_oo)


### sample depth and diversity
## rarefaction curves
# Fungi
pdf("output/rarefaction_fungi.pdf", w = 9, h = 9, pointsize = 12)
par(mfrow = c(3, 2), mar = c(4, 4, 1, 1), las = 1, cex = 1.25, lwd = 1.5)
rarecurve(t(colSums(cdm_fu)), step = 100, label = F, ylab = "Number of OTUs",
  xlab = "Reads", axes = F, main = "overall")
axis(1, pos = 0, lwd = 1.5)
axis(2, pos = 0, lwd = 1.5)
legend("bottomright", levels(sample_fu$sp_code), col = col_spec, lty = 1,
  lwd = 1.5, bty = "n")
for(i in levels(sample_fu$compartment)[c(1, 5, 3, 4, 2)]) {
  rarecurve(cdm_fu[sample_fu$compartment == i, ], step = 100,
    label = F, ylab = "Number of OTUs", xlab = "Reads", axes = F,
    ylim = c(0, max(apply(cdm_fu[sample_fu$compartment == i, ],
      1, function(x) sum(x > 0))) * 1.05),
    xlim = c(0, max(rowSums(cdm_fu[sample_fu$compartment == i, ])) *
      1.05), col = col_spec[sample_fu$sp_code[sample_fu$compartment == i]],
      main = i)
  axis(1, pos = 0, lwd = 1.5)
  axis(2, pos = 0, lwd = 1.5)
}
dev.off()

# Oomycota
pdf("output/rarefaction_oomycota.pdf", w = 9, h = 9, pointsize = 12)
par(mfrow = c(3, 2), mar = c(4, 4, 1, 1), las = 1, cex = 1.25, lwd = 1.5)
rarecurve(t(colSums(cdm_oo)), step = 100, label = F, ylab = "Number of OTUs",
  xlab = "Reads", axes = F, main = "overall")
axis(1, pos = 0, lwd = 1.5)
axis(2, pos = 0, lwd = 1.5)
legend("bottomright", levels(sample_oo$sp_code), col = col_spec, lty = 1,
  lwd = 1.5, bty = "n")
for(i in levels(sample_oo$compartment)[c(1, 5, 3, 4, 2)]) {
  rarecurve(cdm_oo[sample_oo$compartment == i, ], step = 100,
    label = F, ylab = "Number of OTUs", xlab = "Reads", axes = F,
    ylim = c(0, max(apply(cdm_oo[sample_oo$compartment == i, ],
      1, function(x) sum(x > 0))) * 1.05),
    xlim = c(0, max(rowSums(cdm_oo[sample_oo$compartment == i, ])) *
      1.05), col = col_spec[sample_oo$sp_code[sample_oo$compartment == i]],
      main = i)
  axis(1, pos = 0, lwd = 1.5)
  axis(2, pos = 0, lwd = 1.5)
}
dev.off()

## sample depth
sample_fu$reads <- rowSums(cdm_fu)
sample_oo$reads <- rowSums(cdm_oo)

# plot sample depth
pdf("output/reads_x_treatment.pdf", w = 8, h = 5, pointsize = 12)
layout(matrix(c(rep(1, 5), rep(2, 4), rep(3, 3),
  rep(4, 5), rep(5, 4), rep(6, 3)), byrow = T, nrow = 2))
boxplot_pt(sample_fu$compartment, sample_fu$reads, log = "y",richness las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Compartment", ylab = "Reads [log(x)]", lwd = 1.25)
boxplot_pt(sample_fu$species, sample_fu$reads, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Host species", ylab = "Reads [log(x)]", lwd = 1.25)
boxplot_pt(sample_fu$date, sample_fu$reads, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Date", ylab = "Reads [log(x)]", lwd = 1.25)
boxplot_pt(sample_oo$compartment, sample_oo$reads + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Compartment", ylab = "Reads [log(x)]", lwd = 1.25)
boxplot_pt(sample_oo$species, sample_oo$reads + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Host species", ylab = "Reads [log(x)]", lwd = 1.25)
boxplot_pt(sample_oo$date, sample_oo$reads + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Date", ylab = "Reads [log(x)]", lwd = 1.25)
dev.off()

# test sample depth differences
kruskal.test(sample_fu$reads ~ sample_fu$compartment)
kruskal.test(sample_fu$reads ~ sample_fu$species)
wilcox.test(sample_fu$reads ~ sample_fu$date)
kruskal.test(sample_oo$reads ~ sample_oo$compartment)
kruskal.test(sample_oo$reads ~ sample_oo$species)
wilcox.test(sample_oo$reads ~ sample_oo$date)


# richness vs reads
plot(sample_fu$reads ~ sample_fu$Sobs)
cor.test(sample_fu$reads, sample_fu$Sobs)
plot(sample_fu$reads ~ sample_fu$Sha)
cor.test(sample_fu$reads, sample_fu$Sha)

plot(sample_oo$reads ~ sample_oo$Sobs)
cor.test(sample_oo$reads, sample_oo$Sobs)
plot(sample_oo$reads ~ sample_oo$Sha)
cor.test(sample_oo$reads, sample_oo$Sha)

## richness
sample_fu$Sobs <- apply(cdm_fu, 1, function(x) sum(x > 0))
sample_oo$Sobs <- apply(cdm_oo, 1, function(x) sum(x > 0))

# plot richness
pdf("output/richness_x_treatment.pdf", w = 8, h = 5, pointsize = 12)
layout(matrix(c(rep(1, 5), rep(2, 4), rep(3, 3),
  rep(4, 5), rep(5, 4), rep(6, 3)), byrow = T, nrow = 2))
boxplot_pt(sample_fu$compartment, sample_fu$Sobs, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Compartment", ylab = "Richness [log(S)]", lwd = 1.25)
boxplot_pt(sample_fu$species, sample_fu$Sobs, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Host species", ylab = "Richness [log(S)]", lwd = 1.25)
boxplot_pt(sample_fu$date, sample_fu$Sobs, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Date", ylab = "Richness [log(S)]", lwd = 1.25)
boxplot_pt(sample_oo$compartment, sample_oo$Sobs + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Compartment",
  ylab = "Richness [log(S + 1)]", lwd = 1.25)
boxplot_pt(sample_oo$species, sample_oo$Sobs + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Host species",
  ylab = "Richness [log(S + 1)]", lwd = 1.25)
boxplot_pt(sample_oo$date, sample_oo$Sobs + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Date",
  ylab = "Richness [log(S + 1)]", lwd = 1.25)
dev.off()

# test sample depth differences
kruskal.test(sample_fu$Sobs ~ sample_fu$compartment)
kruskal.test(sample_fu$Sobs ~ sample_fu$species)
wilcox.test(sample_fu$Sobs ~ sample_fu$date)
kruskal.test(sample_oo$Sobs ~ sample_oo$compartment)
kruskal.test(sample_oo$Sobs ~ sample_oo$species)
wilcox.test(sample_oo$Sobs ~ sample_oo$date)

## Shannon diversity (effective species)
sample_fu$Shan <- exp(diversity(cdm_fu, index = "shannon"))
sample_oo$Shan <- exp(diversity(cdm_oo, index = "shannon"))

# plot Shannon diversity
pdf("output/shannon_x_treatment.pdf", w = 8, h = 5, pointsize = 12)
layout(matrix(c(rep(1, 5), rep(2, 4), rep(3, 3),
  rep(4, 5), rep(5, 4), rep(6, 3)), byrow = T, nrow = 2))
boxplot_pt(sample_fu$compartment, sample_fu$Shan, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Compartment", ylab = "Shannon diversity [log(ES)]", lwd = 1.25)
boxplot_pt(sample_fu$species, sample_fu$Shan, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Host species", ylab = "Shannon diversity [log(ES)]", lwd = 1.25)
boxplot_pt(sample_fu$date, sample_fu$Shan, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
  xlab = "Date", ylab = "Shannon diversity [log(ES)]", lwd = 1.25)
boxplot_pt(sample_oo$compartment, sample_oo$Shan + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Compartment",
  ylab = "Shannon diversity [log(ES) + 1]", lwd = 1.25)
boxplot_pt(sample_oo$species, sample_oo$Shan + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Host species",
  ylab = "Shannon diversity [log(ES) + 1]", lwd = 1.25)
boxplot_pt(sample_oo$date, sample_oo$Sobs + 1, log = "y", las = 1,
  border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
  main = "Oomycota", xlab = "Date",
  ylab = "Shannon diversity [log(ES) + 1]", lwd = 1.25)
dev.off()

# test Shannon diversity differences
kruskal.test(sample_fu$Shan ~ sample_fu$compartment)
kruskal.test(sample_fu$Shan ~ sample_fu$species)
wilcox.test(sample_fu$Shan ~ sample_fu$date, paired = T)
kruskal.test(sample_oo$Shan ~ sample_oo$compartment)
kruskal.test(sample_oo$Shan ~ sample_oo$species)
wilcox.test(sample_oo$Shan ~ sample_oo$date)

## table with reads and richness
tab_fu <- data.frame(
  n = as.character(table(sample_fu$sample)),
  reads_t = rowSums(apply(cdm_fu, 2,
    function(x) tapply(x, sample_fu$sample, sum))),
  sobs_t = specnumber(cdm_fu, sample_fu$sample),
  reads = paste(round(tapply(sample_fu$reads, sample_fu$sample, mean), 1), "±",
    round(tapply(sample_fu$reads, sample_fu$sample, sd), 1)),
  sobs_s = paste(round(tapply(sample_fu$Sobs, sample_fu$sample, mean), 1), "±",
    round(tapply(sample_fu$Sobs, sample_fu$sample, sd), 1))
)
rownames(tab_fu) <- names(tapply(sample_fu$reads, sample_fu$sample, mean))
write.table(tab_fu, "output/diversity_fu.csv", sep = ";", col.names= NA)

tab_oo <- data.frame(
  n = as.character(table(sample_oo$sample)),
  reads_t = rowSums(apply(cdm_oo, 2,
    function(x) tapply(x, sample_oo$sample, sum))),
  sobs_t = specnumber(cdm_oo, sample_oo$sample),
  reads = paste(round(tapply(sample_oo$reads, sample_oo$sample, mean), 1), "±",
    round(tapply(sample_oo$reads, sample_oo$sample, sd), 1)),
  sobs_s = paste(round(tapply(sample_oo$Sobs, sample_oo$sample, mean), 1), "±",
    round(tapply(sample_oo$Sobs, sample_oo$sample, sd), 1))
)
rownames(tab_oo) <- names(tapply(sample_oo$reads, sample_oo$sample, mean))
write.table(tab_oo, "output/diversity_oo.csv", sep = ";", col.names= NA)


### comparing richness across dates
wilcox.test(sample_fu$Sobs ~ sample_fu$date, paired = T)
wilcox.test(sample_fu$reads ~ sample_fu$date, paired = T)
wilcox.test(sample_fu$Shan ~ sample_fu$date, paired = T)

x <- rownames(sample_fu)[!(rownames(sample_fu) %in% rownames(sample_oo))]
miss <- data.frame(reads = rep(0, length(x)), Sobs = rep(0, length(x)),
  Shan = rep(0, length(x)), date = c(rep("April", 5), rep("February", 4)))
rownames(miss) <- x
rm(x)
x <- rbind(sample_oo[, c("reads", "Sobs", "Shan", "date")], miss)
x <- x[order(x$date, rownames(x)),]
wilcox.test(x$Sobs ~ x$date, paired = T)
wilcox.test(x$reads ~ x$date, paired = T)
wilcox.test(x$Shan ~ x$date, paired = T)
rm(x)

# relationship between reads and richness / diversity
cor.test(sample_fu$Sobs, sample_fu$reads)
cor.test(sample_fu$Shan, sample_fu$reads)

cor.test(sample_oo$Sobs, sample_oo$reads)
cor.test(sample_oo$Shan, sample_oo$reads)


## diversity comparison plots
# Function for comparisons
plot_diversity <- function(sp, var, fac) {
  x <- fac[fac$sp_code == sp, ]
  t0 <- x$date == "February"
  t1 <- x$date == "April"
  comp_col <- brewer.pal(length(levels(x$compartment)), "Set2")
  names(comp_col) <- levels(x$compartment)
  plot(0, 0, type = "n", xlim = c(0.9, 2.1), ylim = c(0, max(fac[, var])),
    axes = F, xlab = NA, ylab = NA, main = sp)
  #points(jitter(c(rep(1, sum(t0)), rep(2, sum(t1))), amount = 0.1),
  #  c(x[t0, var], x[t1, var]), pch = 16,
  #  col = alpha(comp_col[x$compartment], 0.5))
  for(i in levels(x$compartment)) {
    y <- tapply(x[x$compartment == i, var],
      x[x$compartment == i, "date"], mean)[c(2, 1)]
    ysd <- tapply(x[x$compartment == i, var],
      x[x$compartment == i, "date"], sd)[c(2, 1)]
    polygon(c(1, 2, 2, 1, 1), c(y[1] - ysd[1], y[2] - ysd[2],
      y[2] + ysd[2], y[1] + ysd[1], y[1] - ysd[1]),
      col = alpha(comp_col[i], 0.25), border = F)
    #y <- tapply(x[x$compartment == i, var],
    #  x[x$compartment == i, "date"], median)[c(2, 1)]
    #iq25 <- tapply(x[x$compartment == i, var],
    #  x[x$compartment == i, "date"],
    #  function(x) as.vector(quantile(x, 0.25)))
    #iq75 <- tapply(x[x$compartment == i, var],
    #  x[x$compartment == i, "date"],
    #  function(x) quantile(x, 0.75))
    #polygon(c(1, 2, 2, 1, 1), c(y[1] - iq25[1], y[2] - iq75[1],
    #  y[2] + iq75[2], y[1] + iq25[2], y[1] - iq25[1]),
    #  col = alpha(comp_col[i], 0.25), border = F)
    lines(c(1, 2), y, type = "o", pch = 16, col = comp_col[i],
      cex = 1.2, lwd = 2)
  }
  axis(1, at = c(1, 2), labels = c("Feb.", "Ap."), lwd = 1.5)
  axis(2, lwd = 1.5)
}

# richness
pdf("output/richness_comparison.pdf", w = 7, h = 9, pointsize = 12)
layout(matrix(c(1:6, 7, 7, 7), ncol = 3, byrow = T))
par(mar = c(4, 3.5, 1, 0), cex = 1.25, las = 1)
plot_diversity("At", "Sobs", sample_fu)
mtext("Fungi", side = 2, line = 2.5, las = 0, cex = 1.25)
plot_diversity("Ch", "Sobs", sample_fu)
plot_diversity("Dv", "Sobs", sample_fu)
plot_diversity("At", "Sobs", sample_oo)
mtext("Oomycota", side = 2, line = 2.5, las = 0, cex = 1.25)
plot_diversity("Ch", "Sobs", sample_oo)
plot_diversity("Dv", "Sobs", sample_oo)
plot(0, 0, type = "n", axes = F, xlab = NA, ylab = NA)
legend("top", levels(sample$compartment), pch = 16, lty = 1, col = col_comp,
  bty = "n", pt.cex = 1.25, lwd = 2, ncol = 3, inset = -0.1)
dev.off()

# Shannon diversity
pdf("output/shannon_comparison.pdf", w = 7, h = 9, pointsize = 12)
layout(matrix(c(1:6, 7, 7, 7), ncol = 3, byrow = T))
par(mar = c(4, 3.5, 1, 0), cex = 1.25, las = 1)
plot_diversity("At", "Shan", sample_fu)
mtext("Fungi", side = 2, line = 2.5, las = 0, cex = 1.25)
plot_diversity("Ch", "Shan", sample_fu)
plot_diversity("Dv", "Shan", sample_fu)
plot_diversity("At", "Shan", sample_oo)
mtext("Oomycota", side = 2, line = 2.5, las = 0, cex = 1.25)
plot_diversity("Ch", "Shan", sample_oo)
plot_diversity("Dv", "Shan", sample_oo)
plot(0, 0, type = "n", axes = F, xlab = NA, ylab = NA)
legend("top", levels(sample$compartment), pch = 16, lty = 1, col = col_comp,
  bty = "n", pt.cex = 1.25, lwd = 2, ncol = 3, inset = -0.1)
dev.off()


### data normalization with DSeq, adapted from McMurdie and Holmes 2014
### (Protocol S1 - simulation-cluster-accuracy)
## estimate mean-variance relationships
# Fungi
mean_var_fu <- data.frame(mean = apply(cdm_fu, 2, mean),
  var = apply(cdm_fu, 2, function(x) sd(x) * sd(x)))
plot(mean_var_fu$var ~ mean_var_fu$mean, log = "xy", xlab = "mean",
  ylab = "variance")
abline(a = 0, b = 1, lty = 2)

# Oomycota
mean_var_oo <- data.frame(mean = apply(cdm_oo, 2, mean),
  var = apply(cdm_oo, 2, function(x) sd(x) * sd(x)))
plot(mean_var_oo$var ~ mean_var_oo$mean, log = "xy", xlab = "mean",
  ylab = "variance")
abline(a = 0, b = 1, lty = 2)


## normalize reads
# Fungi
dds_fu <- newCountDataSet(t(cdm_fu) + 1, conditions = sample_fu,
  featureData = AnnotatedDataFrame(tax_fu)) # x + 1 to avoid errors 
dds_fu <- estimateSizeFactors(dds_fu)
sizeFactors(dds_fu)
dds_fu <- estimateDispersions(dds_fu, method = "blind",
  sharingMode = "maximum", fitType = "local") # change sharingMode to gene-est-only if n > 7
dispcol <- grep("disp\\_", colnames(fData(dds_fu)))
if (any(!is.finite(fData(dds_fu)[, dispcol]))) {
  fData(dds_fu)[which(!is.finite(fData(dds_fu)[, dispcol])), dispcol] <- 0
}
rm(dispcol)
cdm_vs_fu <- t(exprs(varianceStabilizingTransformation(dds_fu)))

# Oomycota
dds_oo <- newCountDataSet(t(cdm_oo) + 1, conditions = sample_oo,
  featureData = AnnotatedDataFrame(tax_oo)) # x + 1 to avoid errors 
dds_oo <- estimateSizeFactors(dds_oo)
sizeFactors(dds_oo)
dds_oo <- estimateDispersions(dds_oo, method = "blind",
  sharingMode = "maximum", fitType = "local") # change sharingMode to gene-est-only if n > 7
dispcol <- grep("disp\\_", colnames(fData(dds_oo)))
if (any(!is.finite(fData(dds_oo)[, dispcol]))) {
  fData(dds_oo)[which(!is.finite(fData(dds_oo)[, dispcol])), dispcol] <- 0
}
rm(dispcol)
cdm_vs_oo <- t(exprs(varianceStabilizingTransformation(dds_oo)))

# set negative values (corresponding to zeroes in raw data set) to zero
length(cdm_vs_fu[cdm_vs_fu < 0])
length(cdm_fu[cdm_fu == 0])
cdm_vs_fu[cdm_vs_fu < 0.0] <- 0.0
cdm_vs_fu <- cdm_vs_fu[, colSums(cdm_vs_fu) > 0]
tax_vs_fu <- droplevels(tax_fu[colnames(cdm_vs_fu), ])
cdm_vs_fu <- cdm_vs_fu[rowSums(cdm_vs_fu) > 0, ]
sample_vs_fu <- droplevels(sample_fu[rownames(cdm_vs_fu), ])

length(cdm_vs_oo[cdm_vs_oo < 0])
length(cdm_oo[cdm_oo == 0])
cdm_vs_oo[cdm_vs_oo < 0.0] <- 0.0
cdm_vs_oo <- cdm_vs_oo[, colSums(cdm_vs_oo) > 0]
tax_vs_oo <- droplevels(tax_oo[colnames(cdm_vs_oo), ])
cdm_vs_oo <- cdm_vs_oo[rowSums(cdm_vs_oo) > 0, ]
sample_vs_oo <- droplevels(sample_oo[rownames(cdm_vs_oo), ])


### effects of factors on community structure
## calculate ecological distances
# Fungi
cdm_dist_fu <- vegdist(cdm_vs_fu, method = "bray")
#cdm_dist_fu <- vegdist(decostand(cdm_vs_fu, method = "total"), method = "bray")
comp_dists_fu <- vector("list", length(levels(sample_vs_fu$compartment)))
names(comp_dists_fu) <- levels(sample_vs_fu$compartment)
for(i in names(comp_dists_fu)) {
  comp_dists_fu[[i]] <- vegdist(cdm_vs_fu[sample_vs_fu$compartment == i, ],
    method = "bray")
}

# Oomycota
cdm_dist_oo <- vegdist(cdm_vs_oo, method = "bray")
comp_dists_oo <- vector("list", length(levels(sample_vs_oo$compartment)))
names(comp_dists_oo) <- levels(sample_vs_oo$compartment)
for(i in names(comp_dists_oo)) {
  comp_dists_oo[[i]] <- vegdist(cdm_vs_oo[sample_vs_oo$compartment == i, ],
    method = "bray")
}

## calculate PCoAs
# Fungi
pcoa_fu <- cmdscale(cdm_dist_fu, eig = T, add = T)
pcoa_fu_expl <- round(eigenvals(pcoa_fu) / sum(eigenvals(pcoa_fu)) * 100, 1)
for(i in names(comp_dists_fu)) {
  name <- paste0("pcoa_", i, "_fu")
  pcoa <- cmdscale(comp_dists_fu[[i]], eig = T, add = T)
  expl <- round(eigenvals(pcoa) / sum(eigenvals(pcoa)) * 100, 1)
  assign(name, pcoa)
  assign(paste(name, "expl", sep = "_"), expl)
  rm(name, pcoa, expl)
}

# Oomycota
pcoa_oo <- cmdscale(cdm_dist_oo, eig = T, add = T)
pcoa_oo_expl <- round(eigenvals(pcoa_oo) / sum(eigenvals(pcoa_oo)) * 100, 1)
for(i in names(comp_dists_oo)) {
  name <- paste0("pcoa_", i, "_oo")
  pcoa <- cmdscale(comp_dists_oo[[i]], eig = T, add = T)
  expl <- round(eigenvals(pcoa) / sum(eigenvals(pcoa)) * 100, 1)
  assign(name, pcoa)
  assign(paste(name, "expl", sep = "_"), expl)
  rm(name, pcoa, expl)
}

# plot PCoAs
for(i in ls(pattern = "pcoa.*_expl")) {
  pcoa <- gsub("_expl", "", i)
  pdf(paste("output/", pcoa, ".pdf", sep = ""), w = 3, h = 3)
  par(mar = c(4, 4, 1, 1), las = 1)
  plot_pcoa(get(pcoa), get(i))
  dev.off()
  rm(pcoa)
}


## calculate dbRDA (when controlling for total abundance
## (adding  + Condition(rowSums(cdm_vs_fu))), very little is explained)
# Fungi
dbrda_fu <- dbrda(cdm_dist_fu ~ sample_vs_fu$sp_code +
  sample_vs_fu$compartment + sample_vs_fu$date,
  comm = cdm_vs_fu)
dbrda_fu_expl <- round(eigenvals(dbrda_fu) / sum(eigenvals(dbrda_fu)) * 100, 1)
(dbrda_fu_aov <- anova(dbrda_fu))
(dbrda_fu_aov_axis <- anova(dbrda_fu, by = "axis"))
(dbrda_fu_aov_margin <- anova(dbrda_fu, by = "margin"))
dbrda_fu_aov_var_expl <- (dbrda_fu_aov_margin$SumOfSqs /
  with(dbrda_fu, tot.chi)) * 100
names(dbrda_fu_aov_var_expl) <- c("Species", "Compartment", "Date", "Residuals")

# Oomycota
dbrda_oo <- dbrda(cdm_dist_oo ~ sample_vs_oo$sp_code +
  sample_vs_oo$compartment + sample_vs_oo$date,
  comm = cdm_vs_oo)
dbrda_oo_e <- round(eigenvals(dbrda_oo) / sum(eigenvals(dbrda_oo)) * 100, 1)
(dbrda_oo_aov <- anova(dbrda_oo))
(dbrda_oo_aov_axis <- anova(dbrda_oo, by = "axis"))
(dbrda_oo_aov_margin <- anova(dbrda_oo, by = "margin"))
dbrda_oo_aov_var_expl <- (dbrda_oo_aov_margin$SumOfSqs /
  with(dbrda_oo, tot.chi)) * 100
names(dbrda_oo_aov_var_expl) <- c("Species", "Compartment", "Date", "Residuals")


## dbRDA per compartment
# calculate dbRDAs per compartment
for(i in names(comp_dists_fu)) {
  s  <- droplevels(sample_vs_fu[sample_vs_fu$compartment == i, ])
  d  <- cdm_vs_fu[rownames(s), ]
  dd <- comp_dists_fu[[i]]
  #cp <- dbrda(dd ~ s$sp_code + s$date + Condition(s$dbrda_n), comm = d)
  cp <- dbrda(dd ~ s$sp_code + s$date, comm = d)
  cp_expl <- round(eigenvals(cp) / sum(eigenvals(cp)) * 100, 1)
  cp_aov <- anova(cp)
  cp_aov_axis <- anova(cp, by = "axis")
  cp_aov_margin <- anova(cp, by = "margin")
  if(is.null(cp$CA$imaginary.chi)) {
    cp_aov_var_expl <- (cp_aov_margin$SumOfSqs / 
      with(cp, tot.chi)) * 100
  } else {
    cp_aov_var_expl <- (cp_aov_margin$SumOfSqs / 
      with(cp, tot.chi - CA$imaginary.chi)) * 100
  }
  names(cp_aov_var_expl) <- c("Species", "Date", "Residuals")
  assign(paste("dbrda_fu", i, sep = "_"), cp)
  assign(paste("dbrda_fu_expl", i, sep = "_"), cp_expl)
  assign(paste("dbrda_fu_aov", i, sep = "_"), cp_aov)
  assign(paste("dbrda_fu_aov_axis", i, sep = "_"), cp_aov_axis)
  assign(paste("dbrda_fu_aov_margin", i, sep = "_"), cp_aov_margin)
  assign(paste("dbrda_fu_aov_var_expl", i, sep = "_"), cp_aov_var_expl)
  rm(cp, cp_expl, cp_aov, cp_aov_axis, cp_aov_margin, cp_aov_var_expl)
}


for(i in names(comp_dists_oo)) {
  s  <- droplevels(sample_vs_oo[sample_vs_oo$compartment == i, ])
  d  <- cdm_vs_oo[rownames(s), ]
  dd <- comp_dists_oo[[i]]
  #cp <- dbrda(dd ~ s$sp_code + s$date + Condition(s$dbrda_n), comm = d)
  cp <- dbrda(dd ~ s$sp_code + s$date, comm = d)
  cp_expl <- round(eigenvals(cp) / sum(eigenvals(cp)) * 100, 1)
  cp_aov <- anova(cp)
  cp_aov_axis <- anova(cp, by = "axis")
  cp_aov_margin <- anova(cp, by = "margin")
  if(is.null(cp$CA$imaginary.chi)) {
    cp_aov_var_expl <- (cp_aov_margin$SumOfSqs / 
      with(cp, tot.chi)) * 100
  } else {
    cp_aov_var_expl <- (cp_aov_margin$SumOfSqs / 
      with(cp, tot.chi - CA$imaginary.chi)) * 100
  }
  names(cp_aov_var_expl) <- c("Species", "Date", "Residuals")
  assign(paste("dbrda_oo", i, sep = "_"), cp)
  assign(paste("dbrda_oo_expl", i, sep = "_"), cp_expl)
  assign(paste("dbrda_oo_aov", i, sep = "_"), cp_aov)
  assign(paste("dbrda_oo_aov_axis", i, sep = "_"), cp_aov_axis)
  assign(paste("dbrda_oo_aov_margin", i, sep = "_"), cp_aov_margin)
  assign(paste("dbrda_oo_aov_var_expl", i, sep = "_"), cp_aov_var_expl)
  rm(cp, cp_expl, cp_aov, cp_aov_axis, cp_aov_margin, cp_aov_var_expl)
}

for(i in ls(pattern = "dbrda.*fu_expl")) {
  dbrda <- gsub("_expl", "", i)
  pdf(paste("output/", dbrda, ".pdf", sep = ""), w = 8, h = 8, pointsize = 36)
  par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
  plot_cap(get(dbrda), get(i))
  dev.off()
  rm(dbrda)
}

for(i in ls(pattern = "dbrda.*oo_expl")) {
  dbrda <- gsub("_expl", "", i)
  pdf(paste("output/", dbrda, ".pdf", sep = ""), w = 8, h = 8, pointsize = 36)
  par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
  plot_cap(get(dbrda), get(i))
  dev.off()
  rm(dbrda)
}


# plot variance explained by dbRDAs per compartment
pdf("output/dbRDA_var_fu.pdf", w = 6, h = 8, pointsize = 36)
par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
x <- barplot(rbind(dbrda_fu_aov_var_expl[1:3]), border = F,
  beside = T, las = 1, ylim = c(0, 25), ylab = "Variance explained (%)",
  xaxt = "n", yaxt = "n", space = c(rep(0.2, 3)))
grid(nx = 0, ny = 5, col = gray(0.5))
axis(2, lwd = 3)
text(x = x, y = -1, names(dbrda_fu_aov_var_expl[1:3]), xpd = TRUE, srt = 45,
  adj = 1)
box()
dev.off()

compartments <- c("Bulk soil", "Root zone", "Rhizosphere", "Root", "Endosphere")
pdf("output/dbRDA_var_compartment_fu.pdf", w = 8, h = 8, pointsize = 36)
par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
x <- barplot(cbind(dbrda_fu_aov_var_expl_bare_soil[1:2],
  dbrda_fu_aov_var_expl_root_soil[1:2],
  dbrda_fu_aov_var_expl_rhizosphere[1:2],
  dbrda_fu_aov_var_expl_root[1:2],
  dbrda_fu_aov_var_expl_endosphere[1:2]), beside = T, border = F,
  ylim = c(0, 25), xaxt = "n", yaxt = "n", ylab = "Variance explained (%)")
grid(nx = 0, ny = 5, col = gray(0.5))
axis(2, lwd = 3)
text(x = apply(x, 2, mean), y = -1, compartments, xpd = TRUE, srt = 45,
  adj = 1)
box()
dev.off()

pdf("output/dbRDA_var_oo.pdf", w = 6, h = 8, pointsize = 36)
par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
x <- barplot(rbind(dbrda_oo_aov_var_expl[1:3]), border = F,
  beside = T, las = 1, ylim = c(0, 25), ylab = "Variance explained (%)",
  xaxt = "n", yaxt = "n", space = c(rep(0.2, 3)))
grid(nx = 0, ny = 5, col = gray(0.5))
axis(2, lwd = 3)
text(x = x, y = -1, names(dbrda_oo_aov_var_expl[1:3]), xpd = TRUE, srt = 45,
  adj = 1)
box()
dev.off()

compartments <- c("Bulk soil", "Root zone", "Rhizosphere", "Root", "Endosphere")
pdf("output/dbRDA_var_compartment_oo.pdf", w = 8, h = 8, pointsize = 36)
par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
x <- barplot(cbind(dbrda_oo_aov_var_expl_bare_soil[1:2],
  dbrda_oo_aov_var_expl_root_soil[1:2],
  dbrda_oo_aov_var_expl_rhizosphere[1:2],
  dbrda_oo_aov_var_expl_root[1:2],
  dbrda_oo_aov_var_expl_endosphere[1:2]), beside = T, border = F,
  ylim = c(0, 25), xaxt = "n", yaxt = "n", ylab = "Variance explained (%)")
grid(nx = 0, ny = 5, col = gray(0.5))
axis(2, lwd = 3)
text(x = apply(x, 2, mean), y = -1, compartments, xpd = TRUE, srt = 45,
  adj = 1)
box()
dev.off()


### taxonomic summaries
# see https://grunwaldlab.github.io/metacoder_documentation/example.html
## create taxmaps
# Fungi
taxmap_fu <- parse_tax_data(cbind(tax_vs_fu, asv = rownames(tax_vs_fu)),
  class_cols = 1:4, named_by_rank = T,
  datasets = list(norm_abundance = as.data.frame(t(decostand(cdm_vs_fu,
  "total")))), mapping = c("{{name}}" = "{{name}}"))
taxmap_fu$data$tax_abund <- calc_taxon_abund(taxmap_fu, "norm_abundance")

# Oomycota
taxmap_oo <- parse_tax_data(cbind(tax_vs_oo, asv = rownames(tax_vs_oo)),
  class_cols = 1:4, named_by_rank = T,
  datasets = list(norm_abundance = as.data.frame(t(decostand(cdm_vs_oo,
  "total")))), mapping = c("{{name}}" = "{{name}}"))
taxmap_oo$data$tax_abund <- calc_taxon_abund(taxmap_oo, "norm_abundance")


## create differences between treatments
sample_vs_fu$diff_comp <- sample_vs_fu$compartment
levels(sample_vs_fu$diff_comp) <- c(levels(sample_vs_fu$diff_comp), "soil")
sample_vs_fu$diff_comp[sample_vs_fu$diff_comp == "root_soil" |
  sample_vs_fu$diff_comp == "bare_soil" |
  sample_vs_fu$diff_comp == "rhizosphere"] <- "soil"
taxmap_fu$data$diff_table <- compare_groups(taxmap_fu,
  data = "tax_abund", cols = rownames(sample_vs_fu),
  groups = sample_vs_fu$diff_comp)

sample_vs_oo$diff_comp <- sample_vs_oo$compartment
levels(sample_vs_oo$diff_comp) <- c(levels(sample_vs_oo$diff_comp), "soil")
sample_vs_oo$diff_comp[sample_vs_oo$diff_comp == "root_soil" |
  sample_vs_oo$diff_comp == "bare_soil" |
  sample_vs_oo$diff_comp == "rhizosphere"] <- "soil"
taxmap_oo$data$diff_table <- compare_groups(taxmap_oo,
  data = "tax_abund", cols = rownames(sample_vs_oo),
  groups = sample_vs_oo$diff_comp)


# plot overall heat tree
heat_tree_matrix(taxmap_fu, data = "diff_table", node_size = n_obs,
  node_label = taxon_names, node_color = log2_median_ratio,
  node_color_range = diverging_palette(), node_color_trans = "linear",
  node_color_interval = c(-3, 3), edge_color_interval = c(-3, 3),
  node_size_axis_label = "Number of OTUs", node_label_max = 25,
  node_size_range = c(0.01, 0.05), node_label_size_range = c(0.025, 0.05), 
  node_color_axis_label = "Log2 ratio median proportions",
  row_label_size = 25, col_label_size = 25, key_size = 0.6,
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  make_node_legend = F, make_edge_legend = F, 
  output_file = "output/heat_tree_compartment_fungi.png") 

heat_tree_matrix(taxmap_oo, data = "diff_table", node_size = n_obs,
  node_label = taxon_names, node_color = log2_median_ratio,
  node_color_range = diverging_palette(), node_color_trans = "linear",
  node_color_interval = c(-3, 3), edge_color_interval = c(-3, 3),
  node_size_axis_label = "Number of OTUs", node_label_max = 25,
  node_size_range = c(0.01, 0.05), node_label_size_range = c(0.025, 0.05), 
  node_color_axis_label = "Log2 ratio median proportions",
  row_label_size = 25, col_label_size = 25, key_size = 0.6,
  layout = "davidson-harel",
  initial_layout = "reingold-tilford",
  make_node_legend = F, make_edge_legend = F, 
  output_file = "output/heat_tree_compartment_oomycota.png") 




## barplots
# barplots per order and sample
order_tab_fu <- tax_prop_table(cdm_vs_fu, tax_vs_fu$order, sample_vs_fu$sample2,
  n = 15)
order_tab_fu <- order_tab_fu[c("Atba", "Atrs", "Atrh", "Atro", "Aten", "Chba",
  "Chrs", "Chrh", "Chro", "Chen", "Dvba", "Dvrs", "Dvrh", "Dvro", "Dven"), ]
#col_order_fu <- c(brewer.pal(ncol(order_tab_fu) - 1, "Paired"), gray(0.9))
col_order_fu <- color(ncol(order_tab_fu), start = 17)
names(col_order_fu) <- colnames(order_tab_fu)
pdf("output/orders_sample_fu.pdf", w = 6, h = 2.5, pointsize = 12)
par(mar = c(6, 4, 1, 6), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, -0.5))
x <- barplot(t(order_tab_fu), space = c(rep(c(0.3, rep(0.1, 4)), 3)), las = 1,
  col = col_order_fu, names.arg = rep(NA, nrow(order_tab_fu)),
  ylab = "Proportion", border = NA)
text(x, y = -0.1, rep(c("Bulk soil", "Root zone", "Rhizosphere", "Root",
  "Endosphere"), 3), xpd = TRUE, srt = 45, adj = 1)
lines(c(x[1], x[5]), c(-0.75, -0.75))
text((x[5] + x[1]) / 2, -0.85, "A. thaliana", font = 3, adj = 0.5)
lines(c(x[6], x[10]), c(-0.75, -0.75))
text((x[6] + x[10]) / 2, -0.85, "C. hirsuta", font = 3, adj = 0.5)
lines(c(x[11], x[15]), c(-0.75, -0.75))
text((x[11] + x[15]) / 2, -0.85, "D. verna", font = 3, adj = 0.5)
legend("topright", legend = colnames(order_tab_fu), fill = col_order_fu,
  bty = "n", inset = c(-0.28, 0), cex = 0.8, border = NA)
dev.off()

order_tab_oo <- tax_prop_table(cdm_vs_oo, tax_vs_oo$order, sample_vs_oo$sample2)
order_tab_oo <- order_tab_oo[c("Atba", "Atrs", "Atrh", "Atro", "Aten", "Chba",
  "Chrs", "Chrh", "Chro", "Chen", "Dvba", "Dvrs", "Dvrh", "Dvro", "Dven"), ]
col_order_oo <- c(alpha(brewer.pal(ncol(order_tab_oo) - 1, "Paired"), 0.75), gray(0.9))
names(col_order_oo) <- colnames(order_tab_oo)
pdf("output/orders_sample_oo.pdf", w = 6, h = 2.5, pointsize = 12)
par(mar = c(6, 4, 1, 6), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, -0.5))
x <- barplot(t(order_tab_oo), space = c(rep(c(0.3, rep(0.1, 4)), 3)), las = 1,
  col = col_order_oo, names.arg = rep(NA, nrow(order_tab_oo)),
  ylab = "Proportion", border = NA)
text(x, y = -0.1, rep(c("Bulk soil", "Root zone", "Rhizosphere", "Root",
  "Endosphere"), 3), xpd = TRUE, srt = 45, adj = 1)
lines(c(x[1], x[5]), c(-0.75, -0.75))
text((x[5] + x[1]) / 2, -0.85, "A. thaliana", font = 3, adj = 0.5)
lines(c(x[6], x[10]), c(-0.75, -0.75))
text((x[6] + x[10]) / 2, -0.85, "C. hirsuta", font = 3, adj = 0.5)
lines(c(x[11], x[15]), c(-0.75, -0.75))
text((x[11] + x[15]) / 2, -0.85, "D. verna", font = 3, adj = 0.5)
legend("topright", legend = colnames(order_tab_oo), fill = col_order_oo,
  bty = "n", inset = c(-0.28, 0), cex = 0.8, border = NA)
dev.off()


# barplots per phylum and sample
phylum_tab_fu <- tax_prop_table(cdm_vs_fu, tax_vs_fu$phylum,
  sample_vs_fu$sample2, n = 10)
phylum_tab_fu <- phylum_tab_fu[c("Atba", "Atrs", "Atrh", "Atro", "Aten", "Chba",
  "Chrs", "Chrh", "Chro", "Chen", "Dvba", "Dvrs", "Dvrh", "Dvro", "Dven"), ]
col_phylum_fu <- c(brewer.pal(ncol(phylum_tab_fu) - 1, "Set3"), gray(0.9))
pdf("output/phylum_sample_fu.pdf", w = 7, h = 3, pointsize = 12)
par(mar = c(6, 4, 1, 6), xpd = T)
x <- barplot(t(phylum_tab_fu), space = c(rep(c(0.3, rep(0.1, 4)), 3)), las = 1,
  col = col_phylum_fu, names.arg = rep(NA, nrow(phylum_tab_fu)),
  ylab = "Proportion", border = NA)
text(x, y = -0.1, rep(c("Bulk soil", "Root zone", "Rhizosphere", "Root",
  "Endosphere"), 3), xpd = TRUE, srt = 45, adj = 1)
lines(c(x[1], x[5]), c(-0.6, -0.6))
text((x[5] + x[1]) / 2, -0.7, "A. thaliana", font = 3, adj = 0.5)
lines(c(x[6], x[10]), c(-0.6, -0.6))
text((x[6] + x[10]) / 2, -0.7, "C. hirsuta", font = 3, adj = 0.5)
lines(c(x[11], x[15]), c(-0.6, -0.6))
text((x[11] + x[15]) / 2, -0.7, "D. verna", font = 3, adj = 0.5)
legend("topright", legend = colnames(phylum_tab_fu), fill = col_phylum_fu,
  bty = "n", inset = c(-0.28, 0), cex = 0.8, border = NA)
dev.off()

### distributon of orders
cdm_order_vs_fu <- t(apply(cdm_vs_fu, 1, function(x) tapply(x, tax_vs_fu$order,
  FUN = sum)))
cdm_order_vs_oo <- t(apply(cdm_vs_oo, 1, function(x) tapply(x, tax_vs_oo$order,
  FUN = sum)))


pdf("output/orders_overall_fu.pdf", h = 6, w = 15, pointsize = 12)
par(mar = c(15, 4, 1, 1))
barplot(sort(colSums(cdm_order_vs_fu), decreasing = T), las = 2)
dev.off()

pdf("output/orders_overall_oo.pdf", h = 6, w = 15, pointsize = 12)
par(mar = c(15, 4, 1, 1))
barplot(sort(colSums(cdm_order_vs_oo), decreasing = T), las = 2)
dev.off()



# number of orders
x <- names(table(tax_vs_fu$order))
length(x)
x <- x[-(grep("^unclassified", x))]
x <- x[-(grep("_Incertae_sedis", x))]
length(x)

x <- names(table(tax_vs_oo$order))
length(x)
x <- x[-(grep("^unclassified", x))]
x <- x[-(grep("_Incertae_sedis", x))]
length(x)

sort(colSums(cdm_order_vs_oo) * 100 / sum(cdm_order_vs_oo), decreasing = T)

### distributon of phyla
cdm_phylum_vs_fu <- t(apply(cdm_vs_fu, 1, function(x) tapply(x, tax_vs_fu$phylum,
  FUN = sum)))
cdm_order_vs_oo <- t(apply(cdm_vs_oo, 1, function(x) tapply(x, tax_oo$order,
  FUN = sum)))

colSums(cdm_phylum_vs_fu) * 100 / sum(cdm_phylum_vs_fu)

sort(table(tax_vs_fu$phylum) * 100 / sum(table(tax_vs_fu$phylum)),
  decreasing = T)

pdf("output/phylum_overall.pdf", h = 6, w = 15, pointsize = 12)
par(mar = c(15, 4, 1, 1))
barplot(sort(colSums(cdm_phylum_vs_fu), decreasing = T), las = 2)
dev.off()



# significance of orders across compartments
order_sig_fu <- as.data.frame(matrix(ncol = 4, nrow = ncol(cdm_order_vs_fu),
  dimnames = list(colnames(cdm_order_vs_fu), c("H", "df", "p", "p.adj"))))
order_sig_fu$H <- apply(cdm_order_vs_fu, 2,
  function(x) kruskal.test(x, sample_fu$compartment)$statistic)
order_sig_fu$df <- apply(cdm_order_vs_fu, 2,
  function(x) kruskal.test(x, sample_fu$compartment)$parameter)
order_sig_fu$p <- apply(cdm_order_vs_fu, 2,
  function(x) kruskal.test(x, sample_fu$compartment)$p.value)
order_sig_fu$p.adj <- p.adjust(order_sig_fu$p, method = "bon")
order_sig_fu[order_sig_fu$p.adj < 0.05, ]

cdm_order_vs_oo <- t(apply(cdm_vs_oo, 1, function(x) tapply(x, tax_oo$order,
  FUN = sum)))
cdm_order_vs_oo <- t(apply(cdm_vs_oo, 1, function(x) tapply(x, tax_oo$order,
  FUN = sum)))
order_sig_oo <- as.data.frame(matrix(ncol = 4, nrow = ncol(cdm_order_vs_oo),
  dimnames = list(colnames(cdm_order_vs_oo), c("H", "df", "p", "p.adj"))))
order_sig_oo$H <- apply(cdm_order_vs_oo, 2,
  function(x) kruskal.test(x, sample_oo$compartment)$statistic)
order_sig_oo$df <- apply(cdm_order_vs_oo, 2,
  function(x) kruskal.test(x, sample_oo$compartment)$parameter)
order_sig_oo$p <- apply(cdm_order_vs_oo, 2,
  function(x) kruskal.test(x, sample_oo$compartment)$p.value)
order_sig_oo$p.adj <- p.adjust(order_sig_oo$p, method = "bon")




### Enrichment analysis



plot(sort(colSums(cdm_vs_fu), decreasing = T), log = "y")

comp <- sample_fu$compartment
comp[comp == "root_soil"] <- "bare_soil"
comp <- droplevels(comp)

mva <- mvabund(cdm_fu[, colSums(cdm_fu) >= 100])
#mva <- mvabund(cdm_fu)
mva_model <- manyglm(mva ~ sample_fu$reads +
  comp, family = "negative.binomial")
saveRDS(mva_model, "output/mva_model_fu.rds")
plot(mva_model)
#mva_anova <- anova(mva_model, show.time = "all", p.uni = "adjusted")
# Time elapsed: 68 hr 2 min 33 sec
#saveRDS(mva_anova, "output/mva_anova.rds")
#mva_anova <- readRDS("output/mva_anova.rds")

mva_coeff <- t(mva_model$coefficients[c("comprhizosphere",
  "comproot", "compendosphere"), ])
colnames(mva_coeff) <- c("rhizosphere", "root", "endosphere")
mva_stder <- t(mva_model$stderr.coefficients[c("comprhizosphere",
  "comproot", "compendosphere"), ])
colnames(mva_stder) <- c("rhizosphere", "root", "endosphere")
mva_cnfin <- mva_stder * 1.96
mva_top <- as.matrix(mva_coeff - mva_cnfin > 0)
mva_low <- as.matrix(mva_coeff + mva_cnfin < 0)
mva_sig <- matrix("gray", nrow = nrow(mva_coeff), ncol = ncol(mva_coeff),
  dimnames = list(rownames(mva_coeff), colnames(mva_coeff)))
mva_sig[mva_top] <- "green"
mva_sig[mva_low] <- "red"
mva_sig <- as.data.frame(mva_sig)
mva_low <- as.data.frame(mva_low)
mva_top <- as.data.frame(mva_top)
for(i in colnames(mva_sig)) mva_sig[, i] <- as.character(mva_sig[, i])

write.table(mva_coeff, "output/mvs_coeff.csv", sep = ";", col.names = NA)
write.table(mva_cnfin, "output/mva_cnfin.csv", sep = ";", col.names = NA)

library(Hmisc)
for(i in rownames(mva_coeff)) {
  pdf(paste0("coefficients/", i, ".pdf"), h = 3, w = 4)
  errbar(1:3, mva_coeff[i, ], mva_coeff[i, ] - mva_cnfin[i, ],
    mva_coeff[i, ] + mva_cnfin[i, ], cap = 0, ylab = "coeff",
    xlab = "compartment")
  mtext(paste(tax_fu[i, "genus"], i), side = 3, line = 1)
  dev.off()
}



mva_coeff <- as.data.frame(mva_coeff)
mva_coeff$abundance <- rep(NA, nrow(mva_coeff))
for(i in rownames(mva_coeff)) {
  mva_coeff[i, "abundance"] <- sum(cdm_vs_fu[, i])
}

pdf("output/differential_abundance_fu.pdf", w = 8, h = 3)
par(mar = c(3, 3, 1, 1), mfrow = c(1, 3))
plot(log10(mva_coeff$abundance), log10(mva_coeff$endosphere + 1000) - 3,
  col = mva_sig$endosphere, ylim = c(-0.05, 0.01))
mtext(sum(mva_sig$endosphere == "green"), side = 3, line = -1.5, col = "green",
  adj = 0.95)
mtext(sum(mva_sig$endosphere == "red"), side = 1, line = -1.5, col = "red",
  adj = 0.95)
lines(c(0, 3), c(2, 2))
plot(log10(mva_coeff$abundance), log10(mva_coeff$root + 1000) - 3,
  col = mva_sig$root, ylim = c(-0.05, 0.01))
mtext(sum(mva_sig$root == "green"), side = 3, line = -1.5, col = "green",
  adj = 0.95)
mtext(sum(mva_sig$root == "red"), side = 1, line = -1.5, col = "red",
  adj = 0.95)
plot(log10(mva_coeff$abundance), log10(mva_coeff$rhizosphere + 1000) - 3,
  col = mva_sig$rhizosphere, ylim = c(-0.05, 0.01))
mtext(sum(mva_sig$rhizosphere == "green"), side = 3, line = -1.5, col = "green",
  adj = 0.95)
mtext(sum(mva_sig$rhizosphere == "red"), side = 1, line = -1.5, col = "red",
  adj = 0.95)
dev.off()


pdf("output/venn_fu.pdf", w = 4, h = 2, pointsize = 12)
par(mfrow = c(1, 2), mar = rep(1, 4))
venn(mva_top)
venn(mva_low)
dev.off()


rownames(mva_top)[mva_top$endosphere & mva_top$root]
rownames(mva_top)[mva_top$endosphere & !(mva_top$root)]


enriched_endo <- rownames(mva_top)[mva_top$endosphere]
enriched_root <- rownames(mva_top)[mva_top$root]
enriched_rhizo <- rownames(mva_top)[mva_top$rhizosphere]
depleted_endo <- rownames(mva_low)[mva_low$endosphere]
depleted_root <- rownames(mva_low)[mva_low$root]
depleted_rhizo <- rownames(mva_low)[mva_low$rhizosphere]

xe <- table(tax_fu[enriched_endo, "order"])
xe <- xe * 100 / sum(xe)
xr <- table(tax_fu[enriched_root, "order"])
xr <- xr * 100 / sum(xr)
xrh <- table(tax_fu[enriched_rhizo, "order"])
xrh <- xrh * 100 / sum(xrh) 
ye <- table(tax_fu[depleted_endo, "order"])
ye <- ye * 100 / sum(ye)
yr <- table(tax_fu[depleted_root, "order"])
yr <- yr * 100 / sum(yr)
yrh <- table(tax_fu[depleted_rhizo, "order"])
yrh <- yrh * 100 / sum(yrh)
z <- cbind(Erh = xrh, Er = xr, Ee = xe, Drh = yrh, Dr = yr, De = ye)

tab <- z[rownames(z) %in% colnames(order_tab_fu), ] 
tab <- rbind(tab, others = colSums(z[!(rownames(z) %in%
  colnames(order_tab_fu)), ]))
tab <- tab[colnames(order_tab_fu), ]

pdf("output/orders_enrich_fu.pdf", w = 2.4, h = 2.1, pointsize = 12)
par(mar = c(4, 4, 1, 1), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, 0))
barplot(tab, las = 1, col = col_order_fu, ylab = "Proportion", border = NA,
  space = c(0.1, 0.1, 0.1, 0.3, 0.1, 0.1))
dev.off()


enriched_plant <- unique(c(enriched_endo, enriched_root))
endo_tab <- cbind(order = as.character(tax_fu[enriched_plant, "order"]),
  ASV = paste(gsub("_", " ", tax_fu[enriched_plant, "species"]), enriched_plant))
x <- cdm_fu[, enriched_plant]
x <- apply(x, 2, function(x) tapply(x, sample_fu$species, sum))
x[x > 0] <- 1
x <- t(x)
endo_tab <- cbind(endo_tab, x)
endo_tab <- as.data.frame(endo_tab)
endo_tab$in_rhizoplane <- as.numeric(enriched_plant %in% enriched_root)
endo_tab$in_endosphere <- as.numeric(enriched_plant %in% enriched_endo)
write.table(endo_tab, "output/endo_tab_fu.csv", col.names = NA, sep = ";")



# manual glm with LR test
x <- glm(cdm_vs_fu[, 2] ~ sample_fu$reads + sample_fu$compartment,
  family = tweedie(var.power = 1.2, link.power = 1.2))
library(MASS)
x <- glm.nb(cdm_fu[, 2] ~ sample_fu$reads + sample_fu$compartment)
plot(x)

summary(x)
x$coefficients


# LR tests removing entire variable
dif_abund <- function(cdm, data) {
  require(MASS)
  data[data$compartment == "root_soil", "compartment"] <- "bare_soil"
  data <- droplevels(data)

  result <- data.frame(matrix(NA, ncol = 9, nrow = ncol(cdm) * 3,
    dimnames = list(1:(ncol(cdm) * 3), c("ASV", "comp", "coef", "logcoef",
    "reads", "logreads", "LR", "P", "Padj"))))
  c <- 0
  for(comp in c("rhizosphere", "root", "endosphere")) {
    data_sub <- droplevels(data[data$compartment == comp |
      data$compartment == "bare_soil", ])
    cdm_sub <- cdm[rownames(data_sub), ]
    for(i in 1:ncol(cdm_sub)) {
      c <- c + 1
      if(c %% 100 == 0) message(paste("process: ", c))
      try(m1 <- glm.nb(cdm_sub[, i] ~ data_sub$reads + data_sub$compartment),
        silent = T)
      #try(m1 <- glm(cdm_sub[, i] ~ data_sub$reads + data_sub$compartment,
      #  family = tweedie()), silent = T)
      try(m0 <- glm.nb(cdm_sub[, i] ~ data_sub$reads), silent = T)
      #try(m0 <- glm(cdm_sub[, i] ~ data_sub$reads, family = tweedie()),
      #  silent = T)
      if(!exists("m1")) m1 <- data.frame(coefficients = c(NA, NA, NA))
      if(exists("m1") & exists("m0")) {
        try(lr <- anova(m0, m1), silent = T)
      } else {
        lr <- data.frame(LR = c(NA, NA), Pr = c(NA, NA))
      }
      if(!exists("lr")) lr <- data.frame(LR = c(NA, NA), Pr = c(NA, NA))
      result[c, ] <- c(colnames(cdm_sub)[i], comp, m1$coefficients[3],
        log10(m1$coefficients[3]), sum(cdm_sub[, i]), log10(sum(cdm_sub[, i])),
        lr$LR[2], lr$Pr[2], NA)
      rm(m1, m0, lr)
    }
  }
  result$Padj <- p.adjust(result$P, method = "holm")
  return(result)   
}
dif_abund_fu <- dif_abund(cdm_fu, sample_fu)
dif_abund_fu$logcoef <- log10(as.numeric(dif_abund_fu$coef) + 1000) - 3

sel <- dif_abund_fu[dif_abund_fu$comp == "root", ]
sel <- na.omit(sel)
col <- rep("gray", nrow(sel))
col[dif_abund_fu$Padj < 0.01] <- "red"
col[is.na(col)] <- "gray"

plot(sel$logcoef ~ sel$logreads, pch = 16, cex = 0.3, col = col)



result[result$Padj < 0.001, ]
x <- na.omit(x)
nrow(x)
nrow(x[x$comp == "endosphere", ])
nrow(x[x$comp == "endosphere" & as.numeric(x$coef) > 0 & as.numeric(x$Padj) < 0.01, ])
nrow(x[x$comp == "endosphere" & as.numeric(x$coef) < 0 & as.numeric(x$Padj) < 0.01, ])
nrow(x[x$comp == "root" & as.numeric(x$coef) > 0 & as.numeric(x$Padj) < 0.01, ])
nrow(x[x$comp == "root" & as.numeric(x$coef) < 0 & as.numeric(x$Padj) < 0.01, ])
nrow(x[x$comp == "rhizosphere" & as.numeric(x$coef) > 0 & as.numeric(x$Padj) < 0.01, ])
nrow(x[x$comp == "rhizosphere" & as.numeric(x$coef) < 0 & as.numeric(x$Padj) < 0.01, ])










# LR tests removing levels of the variable
dif_abund2 <- function(cdm, data) {
  require(MASS)
  data[data$compartment == "root_soil", "compartment"] <- "bare_soil"
  data <- droplevels(data)

  result <- data.frame(matrix(NA, ncol = 9, nrow = ncol(cdm) * 3,
    dimnames = list(1:(ncol(cdm) * 3), c("ASV", "comp", "coef", "logcoef",
    "reads", "logreads", "LR", "P", "Padj"))))

  c <- 0 
  for(i in 1:ncol(cdm)) {
      try(m1 <- glm.nb(cdm[, i] ~ data$reads + data$compartment,
        maxit = 100), silent = T)
    for(comp in c("rhizosphere", "root", "endosphere")) {
      c <- c + 1
      if(c %% 100 == 0) message(paste("process: ", c))
      data_sub <- droplevels(data[data$compartment != comp, ])
      cdm_sub <- cdm[rownames(data_sub), ]
      try(m0 <- glm.nb(cdm_sub[, i] ~ data_sub$reads, maxit = 100), silent = T)
      if(!exists("m1")) m1 <- data.frame(coefficients = c(NA, NA, NA))
      if(exists("m1") & exists("m0")) {
        try(lr <- anova(m0, m1), silent = T)
      } else {
        lr <- data.frame(LR = c(0, 0), Pr = c(1, 1))
      }
      if(!exists("lr")) lr <- data.frame(LR = c(NA, NA), Pr = c(NA, NA))
      result[c, ] <- c(colnames(cdm_sub)[i], comp,
        m1$coefficients[paste0("data$compartment", comp)], NA,
        sum(cdm_sub[, i]), NA, lr$LR[2], lr$Pr[2], NA)
      rm(m0, lr)
    }
  }
  result$Padj <- p.adjust(result$P, method = "holm")
  result$logcoef <- log10(as.numeric(result$coef) + 1000) - 3
  result$logreads <- log10(as.numeric(result$reads))
  return(result)
}


#dif_abund_fu <- dif_abund2(cdm_fu, sample_fu)
write.table(dif_abund_fu, "output/dif_abund_fu.csv", sep = ";", col.names = NA)

sel <- dif_abund_fu[dif_abund_fu$comp == "endosphere", ]
sel <- na.omit(sel)
col <- rep("gray", nrow(sel))
col[dif_abund_fu$Padj < 0.01] <- "red"
col[is.na(col)] <- "gray"

plot(sel$logcoef ~ sel$logreads, pch = 16, cex = 0.5, col = col)
plot(sel$coef ~ sel$reads, pch = 16, cex = 0.5, col = col)


dif_abund_fu[dif_abund_fu$Padj < 0.001, ]
#x <- na.omit(x)
nrow(dif_abund_fu)
nrow(dif_abund_fu[dif_abund_fu$comp == "endosphere", ])
nrow(dif_abund_fu[dif_abund_fu$comp == "endosphere" & as.numeric(dif_abund_fu$coef) > 0, ])
nrow(dif_abund_fu[dif_abund_fu$comp == "root" & as.numeric(dif_abund_fu$coef) > 0, ])
nrow(dif_abund_fu[dif_abund_fu$comp == "endosphere" & as.numeric(dif_abund_fu$coef) < 0, ])
nrow(dif_abund_fu[dif_abund_fu$comp == "root" & as.numeric(dif_abund_fu$coef) < 0, ])



p_cutoff <- 0.01
x <- dif_abund_fu[dif_abund_fu$Padj <= p_cutoff, ]
x_endo_h <- x[x$comp == "endosphere" & x$coef > 0, ]
x_endo_l <- x[x$comp == "endosphere" & x$coef < 0, ]
x_root_h <- x[x$comp == "root" & x$coef > 0, ]
x_root_l <- x[x$comp == "root" & x$coef < 0, ]
x_rhizo_h <- x[x$comp == "rhizosphere" & x$coef > 0, ]
x_rhizo_l <- x[x$comp == "rhizosphere" & x$coef < 0, ]



dif_coeff <- data.frame(
  rhizosphere = as.numeric(dif_abund_fu[dif_abund_fu$comp == "rhizosphere",
  "coef"]), root = as.numeric(dif_abund_fu[dif_abund_fu$comp == "root",
  "coef"]), endosphere = as.numeric(dif_abund_fu[dif_abund_fu$comp ==
  "endosphere", "coef"]))
dif_coeff <- as.matrix(dif_coeff)
rownames(dif_coeff) <- unique(dif_abund_fu$ASV)

dif_pval <- data.frame(
  rhizosphere = as.numeric(dif_abund_fu[dif_abund_fu$comp == "rhizosphere",
  "Padj"]), root = as.numeric(dif_abund_fu[dif_abund_fu$comp == "root",
  "Padj"]), endosphere = as.numeric(dif_abund_fu[dif_abund_fu$comp ==
  "endosphere", "Padj"]))
dif_pval <- as.matrix(dif_pval)
rownames(dif_pval) <- unique(dif_abund_fu$ASV)

dif_top <- dif_coeff > 0 & dif_pval <= p_cutoff
dif_low <- dif_coeff < 0 & dif_pval <= p_cutoff
dif_top[is.na(dif_top)] <- F
dif_low[is.na(dif_low)] <- F

pdf("output/venn_fu2.pdf", w = 4, h = 2, pointsize = 12)
par(mfrow = c(1, 2), mar = rep(1, 4))
venn(dif_top)
venn(dif_low)
dev.off()

dif_top <- as.data.frame(dif_top)
dif_low <- as.data.frame(dif_low)

enriched_endo <- na.omit(rownames(dif_top)[dif_top$endosphere])
enriched_root <- na.omit(rownames(dif_top)[dif_top$root])
enriched_rhizo <- na.omit(rownames(dif_top)[dif_top$rhizosphere])
depleted_endo <- na.omit(rownames(dif_low)[dif_low$endosphere])
depleted_root <- na.omit(rownames(dif_low)[dif_low$root])
depleted_rhizo <- na.omit(rownames(dif_low)[dif_low$rhizosphere])

xe <- table(tax_fu[enriched_endo, "order"])
xe <- xe * 100 / sum(xe)
xr <- table(tax_fu[enriched_root, "order"])
xr <- xr * 100 / sum(xr)
xrh <- table(tax_fu[enriched_rhizo, "order"])
xrh <- xrh * 100 / sum(xrh) 
ye <- table(tax_fu[depleted_endo, "order"])
ye <- ye * 100 / sum(ye)
yr <- table(tax_fu[depleted_root, "order"])
yr <- yr * 100 / sum(yr)
yrh <- table(tax_fu[depleted_rhizo, "order"])
yrh <- yrh * 100 / sum(yrh)
z <- cbind(Erh = xrh, Er = xr, Ee = xe, Drh = yrh, Dr = yr, De = ye)

tab <- z[rownames(z) %in% colnames(order_tab_fu), ] 
tab <- rbind(tab, others = colSums(z[!(rownames(z) %in%
  colnames(order_tab_fu)), ]))
tab <- tab[colnames(order_tab_fu), ]

pdf("output/orders_enrich_fu2.pdf", w = 2.4, h = 2.1, pointsize = 12)
par(mar = c(4, 4, 1, 1), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, 0))
barplot(tab, las = 1, col = col_order_fu, ylab = "Proportion", border = NA,
  space = c(0.1, 0.1, 0.1, 0.3, 0.1, 0.1))
dev.off()


enriched_plant <- unique(c(enriched_endo, enriched_root))
enriched_plant <- na.omit(enriched_plant)
endo_tab <- cbind(order = as.character(tax_fu[enriched_plant, "order"]),
  ASV = paste(gsub("_", " ", tax_fu[enriched_plant, "species"]), enriched_plant))
x <- cdm_fu[, enriched_plant]
x <- apply(x, 2, function(x) tapply(x, sample_fu$species, sum))
x[x > 0] <- 1
x <- t(x)
endo_tab <- cbind(endo_tab, x)
endo_tab <- as.data.frame(endo_tab)
endo_tab$in_rhizoplane <- as.numeric(enriched_plant %in% enriched_root)
endo_tab$in_endosphere <- as.numeric(enriched_plant %in% enriched_endo)
write.table(endo_tab, "output/endo_tab_fu.csv", col.names = NA, sep = ";")































nrow(x[x$comp == "root" & as.numeric(x$Padj) < 0.01 & as.numeric(x$coef) > 0, ])
nrow(x[x$comp == "endosphere" & as.numeric(x$Padj) < 0.01 & as.numeric(x$coef) > 0, ])
nrow(x[x$comp == "rhizosphere" & as.numeric(x$Padj) < 0.01 & as.numeric(x$coef) > 0, ])


root_e <- x[x$comp == "root" & as.numeric(x$Padj) < 0.01 & as.numeric(x$coef) > 0, ]
endo_e <- x[x$comp == "endosphere" & as.numeric(x$Padj) < 0.01 & as.numeric(x$coef) > 0, ]
rhiz_e <- x[x$comp == "rhizosphere" & as.numeric(x$Padj) < 0.01 & as.numeric(x$coef) > 0, ]

root_e$ASV %in% endo_e$ASV
endo_e$ASV %in% rhiz_e$ASV
root_e$ASV %in% rhiz_e$ASV

tax_fu[root_e$ASV, "species"]



# extract p-values from model
dif_abund3 <- function(cdm, data) {
  require(MASS)
  data[data$compartment == "root_soil", "compartment"] <- "bare_soil"
  data <- droplevels(data)

  models <- vector(mode = "list", length = ncol(cdm))
  names(models) <- colnames(cdm)

  for(i in names(models)) {
    try(models[[i]] <- glm.nb(cdm[, i] ~ data$reads + data$compartment,
        maxit = 1000), silent = T)
  }

  result <- data.frame(matrix(NA, ncol = 9, nrow = ncol(cdm) * 3,
    dimnames = list(1:(ncol(cdm) * 3), c("ASV", "comp", "coef", "logcoef",
    "reads", "logreads", "z", "P", "Padj"))))
  result$ASV <- rep(names(models), 3)
  result$comp <- c(rep("rhizosphere", ncol(cdm)), rep("root", ncol(cdm)),
    rep("endosphere", ncol(cdm)))
  result$reads <- rep(colSums(cdm), 3)
  result$logreads <- log10(result$reads)
  for(i in names(models)) {
    result[result$ASV == i & result$comp == "rhizosphere", "coef"] <-
      as.numeric(models[[i]]$coefficients["data$compartmentrhizosphere"])
    result[result$ASV == i & result$comp == "rhizosphere", "z"] <-
      as.numeric(summary(models[[i]])$coefficients[4, 3])
    result[result$ASV == i & result$comp == "rhizosphere", "P"] <-
      as.numeric(summary(models[[i]])$coefficients[4, 4])

    result[result$ASV == i & result$comp == "root", "coef"] <-
      as.numeric(models[[i]]$coefficients["data$compartmentroot"])
    result[result$ASV == i & result$comp == "root", "z"] <-
      as.numeric(summary(models[[i]])$coefficients[5, 3])
    result[result$ASV == i & result$comp == "root", "P"] <-
      as.numeric(summary(models[[i]])$coefficients[5, 4])

    result[result$ASV == i & result$comp == "endosphere", "coef"] <-
      as.numeric(models[[i]]$coefficients["data$compartmentroot"])
    result[result$ASV == i & result$comp == "endosphere", "z"] <-
      as.numeric(summary(models[[i]])$coefficients[3, 3])
    result[result$ASV == i & result$comp == "endosphere", "P"] <-
      as.numeric(summary(models[[i]])$coefficients[3, 4])    
  }




  return(result)
}









### oomycetes

comp <- sample_oo$compartment
comp[comp == "root_soil"] <- "bare_soil"
comp <- droplevels(comp)

mva <- mvabund(cdm_oo[, colSums(cdm_oo) >= 100])
#mva <- mvabund(cdm_oo)
mva_model <- manyglm(mva ~ sample_oo$reads +
  comp, family = "negative.binomial")
saveRDS(mva_model, "output/mva_model_oo.rds")
plot(mva_model)
#mva_anova <- anova(mva_model)

mva_coeff <- t(mva_model$coefficients[c("comprhizosphere",
  "comproot", "compendosphere"), ])
colnames(mva_coeff) <- c("rhizosphere", "root", "endosphere")
mva_stder <- t(mva_model$stderr.coefficients[c("comprhizosphere",
  "comproot", "compendosphere"), ])
colnames(mva_stder) <- c("rhizosphere", "root", "endosphere")
mva_cnfin <- mva_stder / 1.96
mva_top <- as.matrix(mva_coeff - mva_cnfin > 0)
mva_low <- as.matrix(mva_coeff + mva_cnfin < 0)
mva_sig <- matrix("gray", nrow = nrow(mva_coeff), ncol = ncol(mva_coeff),
  dimnames = list(rownames(mva_coeff), colnames(mva_coeff)))
mva_sig[mva_top] <- "green"
mva_sig[mva_low] <- "red"
mva_sig <- as.data.frame(mva_sig)
mva_low <- as.data.frame(mva_low)
mva_top <- as.data.frame(mva_top)
for(i in colnames(mva_sig)) mva_sig[, i] <- as.character(mva_sig[, i])

write.table(mva_coeff, "output/mvs_coeff.csv", sep = ";", col.names = NA)
write.table(mva_cnfin, "output/mva_cnfin.csv", sep = ";", col.names = NA)

library(Hmisc)
for(i in rownames(mva_coeff)) {
  pdf(paste0("coefficients/", i, ".pdf"), h = 3, w = 4)
  errbar(1:4, mva_coeff[i, ], mva_coeff[i, ] - mva_cnfin[i, ],
    mva_coeff[i, ] + mva_cnfin[i, ], cap = 0, ylab = "coeff",
    xlab = "compartment")
  mtext(paste(tax_oo[i, "genus"], i), side = 3, line = 1)
  dev.off()
}



mva_coeff <- as.data.frame(mva_coeff)
mva_coeff$abundance <- rep(NA, nrow(mva_coeff))
for(i in rownames(mva_coeff)) {
  mva_coeff[i, "abundance"] <- sum(cdm_vs_oo[, i])
}



pdf("output/differential_abundance_oo.pdf", w = 8, h = 3)
par(mar = c(3, 3, 1, 1), mfrow = c(1, 3))
plot(log10(mva_coeff$abundance), log10(mva_coeff$endosphere + 1000) - 3,
  col = mva_sig$endosphere, ylim = c(-0.05, 0.01))
mtext(sum(mva_sig$endosphere == "green"), side = 3, line = -1.5, col = "green",
  adj = 0.95)
mtext(sum(mva_sig$endosphere == "red"), side = 1, line = -1.5, col = "red",
  adj = 0.95)
lines(c(0, 3), c(2, 2))
plot(log10(mva_coeff$abundance), log10(mva_coeff$root + 1000) - 3,
  col = mva_sig$root, ylim = c(-0.05, 0.01))
mtext(sum(mva_sig$root == "green"), side = 3, line = -1.5, col = "green",
  adj = 0.95)
mtext(sum(mva_sig$root == "red"), side = 1, line = -1.5, col = "red",
  adj = 0.95)
plot(log10(mva_coeff$abundance), log10(mva_coeff$rhizosphere + 1000) - 3,
  col = mva_sig$rhizosphere, ylim = c(-0.05, 0.01))
mtext(sum(mva_sig$rhizosphere == "green"), side = 3, line = -1.5, col = "green",
  adj = 0.95)
mtext(sum(mva_sig$rhizosphere == "red"), side = 1, line = -1.5, col = "red",
  adj = 0.95)
dev.off()


rownames(mva_top)[mva_top$endosphere & mva_top$root]
rownames(mva_top)[mva_top$endosphere & !(mva_top$root)]

table(droplevels(tax_oo[rownames(mva_top)[mva_top$endosphere], "order"]))
table(droplevels(tax_oo[rownames(mva_top)[mva_top$endosphere & mva_top$root], "order"]))
table(droplevels(tax_oo[rownames(mva_top)[mva_top$endosphere & !(mva_top$root)], "order"]))

library(gplots)
pdf("output/venn_oo.pdf", w = 4, h = 2, pointsize = 12)
par(mfrow = c(1, 2), mar = rep(1, 4))
venn(mva_top)
venn(mva_low)
dev.off()



enriched_endo <- rownames(mva_top)[mva_top$endosphere]
enriched_root <- rownames(mva_top)[mva_top$root]
enriched_rhizo <- rownames(mva_top)[mva_top$rhizosphere]
depleted_endo <- rownames(mva_low)[mva_low$endosphere]
depleted_root <- rownames(mva_low)[mva_low$root]
depleted_rhizo <- rownames(mva_low)[mva_low$rhizosphere]

xe <- table(tax_oo[enriched_endo, "order"])
xe <- xe * 100 / sum(xe)
xr <- table(tax_oo[enriched_root, "order"])
xr <- xr * 100 / sum(xr)
xrh <- table(tax_oo[enriched_rhizo, "order"])
xrh <- xrh * 100 / sum(xrh) 
ye <- table(tax_oo[depleted_endo, "order"])
ye <- ye * 100 / sum(ye)
yr <- table(tax_oo[depleted_root, "order"])
yr <- yr * 100 / sum(yr)
yrh <- table(tax_oo[depleted_rhizo, "order"])
yrh <- yrh * 100 / sum(yrh)
z <- cbind(Erh = xrh, Er = xr, Ee = xe, Drh = yrh, Dr = yr, De = ye)

tab <- z[rownames(z) %in% colnames(order_tab_oo), ] 
tab <- rbind(tab, others = colSums(z[!(rownames(z) %in%
  colnames(order_tab_oo)), ]))
tab <- tab[colnames(order_tab_oo), ]

pdf("output/orders_enrich_oo.pdf", w = 2.4, h = 2.1, pointsize = 12)
par(mar = c(4, 4, 1, 1), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, 0))
barplot(tab, las = 1, col = col_order_oo, ylab = "Proportion", border = NA,
  space = c(0.1, 0.1, 0.1, 0.3, 0.1, 0.1))
dev.off()


enriched_plant <- unique(c(enriched_endo, enriched_root))
endo_tab <- cbind(order = as.character(tax_oo[enriched_plant, "order"]),
  ASV = paste(gsub("_", " ", tax_oo[enriched_plant, "species"]), enriched_plant))
x <- cdm_oo[, enriched_plant]
x <- apply(x, 2, function(x) tapply(x, sample_oo$species, sum))
x[x > 0] <- 1
x <- t(x)
endo_tab <- cbind(endo_tab, x)
endo_tab <- as.data.frame(endo_tab)
endo_tab$in_rhizoplane <- as.numeric(enriched_plant %in% enriched_root)
endo_tab$in_endosphere <- as.numeric(enriched_plant %in% enriched_endo)
write.table(endo_tab, "output/endo_tab_oo.csv", col.names = NA, sep = ";")







################################################################################


### heatmap
x <- cdm_vs_fu
x <- x[, names(sort(colSums(x), decreasing = T))]
x <- decostand(x, "max")

y <- sample_vs_fu
y$sort <- rep(1, nrow(y)) 
for(i in rownames(y)) if(y[i, "compartment"] == "root_soil") y[i, "sort"] <- 2
for(i in rownames(y)) if(y[i, "compartment"] == "rhizosphere") y[i, "sort"] <- 3
for(i in rownames(y)) if(y[i, "compartment"] == "root") y[i, "sort"] <- 4
for(i in rownames(y)) if(y[i, "compartment"] == "endosphere") y[i, "sort"] <- 5
y <- y[order(y$sort),]
x <- x[rownames(y), ]

pdf("output/hm.pdf", w = 20, h = 15)
heatmap.2(x, Rowv = F, Colv = F, trace = "none")
dev.off()


### Indicator OTUs 
# fungi
indssp_fu <- multipatt(as.data.frame(cdm_vs_fu), sample_vs_fu$compartment,
  control = how(nperm = 999))
#indssp_fu <- multipatt(as.data.frame(decostand(cdm_vs_fu, method = "total")),
#  sample_vs_fu$compartment, control = how(nperm = 999))
summary(indssp_fu, alpha = 0.01, indvalcomp = T)
indssp_fu_tab <- na.omit(indssp_fu$sign[indssp_fu$sign$p.value < 0.01, ])
indssp_fu_tab$OTU <- as.vector(tax_vs_fu[rownames(indssp_fu_tab),
   "species"])

g_fu <- indssp_to_network(indssp_fu)
#g_fu <- graph(edges(g_fu))

#node_col <- brewer.pal(5, "RdYlBu")
#node_col <- viridis(5)
node_col <- c("#C83737", "#E9AFAF", "#DECD87", "#AFC6E9", "#3771C8")
node_lab <- names(V(g_fu))
node_col2 <- node_col[degree(g_fu)]
names(node_col2) <- node_lab
#node_col2[grep("^s.", names(node_col2))] <- alpha(col_comp, 0.75)
node_col2[grep("^s.", names(node_col2))] <- alpha("white", 0.5)

node_cex <- colSums(cdm_vs_fu[, colnames(cdm_vs_fu) %in% names(V(g_fu))])
node_cex <- sqrt(node_cex) / 3
node_cex <- c(node_cex, rep(15, 5))
names(node_cex) <- c(node_lab[grep("^FU", node_lab)],
  node_lab[grep("^s.", node_lab)])
node_cex <- node_cex[node_lab]
node_bty <- node_cex
node_bty[grep("^FU", names(node_bty))] <- "black"
node_bty[grep("^s.", names(node_bty))] <- 0

node_lab[grep("^FU", node_lab)] <- ""
node_lab[grep("s.endosphere", node_lab)] <- "En"
node_lab[grep("s.root$", node_lab)] <- "Ro"
node_lab[grep("s.rhizosphere", node_lab)] <- "Rh"
node_lab[grep("s.root_soil", node_lab)] <- "Rs"
node_lab[grep("s.bare_soil", node_lab)] <- "Bs"

pdf("output/network_fu.pdf", h = 6, w = 6)
par(mar = rep(0, 4))
plot(g_fu, vertex.size = node_cex ,
  vertex.color = node_col2, vertex.shape = g_fu$shape,
  vertex.label = node_lab, vertex.label.family = "Helvetica",
  vertex.label.font = 1, vertex.label.cex = 1.5, vertex.frame.color = node_bty,
  vertex.label.color = "black")
legend("bottomleft", legend = rep("", 5), fill = node_col, bty = "n",
  horiz = T, pt.cex = 1.5, border = NA, x.intersp = 0)
dev.off()


na.omit(indssp_fu_tab[indssp_fu_tab$s.root == 1, ])
na.omit(indssp_fu_tab[indssp_fu_tab$s.endosphere == 1, ])

head(tax_vs_fu[tax_vs_fu$order == "Helotiales", ])

x <- decostand(cdm_vs_fu, method = "total")
boxplot(x[, "FU00002"] ~ sample_vs_fu$compartment)
boxplot(x[, "FU00009"] ~ sample_vs_fu$compartment)


boxplot(cdm_vs_fu[, "FU00160"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00294"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00464"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00160"] ~ sample_vs_fu$sample2)
boxplot(cdm_vs_fu[, "FU00294"] ~ sample_vs_fu$sample2)
boxplot(cdm_vs_fu[, "FU00464"] ~ sample_vs_fu$sample2)


boxplot(cdm_vs_fu[, "FU00004"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00012"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00014"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00017"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00043"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00051"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00134"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00147"] ~ sample_vs_fu$compartment)

specnumber(cdm_fu, sample$compartment)
3 * 100 / (248 + 571)

# oomycetes
indssp_oo <- multipatt(as.data.frame(cdm_vs_oo),
  sample_vs_oo$compartment, control = how(nperm = 999))
#indssp_oo <- multipatt(as.data.frame(decostand(cdm_vs_oo, method = "total")),
#  sample_vs_oo$compartment, control = how(nperm = 999))
summary(indssp_oo, alpha = 0.01, indvalcomp = T)
indssp_oo_tab <- na.omit(indssp_oo$sign[indssp_oo$sign$p.value < 0.01, ])
indssp_oo_tab$OTU <- as.vector(tax_vs_oo[rownames(indssp_oo_tab),
   "species"])

g_oo <- indssp_to_network(indssp_oo)
#g_oo <- graph(edges(g_oo))


#node_col <- brewer.pal(5, "RdYlBu")
#node_col <- viridis(5)
node_col <- c("#C83737", "#E9AFAF", "#DECD87", "#AFC6E9", "#3771C8")
node_lab <- names(V(g_oo))
node_col2 <- node_col[degree(g_oo)]
names(node_col2) <- node_lab
#node_col2[grep("^s.", names(node_col2))] <- alpha(col_comp, 0.75)
node_col2[grep("^s.", names(node_col2))] <- alpha("white", 0.5)

node_cex <- colSums(cdm_vs_oo[, colnames(cdm_vs_oo) %in% names(V(g_oo))])
node_cex <- sqrt(node_cex) / 1.5
node_cex <- c(node_cex, rep(15, 5))
names(node_cex) <- c(node_lab[grep("^OO", node_lab)],
  node_lab[grep("^s.", node_lab)])
node_cex <- node_cex[node_lab]
node_bty <- node_cex
node_bty[grep("^OO", names(node_bty))] <- "black"
node_bty[grep("^s.", names(node_bty))] <- 0

node_lab[grep("^OO", node_lab)] <- ""
node_lab[grep("s.endosphere", node_lab)] <- "En"
node_lab[grep("s.root$", node_lab)] <- "Ro"
node_lab[grep("s.rhizosphere", node_lab)] <- "Rh"
node_lab[grep("s.root_soil", node_lab)] <- "Rs"
node_lab[grep("s.bare_soil", node_lab)] <- "Bs"

pdf("output/network_oo.pdf", h = 3, w = 3)
par(mar = rep(0, 4))
plot(g_oo, vertex.size = node_cex ,
  vertex.color = node_col2, vertex.shape = g_oo$shape,
  vertex.label = node_lab, vertex.label.family = "Helvetica",
  vertex.label.font = 1, vertex.label.cex = 1.5, vertex.frame.color = node_bty,
  vertex.label.color = "black")
legend("bottomleft", legend = rep("", 5), fill = node_col, bty = "n",
  horiz = T, pt.cex = 1.5, border = NA, x.intersp = 0)
dev.off()

na.omit(indssp_oo_tab[indssp_oo_tab$s.root == 1, ])
na.omit(indssp_oo_tab[indssp_oo_tab$s.endosphere == 1, ])









### Indicator OTUs root vs. endosphere
# fungi
cdm_vs_fu <- cdm_vs_fu[sample_vs_fu$compartment == "root" |
  sample_vs_fu$compartment == "endosphere", ]
sample_vs_fu <- droplevels(sample_vs_fu[sample_vs_fu$compartment == "root" |
  sample_vs_fu$compartment == "endosphere", ])
indssp_fu <- multipatt(as.data.frame(cdm_vs_fu), sample_vs_fu$compartment,
  control = how(nperm = 999))
#indssp_fu <- multipatt(as.data.frame(decostand(cdm_vs_fu, method = "total")),
#  sample_vs_fu$compartment, control = how(nperm = 999))
summary(indssp_fu, alpha = 0.05, indvalcomp = T)
indssp_fu_tab <- na.omit(indssp_fu$sign[indssp_fu$sign$p.value < 0.01, ])
indssp_fu_tab$OTU <- as.vector(tax_vs_fu[rownames(indssp_fu_tab),
   "species"])

g_fu <- indssp_to_network(indssp_fu)
#g_fu <- graph(edges(g_fu))

#node_col <- brewer.pal(5, "RdYlBu")
#node_col <- viridis(5)
node_col <- c("#C83737", "#E9AFAF", "#DECD87", "#AFC6E9", "#3771C8")
node_lab <- names(V(g_fu))
node_col2 <- node_col[degree(g_fu)]
names(node_col2) <- node_lab
#node_col2[grep("^s.", names(node_col2))] <- alpha(col_comp, 0.75)
node_col2[grep("^s.", names(node_col2))] <- alpha("white", 0.5)

node_cex <- colSums(cdm_vs_fu[, colnames(cdm_vs_fu) %in% names(V(g_fu))])
node_cex <- sqrt(node_cex) / 3
node_cex <- c(node_cex, rep(15, 5))
names(node_cex) <- c(node_lab[grep("^FU", node_lab)],
  node_lab[grep("^s.", node_lab)])
node_cex <- node_cex[node_lab]
node_bty <- node_cex
node_bty[grep("^FU", names(node_bty))] <- "black"
node_bty[grep("^s.", names(node_bty))] <- 0

node_lab[grep("^FU", node_lab)] <- ""
node_lab[grep("s.endosphere", node_lab)] <- "En"
node_lab[grep("s.root$", node_lab)] <- "Ro"
node_lab[grep("s.rhizosphere", node_lab)] <- "Rh"
node_lab[grep("s.root_soil", node_lab)] <- "Rs"
node_lab[grep("s.bare_soil", node_lab)] <- "Bs"

pdf("output/network_fu.pdf", h = 6, w = 6)
par(mar = rep(0, 4))
plot(g_fu, vertex.size = node_cex ,
  vertex.color = node_col2, vertex.shape = g_fu$shape,
  vertex.label = node_lab, vertex.label.family = "Helvetica",
  vertex.label.font = 1, vertex.label.cex = 1.5, vertex.frame.color = node_bty,
  vertex.label.color = "black")
legend("bottomleft", legend = rep("", 5), fill = node_col, bty = "n",
  horiz = T, pt.cex = 1.5, border = NA, x.intersp = 0)
dev.off()


na.omit(indssp_fu_tab[indssp_fu_tab$s.root == 1, ])
na.omit(indssp_fu_tab[indssp_fu_tab$s.endosphere == 1, ])










### indicator OTUs per host

# fungi
cdm_at <- cdm_vs_fu[sample_vs_fu$species == "Arabidopsis thaliana", ]
sam_at <- droplevels(sample_vs_fu[sample_vs_fu$species ==
  "Arabidopsis thaliana", ])
indssp_at_fu <- multipatt(as.data.frame(cdm_at), sam_at$compartment,
  control = how(nperm = 999))
summary(indssp_at_fu, alpha = 0.01, indvalcomp = T)
indssp_fu_at_tab <- na.omit(indssp_at_fu$sign[indssp_at_fu$sign$p.value < 0.01, ])
indssp_fu_at_tab$OTU <- as.vector(tax_vs_fu[rownames(indssp_fu_at_tab),
   "species"])

g_at_fu <- indssp_to_network(indssp_at_fu)
#g_fu <- graph(edges(g_fu))

#node_col <- brewer.pal(5, "RdYlBu")
#node_col <- viridis(5)
node_col <- c("#C83737", "#E9AFAF", "#DECD87", "#AFC6E9", "#3771C8")
node_lab <- names(V(g_at_fu))
node_col2 <- node_col[degree(g_at_fu)]
names(node_col2) <- node_lab
#node_col2[grep("^s.", names(node_col2))] <- alpha(col_comp, 0.75)
node_col2[grep("^s.", names(node_col2))] <- alpha("white", 0.5)

node_cex <- colSums(cdm_at[, colnames(cdm_at) %in% names(V(g_at_fu))])
node_cex <- sqrt(node_cex) / 3
node_cex <- c(node_cex, rep(15, 5))
names(node_cex) <- c(node_lab[grep("^FU", node_lab)],
  node_lab[grep("^s.", node_lab)])
node_cex <- node_cex[node_lab]
node_bty <- node_cex
node_bty[grep("^FU", names(node_bty))] <- "black"
node_bty[grep("^s.", names(node_bty))] <- 0

node_lab[grep("^FU", node_lab)] <- ""
node_lab[grep("s.endosphere", node_lab)] <- "En"
node_lab[grep("s.root$", node_lab)] <- "Ro"
node_lab[grep("s.rhizosphere", node_lab)] <- "Rh"
node_lab[grep("s.root_soil", node_lab)] <- "Rs"
node_lab[grep("s.bare_soil", node_lab)] <- "Bs"

pdf("output/network_fu.pdf", h = 6, w = 6)
par(mar = rep(0, 4))
plot(g_at_fu, vertex.size = node_cex ,
  vertex.color = node_col2, vertex.shape = g_at_fu$shape,
  vertex.label = node_lab, vertex.label.family = "Helvetica",
  vertex.label.font = 1, vertex.label.cex = 1.5, vertex.frame.color = node_bty,
  vertex.label.color = "black")
legend("bottomleft", legend = rep("", 5), fill = node_col, bty = "n",
  horiz = T, pt.cex = 1.5, border = NA, x.intersp = 0)
dev.off()


na.omit(indssp_fu_at_tab[indssp_fu_at_tab$s.root == 1, ])
na.omit(indssp_fu_at_tab[indssp_fu_at_tab$s.endosphere == 1, ])

head(tax_vs_fu[tax_vs_fu$order == "Helotiales", ])

x <- decostand(cdm_vs_fu, method = "total")
boxplot(x[, "FU00002"] ~ sample_vs_fu$compartment)
boxplot(x[, "FU00009"] ~ sample_vs_fu$compartment)


boxplot(cdm_vs_fu[, "FU00160"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00294"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00464"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00160"] ~ sample_vs_fu$sample2)
boxplot(cdm_vs_fu[, "FU00294"] ~ sample_vs_fu$sample2)
boxplot(cdm_vs_fu[, "FU00464"] ~ sample_vs_fu$sample2)


boxplot(cdm_vs_fu[, "FU00004"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00012"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00014"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00017"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00043"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00051"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00134"] ~ sample_vs_fu$compartment)
boxplot(cdm_vs_fu[, "FU00147"] ~ sample_vs_fu$compartment)

specnumber(cdm_fu, sample$compartment)
3 * 100 / (248 + 571)

# oomycetes
indssp_oo <- multipatt(as.data.frame(cdm_vs_oo),
  sample_vs_oo$compartment, control = how(nperm = 999))
#indssp_oo <- multipatt(as.data.frame(decostand(cdm_vs_oo, method = "total")),
#  sample_vs_oo$compartment, control = how(nperm = 999))
summary(indssp_oo, alpha = 0.01, indvalcomp = T)
indssp_oo_tab <- na.omit(indssp_oo$sign[indssp_oo$sign$p.value < 0.01, ])
indssp_oo_tab$OTU <- as.vector(tax_vs_oo[rownames(indssp_oo_tab),
   "species"])

g_oo <- indssp_to_network(indssp_oo)
#g_oo <- graph(edges(g_oo))


#node_col <- brewer.pal(5, "RdYlBu")
#node_col <- viridis(5)
node_col <- c("#C83737", "#E9AFAF", "#DECD87", "#AFC6E9", "#3771C8")
node_lab <- names(V(g_oo))
node_col2 <- node_col[degree(g_oo)]
names(node_col2) <- node_lab
#node_col2[grep("^s.", names(node_col2))] <- alpha(col_comp, 0.75)
node_col2[grep("^s.", names(node_col2))] <- alpha("white", 0.5)

node_cex <- colSums(cdm_vs_oo[, colnames(cdm_vs_oo) %in% names(V(g_oo))])
node_cex <- sqrt(node_cex) / 1.5
node_cex <- c(node_cex, rep(15, 5))
names(node_cex) <- c(node_lab[grep("^OO", node_lab)],
  node_lab[grep("^s.", node_lab)])
node_cex <- node_cex[node_lab]
node_bty <- node_cex
node_bty[grep("^OO", names(node_bty))] <- "black"
node_bty[grep("^s.", names(node_bty))] <- 0

node_lab[grep("^OO", node_lab)] <- ""
node_lab[grep("s.endosphere", node_lab)] <- "En"
node_lab[grep("s.root$", node_lab)] <- "Ro"
node_lab[grep("s.rhizosphere", node_lab)] <- "Rh"
node_lab[grep("s.root_soil", node_lab)] <- "Rs"
node_lab[grep("s.bare_soil", node_lab)] <- "Bs"

pdf("output/network_oo.pdf", h = 3, w = 3)
par(mar = rep(0, 4))
plot(g_oo, vertex.size = node_cex ,
  vertex.color = node_col2, vertex.shape = g_oo$shape,
  vertex.label = node_lab, vertex.label.family = "Helvetica",
  vertex.label.font = 1, vertex.label.cex = 1.5, vertex.frame.color = node_bty,
  vertex.label.color = "black")
legend("bottomleft", legend = rep("", 5), fill = node_col, bty = "n",
  horiz = T, pt.cex = 1.5, border = NA, x.intersp = 0)
dev.off()

na.omit(indssp_oo_tab[indssp_oo_tab$s.root == 1, ])
na.omit(indssp_oo_tab[indssp_oo_tab$s.endosphere == 1, ])












dif_abund_oo <- dif_abund(cdm_oo, sample_oo)
write.table(dif_abund_oo, "output/dif_abund_oo.csv", sep = ";", col.names = NA)


dif_abund_oo[dif_abund_oo$Padj < 0.001, ]
#x <- na.omit(x)
nrow(dif_abund_oo)
nrow(dif_abund_oo[dif_abund_oo$comp == "endosphere", ])
nrow(dif_abund_oo[dif_abund_oo$comp == "endosphere" & as.numeric(dif_abund_oo$coef) > 0, ])
nrow(dif_abund_oo[dif_abund_oo$comp == "root" & as.numeric(dif_abund_oo$coef) > 0, ])
nrow(dif_abund_oo[dif_abund_oo$comp == "endosphere" & as.numeric(dif_abund_oo$coef) < 0, ])
nrow(dif_abund_oo[dif_abund_oo$comp == "root" & as.numeric(dif_abund_oo$coef) < 0, ])



p_cutoff <- 0.01
x <- dif_abund_oo[dif_abund_oo$Padj <= p_cutoff, ]
x_endo_h <- x[x$comp == "endosphere" & x$coef > 0, ]
x_endo_l <- x[x$comp == "endosphere" & x$coef < 0, ]
x_root_h <- x[x$comp == "root" & x$coef > 0, ]
x_root_l <- x[x$comp == "root" & x$coef < 0, ]
x_rhizo_h <- x[x$comp == "rhizosphere" & x$coef > 0, ]
x_rhizo_l <- x[x$comp == "rhizosphere" & x$coef < 0, ]



dif_coeff <- data.frame(
  rhizosphere = as.numeric(dif_abund_oo[dif_abund_oo$comp == "rhizosphere",
  "coef"]), root = as.numeric(dif_abund_oo[dif_abund_oo$comp == "root",
  "coef"]), endosphere = as.numeric(dif_abund_oo[dif_abund_oo$comp ==
  "endosphere", "coef"]))
dif_coeff <- as.matrix(dif_coeff)
rownames(dif_coeff) <- unique(dif_abund_oo$ASV)

dif_pval <- data.frame(
  rhizosphere = as.numeric(dif_abund_oo[dif_abund_oo$comp == "rhizosphere",
  "Padj"]), root = as.numeric(dif_abund_oo[dif_abund_oo$comp == "root",
  "Padj"]), endosphere = as.numeric(dif_abund_oo[dif_abund_oo$comp ==
  "endosphere", "Padj"]))
dif_pval <- as.matrix(dif_pval)
rownames(dif_pval) <- unique(dif_abund_oo$ASV)

dif_top <- dif_coeff > 0 & dif_pval <= p_cutoff
dif_low <- dif_coeff < 0 & dif_pval <= p_cutoff
dif_top[is.na(dif_top)] <- F
dif_low[is.na(dif_low)] <- F


pdf("output/venn_oo2.pdf", w = 4, h = 2, pointsize = 12)
par(mfrow = c(1, 2), mar = rep(1, 4))
venn(dif_top)
venn(dif_low)
dev.off()


dif_top <- as.data.frame(dif_top)
dif_low <- as.data.frame(dif_low)

enriched_endo <- rownames(dif_top)[dif_top$endosphere]
enriched_root <- rownames(dif_top)[dif_top$root]
enriched_rhizo <- rownames(dif_top)[dif_top$rhizosphere]
depleted_endo <- rownames(dif_low)[dif_low$endosphere]
depleted_root <- rownames(dif_low)[dif_low$root]
depleted_rhizo <- rownames(dif_low)[dif_low$rhizosphere]

xe <- table(tax_oo[enriched_endo, "order"])
xe <- xe * 100 / sum(xe)
xr <- table(tax_oo[enriched_root, "order"])
xr <- xr * 100 / sum(xr)
xrh <- table(tax_oo[enriched_rhizo, "order"])
xrh <- xrh * 100 / sum(xrh) 
ye <- table(tax_oo[depleted_endo, "order"])
ye <- ye * 100 / sum(ye)
yr <- table(tax_oo[depleted_root, "order"])
yr <- yr * 100 / sum(yr)
yrh <- table(tax_oo[depleted_rhizo, "order"])
yrh <- yrh * 100 / sum(yrh)
z <- cbind(Erh = xrh, Er = xr, Ee = xe, Drh = yrh, Dr = yr, De = ye)

tab <- z[rownames(z) %in% colnames(order_tab_oo), ] 
tab <- rbind(tab, others = colSums(z[!(rownames(z) %in%
  colnames(order_tab_oo)), ]))
tab <- tab[colnames(order_tab_oo), ]

pdf("output/orders_enrich_oo2.pdf", w = 2.4, h = 2.1, pointsize = 12)
par(mar = c(4, 4, 1, 1), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, 0))
barplot(tab, las = 1, col = col_order_oo, ylab = "Proportion", border = NA,
  space = c(0.1, 0.1, 0.1, 0.3, 0.1, 0.1))
dev.off()


enriched_plant <- unique(c(enriched_endo, enriched_root))
endo_tab <- cbind(order = as.character(tax_oo[enriched_plant, "order"]),
  ASV = paste(gsub("_", " ", tax_oo[enriched_plant, "species"]), enriched_plant))
x <- cdm_oo[, enriched_plant]
x <- apply(x, 2, function(x) tapply(x, sample_oo$species, sum))
x[x > 0] <- 1
x <- t(x)
endo_tab <- cbind(endo_tab, x)
endo_tab <- as.data.frame(endo_tab)
endo_tab$in_rhizoplane <- as.numeric(enriched_plant %in% enriched_root)
endo_tab$in_endosphere <- as.numeric(enriched_plant %in% enriched_endo)
write.table(endo_tab, "output/endo_tab_oo.csv", col.names = NA, sep = ";")






#### supplementary material
# based on abundance
tax_prop_table(cdm_fu, tax_fu$order, n = 15)


x <- colSums(tax_prop_table(cdm_vs_fu, tax_vs_fu$order, n = 15))
x <- sort(x, decreasing = T)
x <- x[c(2:length(x), 1)]
barplot(x)



x <- colSums(tax_prop_table(cdm_vs_oo, tax_vs_oo$order, n = 5))
x <- sort(x, decreasing = T)
x <- x[c(2:length(x), 1)]
barplot(x)


# based on number of OTUs
sort(table(table_fu$order), decreasing = T)

















### end
