### color key by plant species for plotting ------------------------------------
col_spec <- function(x = NULL) {
  require(RColorBrewer)
  col_spec <- brewer.pal(length(levels(sample$sp_code)), "Set1")
  names(col_spec) <- levels(sample$sp_code)
  sample$col_spec <<- col_spec[sample$sp_code]
  if(!is.null(x)) {
    x <- as.character(x)
    col_spec <- col_spec[x]
  }
  return(col_spec)
}


### symbol key by plant species for plotting -----------------------------------
pch_spec <- function(x = NULL) {
  pch_spec <- unique(as.numeric(sample$sp_code)) + 20
  names(pch_spec) <- levels(sample$sp_code)
  sample$pch_spec <<- pch_spec[sample$sp_code]
  if(!is.null(x)) {
    x <- as.character(x)
    pch_spec <- pch_spec[x]
  }
  return(pch_spec)
}


### color key by compartment ---------------------------------------------------
col_comp <- function(x = NULL) {
  require(RColorBrewer)
  col_comp <- brewer.pal(length(levels(sample$compartment)), "Set2")
  names(col_comp) <- levels(sample$compartment)
  sample$col_comp <<- col_comp[sample$compartment]
  if(!is.null(x)) {
    x <- as.character(x)
    col_comp <- col_comp[x]
  }
  return(col_comp)
}


### function to generate distinct colors ---------------------------------------
# adapted from https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color <- function(n, start = 1, random = F, last_gray = T) {
  require(RColorBrewer)
  
  # 433 colors
  # x <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  # 74 colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  x <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
    rownames(qual_col_pals)))

  result <- x[start:(start + n - 1)]
  if(random) result <- sample(x, n)
  if(last_gray) result[length(result)] <- gray(0.9)
  
  return(result)  
}


### boxplot with data points ---------------------------------------------------
boxplot_pt <- function(x, y, log = "", lty = 1, staplewex = 0,
  border = par("fg"), bp_col = NULL, jitter_amount = 0.25, pch = 21,
  pch_col = 1, pch_bg = "transparent", las = 1, xlab = NA, lwd = 1, main = NA,
  ylab = NA, notch = F, y_ticks = 5, y_range_mult = c(0.8, 1.2)) {
  y_range <- signif(c(min(y) * y_range_mult[1], max(y) * y_range_mult[2]), 0)
  bp <- boxplot(y ~ x, log = log, lty = lty, staplewex = staplewex, outline = F,
    axes = F, main = main, xlab = xlab, ylab = ylab, notch = notch,
    ylim = y_range, border = border, col = bp_col, lwd = lwd)
  points(y ~ jitter(as.numeric(x), amount = jitter_amount), pch = pch,
    col = pch_col, bg = pch_bg)
  axis(1, pos = y_range[1], at = c(0.5, length(bp$names) + 0.5), lwd.ticks = 0,
    labels = F, lwd = lwd)
  axis(1, pos = y_range[1], at = 1:length(bp$names), labels = bp$names,
    tick = F, las = las, lwd = lwd)
  axis(1, pos = y_range[1], at = seq(0.5, length(bp$names) + 1, 1),
    lwd.ticks = lwd, lwd = 0, labels = F)
  axis(2, pos = 0.5, at = y_range, labels = F, lwd.ticks = 0, lwd = lwd)
  if (log == "y") {
    require(emdbook)
    axis(2, pos = 0.5, at =  round(lseq(y_range[1], y_range[2], y_ticks), 0),
      labels = T, las = las, lwd.ticks = lwd, lwd = lwd)
  } else {
    axis(2, pos = 0.5, at = seq(y_range[1], y_range[2], length.out = y_ticks),
      labels = T, lwd.ticks = lwd, lwd = lwd)
  }
}


### calculate richness and diversity -------------------------------------------
calculate_diversity <- function(cdm, sam) {
  require(vegan)
  sam$Sobs <- apply(cdm, 1, function(x) sum(x > 0))
  sam$Shan <- exp(diversity(cdm, index = "shannon"))
  return(sam)
}


### modify names in taxonomy ---------------------------------------------------
modify_taxonomy <- function(x) {
  ncat <- ncol(x)
  for(i in rownames(x)) {
    for(j in 2:(ncat - 1)) {
      if(any(grep("unclassified", x[i, j])) &
        any(grep("unclassified", x[i, j - 1], invert = T))) {
        tag <- paste0("unclassified_", x[i, j - 1])
        if(!(tag %in% levels(x[, j]))) {
          levels(x[, j]) <- c(levels(x[, j]), tag)
          x[i, j] <- tag
        } else {
          x[i, j] <- tag
        }
      } else if(any(grep("unclassified", x[i, j])) &
        any(grep("unclassified", x[i, j - 1]))) {
        tag <- as.character(x[i, j - 1])
        if(!(tag %in% levels(x[, j]))) {
          levels(x[, j]) <- c(levels(x[, j]), tag)
          x[i, j] <- tag
        } else {
          x[i, j] <- tag
        }      
      }
    }
    if(any(grep("unclassified", x[i, ncat - 1]))) {
      tag <- as.character(x[i, ncat - 1])
      if(!(tag %in% levels(x[, ncat]))) {
        levels(x[, ncat]) <- c(levels(x[, ncat]), tag)
        x[i, ncat] <- tag
      } else {
        x[i, ncat] <- tag
      }
    } else if(any(grep("_sp", x[i, ncat], invert = T))) {
      tag <- paste0(x[i, ncat - 1], "_sp")
      if(!(tag %in% levels(x[, ncat]))) {
        levels(x[, ncat]) <- c(levels(x[, ncat]), tag)
        x[i, ncat] <- tag
      }
    }
  }
  return(x)
}


### plot rarefaction curves ----------------------------------------------------
plot_rarefaction_curves <- function(cdm, sam,
  outfile = "rarefaction_curves.pdf") {
  require(vegan)
  options(warn = -1)
  pdf(paste0("output/", outfile), w = 9, h = 9, pointsize = 12)
  par(mfrow = c(3, 2), mar = c(4, 4, 1, 1), las = 1, cex = 1.25, lwd = 1.5)
  rarecurve(t(colSums(cdm)), step = 100, label = F, ylab = "Number of OTUs",
    xlab = "Reads", axes = F, main = "overall")
  axis(1, pos = 0, lwd = 1.5)
  axis(2, pos = 0, lwd = 1.5)
  legend("bottomright", levels(sam$sp_code), col = col_spec(), lty = 1,
    lwd = 1.5, bty = "n")
  for(i in levels(sam$compartment)[c(1, 5, 3, 4, 2)]) {
    rarecurve(cdm[sam$compartment == i, ], step = 100,
      label = F, ylab = "Number of OTUs", xlab = "Reads", axes = F,
      ylim = c(0, max(apply(cdm[sam$compartment == i, ],
        1, function(x) sum(x > 0))) * 1.05),
      xlim = c(0, max(rowSums(cdm[sam$compartment == i, ])) *
        1.05), col = col_spec(sam$sp_code[sam$compartment == i]),
        main = i)
    axis(1, pos = 0, lwd = 1.5)
    axis(2, pos = 0, lwd = 1.5)
  }
  dev.off()
}


### boxplots with sample depth -------------------------------------------------
plot_sample_depth <- function(outfile = "reads_x_treatment.pdf") {
  pdf(paste0("output/", outfile), w = 8, h = 5, pointsize = 12)
  layout(matrix(c(rep(1, 5), rep(2, 4), rep(3, 3),
    rep(4, 5), rep(5, 4), rep(6, 3)), byrow = T, nrow = 2))
  boxplot_pt(sample_fu$compartment, sample_fu$reads, log = "y", las = 1,
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
    main = "Oomycota", xlab = "Compartment", ylab = "Reads [log(x)]",
    lwd = 1.25)
  boxplot_pt(sample_oo$species, sample_oo$reads + 1, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
    main = "Oomycota", xlab = "Host species", ylab = "Reads [log(x)]",
    lwd = 1.25)
  boxplot_pt(sample_oo$date, sample_oo$reads + 1, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
    main = "Oomycota", xlab = "Date", ylab = "Reads [log(x)]", lwd = 1.25)
  dev.off()
}


### boxplots with richness/diversity per sample --------------------------------
diversity_boxplots <- function(div = F, outfile = "richness_x_treatment.pdf") {
  if(!div) {
    y_lab <- "Richness [log(S)]"
  } else {
    y_lab <- "Shannon diversity [log(ES)]"
  }
  pdf(paste0("output/", outfile), w = 8, h = 5, pointsize = 12)
  layout(matrix(c(rep(1, 5), rep(2, 4), rep(3, 3),
    rep(4, 5), rep(5, 4), rep(6, 3)), byrow = T, nrow = 2))
  boxplot_pt(sample_fu$compartment, sample_fu$Sobs, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
    xlab = "Compartment", ylab = y_lab, lwd = 1.25)
  boxplot_pt(sample_fu$species, sample_fu$Sobs, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
    xlab = "Host species", ylab = y_lab, lwd = 1.25)
  boxplot_pt(sample_fu$date, sample_fu$Sobs, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2), main = "Fungi",
    xlab = "Date", ylab = y_lab, lwd = 1.25)
  boxplot_pt(sample_oo$compartment, sample_oo$Sobs + 1, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
    main = "Oomycota", xlab = "Compartment",
    ylab = y_lab, lwd = 1.25)
  boxplot_pt(sample_oo$species, sample_oo$Sobs + 1, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
    main = "Oomycota", xlab = "Host species",
    ylab = y_lab, lwd = 1.25)
  boxplot_pt(sample_oo$date, sample_oo$Sobs + 1, log = "y", las = 1,
    border = gray(0.3), pch_col = gray(0.2), pch_bg = gray(0.2),
    main = "Oomycota", xlab = "Date",
    ylab = y_lab, lwd = 1.25)
  dev.off()
}


### fancy plots for diversity data ---------------------------------------------
fancyplot_diversity <- function(sp, var, fac) {
  require(Hmisc)
  x <- fac[fac$sp_code == sp, ]
  t0 <- x$date == "February"
  t1 <- x$date == "April"
  plot(0, 0, type = "n", xlim = c(0.9, 2.1), ylim = c(0, max(fac[, var])),
    axes = F, xlab = NA, ylab = NA, main = sp)
  for(i in levels(x$compartment)) {
    y <- tapply(x[x$compartment == i, var],
      x[x$compartment == i, "date"], mean)[c(2, 1)]
    ysd <- tapply(x[x$compartment == i, var],
      x[x$compartment == i, "date"], sd)[c(2, 1)]
    polygon(c(1, 2, 2, 1, 1), c(y[1] - ysd[1], y[2] - ysd[2],
      y[2] + ysd[2], y[1] + ysd[1], y[1] - ysd[1]),
      col = alpha(col_comp(i), 0.25), border = F)
    lines(c(1, 2), y, type = "o", pch = 16, col = col_comp(i),
      cex = 1.2, lwd = 2)
  }
  axis(1, at = c(1, 2), labels = c("Feb.", "Ap."), lwd = 1.5)
  axis(2, lwd = 1.5)
}


### diversity comparison plots -------------------------------------------------
plot_diversity <- function(div = F, outfile = "richness_comparison.pdf") {
  if(!div) {
    var <- "Sobs"
  } else {
    var <- "Shan"
  }
  pdf(paste0("output/", outfile), w = 7, h = 9, pointsize = 12)
  layout(matrix(c(1:6, 7, 7, 7), ncol = 3, byrow = T))
  par(mar = c(4, 3.5, 1, 0), cex = 1.25, las = 1)
  fancyplot_diversity("At", var, sample_fu)
  mtext("Fungi", side = 2, line = 2.5, las = 0, cex = 1.25)
  fancyplot_diversity("Ch", var, sample_fu)
  fancyplot_diversity("Dv", var, sample_fu)
  fancyplot_diversity("At", var, sample_oo)
  mtext("Oomycota", side = 2, line = 2.5, las = 0, cex = 1.25)
  fancyplot_diversity("Ch", var, sample_oo)
  fancyplot_diversity("Dv", var, sample_oo)
  plot(0, 0, type = "n", axes = F, xlab = NA, ylab = NA)
  legend("top", levels(sample$compartment), pch = 16, lty = 1, col = col_comp(),
    bty = "n", pt.cex = 1.25, lwd = 2, ncol = 3, inset = -0.1)
  dev.off()
}


### remove rare ASVs -----------------------------------------------------------
remove_rare_ASVs <- function(cdm, n_min = 5) {
  cdm <- cdm[, colSums(cdm) >= n_min]
  cdm <- cdm[rowSums(cdm) > 0, ]
  return(cdm)
}


### summary tables of sampling depth and diversity -----------------------------
summarize_diversity <- function(cdm, sam, outfile = "diversity.csv") {
  require(vegan)
  tab_fu <- data.frame(
    n = as.character(table(sam$sample)),
    reads_t = rowSums(apply(cdm, 2,
      function(x) tapply(x, sam$sample, sum))),
    sobs_t = specnumber(cdm, sam$sample),
    reads = paste(round(tapply(sam$reads, sam$sample, mean), 1), "±",
      round(tapply(sam$reads, sam$sample, sd), 1)),
    sobs_s = paste(round(tapply(sam$Sobs, sam$sample, mean), 1), "±",
      round(tapply(sam$Sobs, sam$sample, sd), 1))
  )
  rownames(tab_fu) <- names(tapply(sam$reads, sam$sample, mean))
  write.table(tab_fu, paste0("output/", outfile), sep = ";", col.names = NA)
}


### account for missinf samples in 'sample_oo', for paired comparisons ---------
pair_oo_samples <- function() {
  x <- rownames(sample_fu)[!(rownames(sample_fu) %in% rownames(sample_oo))]
  miss <- data.frame(reads = rep(0, length(x)), Sobs = rep(0, length(x)),
    Shan = rep(0, length(x)), date = c(rep("April", 5), rep("February", 4)))
  rownames(miss) <- x
  x <- rbind(sample_oo[, c("reads", "Sobs", "Shan", "date")], miss)
  x <- x[order(x$date, rownames(x)),]
  return(x)
}


### estimate mean-variance relationship ----------------------------------------
mean_var_relationship <- function(cdm) {
  mean_var <- data.frame(mean = apply(cdm, 2, mean),
    var = apply(cdm, 2, function(x) sd(x) * sd(x)))
  plot(mean_var$var ~ mean_var$mean, log = "xy", xlab = "mean",
    ylab = "variance")
  abline(a = 0, b = 1, lty = 2)
}


### normalization of reads abundance with DSeq ---------------------------------
DSeq_normalization <- function(cdm, sam, tax) {
  require(DESeq)
  dds <- newCountDataSet(t(cdm) + 1, conditions = sam,
    featureData = AnnotatedDataFrame(tax)) # x + 1 to avoid errors 
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  dds <- estimateDispersions(dds, method = "blind", sharingMode = "maximum",
    fitType = "local")
  # change sharingMode to gene-est-only if n > 7
  dispcol <- grep("disp\\_", colnames(fData(dds)))
  if (any(!is.finite(fData(dds)[, dispcol]))) {
    fData(dds)[which(!is.finite(fData(dds)[, dispcol])), dispcol] <- 0
  }
  cdm_vs <- t(exprs(varianceStabilizingTransformation(dds)))
  # set negative values (corresponding to zeroes in raw data set) to zero
  cdm_vs[cdm_vs < 0.0] <- 0.0
  cdm_vs <- cdm_vs[, colSums(cdm_vs) > 0]
  cdm_vs <- cdm_vs[rowSums(cdm_vs) > 0, ]
  return(cdm_vs)
}


### calculate ecological distances per compartment -----------------------------
distance_per_compartment <- function(cdm, sam) {
  comp_dists <- vector("list", length(levels(sam$compartment)))
  names(comp_dists) <- levels(sam$compartment)
  for(i in names(comp_dists)) {
    comp_dists[[i]] <- vegdist(cdm[sam$compartment == i, ], method = "bray")
  }
  return(comp_dists)
}


### calculate PCoAs ------------------------------------------------------------
calculate_PCoA <- function(dist, outfile = "pcoa.pdf") {
  require(vegan)
  pcoa <- cmdscale(dist, eig = T, add = T)
  pcoa$var_explained <- round(eigenvals(pcoa) / sum(eigenvals(pcoa)) * 100, 1)
  pdf(paste0("output/", outfile), w = 3, h = 3)
  par(mar = c(4, 4, 1, 1), las = 1)
  plot_pcoa(pcoa, pcoa$var_explained)
  dev.off()
  return(pcoa)
}


### calculate PCoAs per compartment --------------------------------------------
calculate_PCoAs_per_compartment <- function(dists, outfile_pattern = "pcoa") {
  require(vegan)
  for(i in names(dists)) {
    name <- paste0("pcoa_", i, "_fu")
    pcoa <- cmdscale(dists[[i]], eig = T, add = T)
    expl <- round(eigenvals(pcoa) / sum(eigenvals(pcoa)) * 100, 1)
    pdf(paste0("output/", outfile_pattern, "_", i, ".pdf"), w = 3, h = 3)
    par(mar = c(4, 4, 1, 1), las = 1)
    plot_pcoa(pcoa, expl)
    dev.off()
    assign(name, pcoa, envir = globalenv())
    assign(paste(name, "expl", sep = "_"), expl, envir = globalenv())
  }
}


### calculate dbRDAs -----------------------------------------------------------
calculate_dbRDA <- function(dist, cdm, sam, outfile = "dbrda.pdf") {
  require(vegan)
  dbrda_out <- dbrda(dist ~ sam$sp_code + sam$compartment + sam$date,
    comm = cdm)
  dbrda_out$expl <- round(eigenvals(dbrda_out) /
    sum(eigenvals(dbrda_out)) * 100, 1)
  dbrda_out$aov <- anova(dbrda_out)
  dbrda_out$aov_axis <- anova(dbrda_out, by = "axis")
  dbrda_out$aov_margin <- anova(dbrda_out, by = "margin")
  dbrda_out$aov_var_expl <- (dbrda_out$aov_margin$SumOfSqs /
    with(dbrda_out, tot.chi)) * 100
  names(dbrda_out$aov_var_expl) <- c("Species", "Compartment", "Date",
    "Residuals")
  pdf(paste0("output/", outfile), w = 3, h = 3)
  par(mar = c(4, 4, 1, 1), las = 1)
  plot_cap(dbrda_out, dbrda_out$expl)
  dev.off()

  return(dbrda_out)
}


### calculate dbRDAs per compartment --------------------------------------------
calculate_dbRDAs_per_compartment <- function(dists, cdm, sam,
  outfile_pattern = "dbrda") {
  require(vegan)
  for(i in names(dists)) {
    s  <- droplevels(sam[sam$compartment == i, ])
    d  <- cdm[rownames(s), ]
    dd <- dists[[i]]
    cp <- dbrda(dd ~ s$sp_code + s$date, comm = d)
    cp$expl <- round(eigenvals(cp) / sum(eigenvals(cp)) * 100, 1)
    cp$aov <- anova(cp)
    cp$aov_axis <- anova(cp, by = "axis")
    cp$aov_margin <- anova(cp, by = "margin")
    if(is.null(cp$CA$imaginary.chi)) {
      cp$aov_var_expl <- (cp$aov_margin$SumOfSqs / 
        with(cp, tot.chi)) * 100
    } else {
      cp$aov_var_expl <- (cp$aov_margin$SumOfSqs / 
        with(cp, tot.chi - CA$imaginary.chi)) * 100
    }
    names(cp$aov_var_expl) <- c("Species", "Date", "Residuals")
    assign(paste(outfile_pattern, i, sep = "_"), cp, envir = globalenv())
    pdf(paste0("output/", outfile_pattern, "_", i, ".pdf"), w = 8, h = 8,
      pointsize = 36)
    par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
    plot_cap(cp, cp$expl)
    dev.off()
  }
}


### plot PCoAs -----------------------------------------------------------------
plot_pcoa <- function(pcoa, expl, index = sample) {
  require(vegan)
  fac <- droplevels(index[rownames(pcoa$points), ])
  date_fill <- col_comp(as.character(fac$compartment))
  date_fill[fac$date == "February"] <- 0
  plot(scores(pcoa), type = "n",
    xlab = paste("PCo1 (", expl[1], "%)", sep = ""),
    ylab = paste("PCo2 (", expl[2], "%)", sep = ""))
  lines(c(0, 0), c(-1, 1), col = gray(0.5), lty = 3)
  lines(c(-1, 1), c(0, 0), col = gray(0.5), lty = 3)
  points(scores(pcoa),
    pch = as.numeric(fac$species) + 20,
    col = col_comp(as.character(fac$compartment)),
    bg = date_fill)
}


### plot dbRDAs ----------------------------------------------------------------
plot_cap <- function(cap, expl, index = sample) {
  fac <- droplevels(index[rownames(cap$CCA$u), ])
  date_fill <- col_comp(as.character(fac$compartment))
  date_fill[fac$date == "February"] <- 0
  plot(scores(cap)$sites,
    type = "n",
    xlab = paste("Axis 1 (", expl[1], "%)", sep = ""),
    ylab = paste("Axis 2 (", expl[2], "%)", sep = ""))
  lines(c(0, 0), c(-10, 10), col = gray(0.5), lty = 3)
  lines(c(-10, 10), c(0, 0), col = gray(0.5), lty = 3)
  points(scores(cap)$sites,
    pch = as.numeric(fac$species) + 20,
    col = col_comp(as.character(fac$compartment)),
    bg = date_fill)
}


### plot dbRDA variance --------------------------------------------------------
plot_dbRDA_var <- function(dbrda_obj, outfile = "dbrda_var.pdf") {
  pdf(paste0("output/", outfile), w = 6, h = 8, pointsize = 36)
  par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
  x <- barplot(rbind(dbrda_obj$aov_var_expl[1:3]), border = F,
    beside = T, las = 1, ylim = c(0, 25), ylab = "Variance explained (%)",
    xaxt = "n", yaxt = "n", space = c(rep(0.2, 3)))
  grid(nx = 0, ny = 5, col = gray(0.5))
  axis(2, lwd = 3)
  text(x = x, y = -1, names(dbrda_obj$aov_var_expl[1:3]), xpd = TRUE, srt = 45,
    adj = 1)
  box()
  dev.off()
}


### plot dbRDA variance per compartment ----------------------------------------
plot_dbRDA_var_per_compartment <- function(dbrda_pattern = "dbrda_fu_",
  outfile = "dbrda_var_compartment.pdf") {
  obj <- ls(pattern = paste0("^", dbrda_pattern), envir = globalenv())
  obj <- obj[c(1, 5, 3, 4, 2)]
  compartments <- c("Bulk soil", "Root zone", "Rhizosphere", "Root",
    "Endosphere")
  pdf(paste0("output/", outfile), w = 8, h = 8, pointsize = 36)
  par(mar = c(4, 4, 1, 1), las = 1, lwd = 3)
  x <- barplot(cbind(get(obj[1])$aov_var_expl[1:2],
    get(obj[2])$aov_var_expl[1:2], get(obj[3])$aov_var_expl[1:2],
    get(obj[4])$aov_var_expl[1:2], get(obj[5])$aov_var_expl[1:2]),
    beside = T, border = F,
    ylim = c(0, 25), xaxt = "n", yaxt = "n", ylab = "Variance explained (%)")
  grid(nx = 0, ny = 5, col = gray(0.5))
  axis(2, lwd = 3)
  text(x = apply(x, 2, mean), y = -1, compartments, xpd = TRUE, srt = 45,
    adj = 1)
  box()
  dev.off()
}


### table with proportions at selected taxonomic levels ------------------------
tax_prop_table <- function(m, t, s = NULL, n = 5) {
  tab <- t(apply(m, 1, function(x) tapply(x, t, sum)))
  if(!is.null(s)) tab <- apply(tab, 2, function(x) tapply(x, s, sum))
  tab <- tab[, names(sort(colSums(tab), decreasing = T))]
  tab <- tab / rowSums(tab)
  tab <- cbind(tab[, 1:n], others = rowSums(tab[, (n + 1):ncol(tab)]))
  return(tab)
}


### barplot with taxon proportions ---------------------------------------------
taxon_barplot <- function(cdm, sam, tax, taxlevel = "order", n = 5,
  outfile = "orders_sample.pdf") {
  options(warn = -1)
  tab <- tax_prop_table(cdm, tax[, taxlevel], sam$sample2, n = n)
  tab <- tab[c("Atba", "Atrs", "Atrh", "Atro", "Aten", "Chba", "Chrs", "Chrh",
    "Chro", "Chen", "Dvba", "Dvrs", "Dvrh", "Dvro", "Dven"), ]
  col_tax <- color(ncol(tab), start = 17)
  names(col_tax) <- colnames(tab)
  pdf(paste0("output/", outfile), w = 6, h = 2.5, pointsize = 12)
  par(mar = c(6, 4, 1, 6), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, -0.5))
  x <- barplot(t(tab), space = c(rep(c(0.3, rep(0.1, 4)), 3)), las = 1,
    col = col_tax, names.arg = rep(NA, nrow(tab)),
    ylab = "Proportion", border = NA)
  text(x, y = -0.1, rep(c("Bulk soil", "Root zone", "Rhizosphere", "Root",
    "Endosphere"), 3), xpd = TRUE, srt = 45, adj = 1)
  lines(c(x[1], x[5]), c(-0.75, -0.75))
  text((x[5] + x[1]) / 2, -0.85, "A. thaliana", font = 3, adj = 0.5)
  lines(c(x[6], x[10]), c(-0.75, -0.75))
  text((x[6] + x[10]) / 2, -0.85, "C. hirsuta", font = 3, adj = 0.5)
  lines(c(x[11], x[15]), c(-0.75, -0.75))
  text((x[11] + x[15]) / 2, -0.85, "D. verna", font = 3, adj = 0.5)
  legend("topright", legend = colnames(tab), fill = col_tax, bty = "n",
    inset = c(-0.28, 0), cex = 0.8, border = NA)
  dev.off()
  return(tab)
}


### plot reapds per taxa -------------------------------------------------------
reads_per_taxa <- function(cdm, tax, taxlevel = "order",
  outfile = "orders_overall.pdf") {
  cdm <- t(apply(cdm, 1, function(x) tapply(x, tax[, taxlevel], FUN = sum)))
  pdf(paste0("output/", outfile), h = 6, w = 15, pointsize = 12)
  par(mar = c(15, 4, 1, 1))
  barplot(sort(colSums(cdm), decreasing = T), las = 2, ylab = "Reads")
  dev.off()
}


### calculate number of taxa ---------------------------------------------------
number_of_taxa <- function(tax, taxlevel = "order") {
  x <- names(table(tax[, taxlevel]))
  x <- x[-(grep("^unclassified", x))]
  x <- x[-(grep("_Incertae_sedis", x))]
  message(paste0("Number of ", taxlevel, "-level taxa: ", length(x)))
}


### significance of taxa occurrence across compartments ------------------------
taxa_across_compartments <- function(cdm, sam, tax, taxlevel = "order",
  outfile = "orders_per_compartment.csv") {
  cdm <- t(apply(cdm, 1, function(x) tapply(x, tax[, taxlevel], FUN = sum)))
  sig <- as.data.frame(matrix(ncol = 4, nrow = ncol(cdm),
    dimnames = list(colnames(cdm), c("H", "df", "p", "p.adj"))))
  sig$H <- apply(cdm, 2,
    function(x) kruskal.test(x, sam$compartment)$statistic)
  sig$df <- apply(cdm, 2,
    function(x) kruskal.test(x, sam$compartment)$parameter)
  sig$p <- apply(cdm, 2,
    function(x) kruskal.test(x, sam$compartment)$p.value)
  sig$p.adj <- p.adjust(sig$p, method = "bon")
  write.table(sig, paste0("output/", outfile), sep = ";", col.names = NA)
  return(sig)
}


### calculate differential ASV abundance across compartments with LR tests -----
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
      try(m0 <- glm.nb(cdm_sub[, i] ~ data_sub$reads), silent = T)
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


### enrichment analysis --------------------------------------------------------
enrichment_analysis <- function(cdm, sam, tax, cdm_order, p_cutoff = 0.01,
  outfile = "dif_abund_fu") {
  require(gplots)
  dif <- dif_abund(cdm, sam)
  write.table(dif, paste0("output/", outfile, ".csv"), sep = ";",
    col.names = NA)
  dif_coeff <- data.frame(
    rhizosphere = as.numeric(dif[dif$comp == "rhizosphere", "coef"]),
    root = as.numeric(dif[dif$comp == "root", "coef"]),
    endosphere = as.numeric(dif[dif$comp == "endosphere", "coef"]))
  dif_coeff <- as.matrix(dif_coeff)
  rownames(dif_coeff) <- unique(dif$ASV)
  dif_pval <- data.frame(
    rhizosphere = as.numeric(dif[dif$comp == "rhizosphere", "Padj"]),
    root = as.numeric(dif[dif$comp == "root", "Padj"]),
    endosphere = as.numeric(dif[dif$comp == "endosphere", "Padj"]))
  dif_pval <- as.matrix(dif_pval)
  rownames(dif_pval) <- unique(dif$ASV)
  dif_top <- dif_coeff > 0 & dif_pval <= p_cutoff
  dif_low <- dif_coeff < 0 & dif_pval <= p_cutoff
  dif_top[is.na(dif_top)] <- F
  dif_low[is.na(dif_low)] <- F
  pdf(paste0("output/venn_", outfile, ".pdf"), w = 4, h = 2, pointsize = 12)
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
  xe <- table(tax[enriched_endo, "order"])
  xe <- xe * 100 / sum(xe)
  xr <- table(tax[enriched_root, "order"])
  xr <- xr * 100 / sum(xr)
  xrh <- table(tax[enriched_rhizo, "order"])
  xrh <- xrh * 100 / sum(xrh) 
  ye <- table(tax[depleted_endo, "order"])
  ye <- ye * 100 / sum(ye)
  yr <- table(tax[depleted_root, "order"])
  yr <- yr * 100 / sum(yr)
  yrh <- table(tax[depleted_rhizo, "order"])
  yrh <- yrh * 100 / sum(yrh)
  z <- cbind(Erh = xrh, Er = xr, Ee = xe, Drh = yrh, Dr = yr, De = ye)
  tab <- z[rownames(z) %in% colnames(cdm_order), ] 
  tab <- rbind(tab, others = colSums(z[!(rownames(z) %in%
    colnames(cdm_order)), ]))
  tab <- tab[colnames(cdm_order), ]
  col_tax <- color(nrow(tab), start = 17)
  names(col_tax) <- colnames(tab)
  pdf(paste0("output/", outfile, "_enriched_orders.pdf"), w = 2.4, h = 2.1,
    pointsize = 12)
  par(mar = c(4, 4, 1, 1), xpd = T, tck = -0.03, mgp = c(1.5, 0.5, 0))
  barplot(tab, las = 1, col = col_tax, ylab = "Proportion", border = NA,
    space = c(0.1, 0.1, 0.1, 0.3, 0.1, 0.1))
  dev.off()
  enriched_plant <- unique(c(enriched_endo, enriched_root))
  endo_tab <- cbind(order = as.character(tax[enriched_plant, "order"]),
    ASV = paste(gsub("_", " ", tax[enriched_plant, "species"]), enriched_plant))
  x <- cdm[, enriched_plant]
  x <- apply(x, 2, function(x) tapply(x, sam$species, sum))
  x[x > 0] <- 1
  x <- t(x)
  endo_tab <- cbind(endo_tab, x)
  endo_tab <- as.data.frame(endo_tab)
  endo_tab$in_rhizoplane <- as.numeric(enriched_plant %in% enriched_root)
  endo_tab$in_endosphere <- as.numeric(enriched_plant %in% enriched_endo)
  write.table(endo_tab, paste0("output/", outfile, "_endo.csv"), col.names = NA,
    sep = ";")
}


### end ------------------------------------------------------------------------
