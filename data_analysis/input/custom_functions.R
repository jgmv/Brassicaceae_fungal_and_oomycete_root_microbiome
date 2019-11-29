### In-house functions

### Plot OTU abundances as needles #############################################
plotNeedle <- function(data, col = gray(0.5), ylab = "Reads", xlab = "OTU") {
  plot(data, type = "n", ylim = c(0, max(data) * 1.2), bty = "n", xlab = xlab, 
    ylab = ylab, las = 1)
  for (i in 1:length(data)) {
    x <- data[i]
    lines(rep(i, 2), cbind(0, x), col = col)
    if (x > max(data) * 0.1) text(i, x + par("usr")[4] * 0.1, names(data)[i])
  }
}
################################################################################

### Automatically perform regressions, corrected by read numbers ###############
plotModel <- function(var, exp, cons) {
  require(car)
  model <- lm(var ~ cons + exp)
  model_pred <- predict(model, type = "response")
    par(mar = c(4, 4, 2, 2), lwd = 2)
    crPlots(model, terms = ~.-cons, grid = F, smoother = F, ylab = NA,
      xlab = NA, axes = F, col.lines = 1, pch = 16,
        col = alpha(gray(0.25), 0.5), lty = 1, cex = 1.5)
    axis(1, lwd = 2)
    axis(2, lwd = 2)
    box()
    mtext(paste("adj. Rsq = ", round(summary(model)$adj.r.squared, 2),
            ", P = ", round(anova(model)[[5]][2], 3), sep = ""),
            adj = 1, cex = 1)
    return(model)
}
################################################################################

### Alternative to 'plotModel' #################################################
plotModel2 <- function(var, exp, cons) {
  require(car)
  model <- lm(var ~ cons + exp)
  model_pred <- predict(model, type = "response")
    par(mar = c(2.5, 2, 1.25, 1.75), lwd = 2, cex = 2)
    crPlots(model, terms = ~.-cons, grid = F, smoother = F, ylab = NA,
      xlab = NA, axes = F, col.lines = 1, pch = 16,
      col = alpha(gray(0.25), 0.5), lty = 1)
    box()
}
################################################################################

### Remove unclassified and undetermined taxa from taxonomic lists #############
replaceUnclassified <- function(tax, rep) {
  levels(tax) <- c(levels(tax), rep)
  tax[grep("n.d.", tax)] <- rep
  tax[grep("Incertae", tax)] <- rep
  tax[grep("unclassified", tax)] <- rep
  tax <- droplevels(tax)
  tax <- sort(table(tax), decreasing = T)
  return(tax)
}
################################################################################

### Plot BLAST identity percentages as polygons ################################
blastPolygon <- function(pident, col) {
  pident[is.na(pident)] <- 0
  pident <- sort(pident, decreasing = T)
  polygon(c(0, pident, 0), col = col, border = NA)
}
################################################################################

### Forward select variables as in Borcard et al. 2011 #########################
fwdSel <- function(cdm, var) {
  require(packfor)
  m <- rda(cdm ~ as.matrix(var), na.action = na.omit)
  m_anova <- anova.cca(m)
  print(m_anova)
  if (m_anova[4]$Pr[1] > 0.05) {
    message("Non-significant RDA model, P > 0.05")
    break    
  } else {
    var <- na.omit(var)
    r2 <- RsquareAdj(m)$adj.r.squared
    fwd_var <- forward.sel(cdm[rownames(var), ], var, adjR2thresh = r2)
    print(fwd_var)
    fwd_var_m <- as.matrix(var[, fwd_var$variables])
    rownames(fwd_var_m) <- rownames(var)
    return(fwd_var_m)
  }
}
################################################################################

### Test significance of variation partition components (3 components) #########
testVP3 <- function(vp, cdm = NULL) {
  # retrieve tables from vp
  y  <- eval(parse(text = vp$call[2]))
  X1 <- as.matrix(eval(parse(text = vp$call[3])))
  X2 <- as.matrix(eval(parse(text = vp$call[4])))
  X3 <- as.matrix(eval(parse(text = vp$call[5])))

  # create an output table
  tab <- rbind(vp$part[[1]][1:4], vp$part[[2]][1:4], vp$part[[3]][1:4])
  tab$percVar <- tab[, "Adj.R.square"] * 100
  tab$P <- rep(NA, nrow(tab))
  #showvarparts(3)

  if(class(y) == "dist") {
    tab[7, "P"]  <- anova.cca(capscale(y ~ X1 + X2 + X3))$Pr[1]
    tab[1, "P"]  <- anova.cca(capscale(y ~ X1))$Pr[1]
    tab[2, "P"]  <- anova.cca(capscale(y ~ X2))$Pr[1]
    tab[3, "P"]  <- anova.cca(capscale(y ~ X3))$Pr[1]
    tab[4, "P"]  <- anova.cca(capscale(y ~ X1 + X2))$Pr[1]
    tab[5, "P"]  <- anova.cca(capscale(y ~ X1 + X3))$Pr[1]
    tab[6, "P"]  <- anova.cca(capscale(y ~ X2 + X3))$Pr[1]
    tab[8, "P"]  <- anova.cca(capscale(y ~ X1 + Condition(X2) +
      Condition(X3)))$Pr[1]
    tab[9, "P"]  <- anova.cca(capscale(y ~ X2 + Condition(X1) +
      Condition(X3)))$Pr[1]
    tab[10, "P"] <- anova.cca(capscale(y ~ X3 + Condition(X1) +
      Condition(X2)))$Pr[1]
    tab[16, "P"] <- anova.cca(capscale(y ~ X1 + Condition(X3)))$Pr[1]
    tab[17, "P"] <- anova.cca(capscale(y ~ X1 + Condition(X2)))$Pr[1]
    tab[18, "P"] <- anova.cca(capscale(y ~ X2 + Condition(X3)))$Pr[1]
    tab[19, "P"] <- anova.cca(capscale(y ~ X2 + Condition(X1)))$Pr[1]
    tab[20, "P"] <- anova.cca(capscale(y ~ X3 + Condition(X1)))$Pr[1]
    tab[21, "P"] <- anova.cca(capscale(y ~ X3 + Condition(X2)))$Pr[1]
  } else {
    tab[7, "P"]  <- anova.cca(rda(y ~ X1 + X2 + X3))$Pr[1]
    tab[1, "P"]  <- anova.cca(rda(y ~ X1))$Pr[1]
    tab[2, "P"]  <- anova.cca(rda(y ~ X2))$Pr[1]
    tab[3, "P"]  <- anova.cca(rda(y ~ X3))$Pr[1]
    tab[4, "P"]  <- anova.cca(rda(y ~ X1 + X2))$Pr[1]
    tab[5, "P"]  <- anova.cca(rda(y ~ X1 + X3))$Pr[1]
    tab[6, "P"]  <- anova.cca(rda(y ~ X2 + X3))$Pr[1]
    tab[8, "P"]  <- anova.cca(rda(y ~ X1 + Condition(X2) + Condition(X3)))$Pr[1]
    tab[9, "P"]  <- anova.cca(rda(y ~ X2 + Condition(X1) + Condition(X3)))$Pr[1]
    tab[10, "P"] <- anova.cca(rda(y ~ X3 + Condition(X1) + Condition(X2)))$Pr[1]
    tab[16, "P"] <- anova.cca(rda(y ~ X1 + Condition(X3)))$Pr[1]
    tab[17, "P"] <- anova.cca(rda(y ~ X1 + Condition(X2)))$Pr[1]
    tab[18, "P"] <- anova.cca(rda(y ~ X2 + Condition(X3)))$Pr[1]
    tab[19, "P"] <- anova.cca(rda(y ~ X2 + Condition(X1)))$Pr[1]
    tab[20, "P"] <- anova.cca(rda(y ~ X3 + Condition(X1)))$Pr[1]
    tab[21, "P"] <- anova.cca(rda(y ~ X3 + Condition(X2)))$Pr[1]
  }
  return(tab)
}
################################################################################

### Plot Euler diagrams from variation partition object ########################
plotVpEuler3 <- function(vp, names = c("X1", "X2", "X3"),
  col = c("brown3", "skyblue3", "orange")) {
  require(venneuler)
  sec <- vp$part$indfract[-8, 3]
  sec <- ifelse(sec < 0, 0, sec)
  res <- vp$part$indfract[8, 3]
  names(sec)<- c(names, paste(names[1], names[2], sep = "&"),
    paste(names[2], names[3], sep = "&"), paste(names[1], names[3], sep = "&"),
    paste(names[1], names[2], names[3], sep = "&"))
  vd <- venneuler(sec, residuals = res)
  plot(vd, col = col)
  mtext(paste("Residuals = ", round(res, 2), sep = ""), 1, cex = 1)
} 
################################################################################

### Plot different components of multispecies-GLM model ########################
plotEffects <- function(x) {
  a <- max(abs(x))
  colort <- colorRampPalette(c("darkred", "white", "darkblue")) 
  plot_spp1 <- levelplot(t(as.matrix(x[,1:3])), xlab = "Climate", ylab = NULL,
    colorkey = F, col.regions = colort(100), at = seq(-a, a, length = 100),
    scales = list(y = list(draw = T), x = list(rot = 0), tck = c(1, 0)))
  plot_spp2 <- levelplot(t(as.matrix(x[, 4:6])), xlab = "Soil", ylab = NULL,
    colorkey = F, col.regions = colort(100), at = seq(-a, a, length=100),
    scales = list(y = list(draw = F), x = list(rot = 0), tck = c(1,0)))
  plot_spp3 <- levelplot(t(as.matrix(x[, 7:11])), xlab="Spatial", ylab = NULL,
    col.regions = colort(100), at = seq(-a, a, length = 100),
    scales = list(y = list(draw = F), x = list(rot = 0), tck = c(1, 0)))
  print(plot_spp1, split = c(1, 1, 3, 1), more = T)
  print(plot_spp2, split = c(2, 1, 3, 1), more = T)
  print(plot_spp3, split = c(3, 1, 3, 1))
}
################################################################################

### Performs PERMANOVA for a list of provided variables ########################
many_adonis <- function(data, var) {
  require(vegan)
  tab <- matrix(NA, ncol = 4, nrow = ncol(var),
    dimnames = list(colnames(var), c("Df", "F", "R2", "P")))
  tab <- as.data.frame(tab)
  for (i in 1:ncol(var)) {
    fac <- na.omit(var[, i])
    if (!is.null(attr(fac, "na.action"))) {
      mat <- data[-attr(fac, "na.action"), ]
    } else {
      mat <- data
    }
    x <- adonis(mat ~ fac)
    tab[i ,"Df"] <- paste(x$aov.tab$Df[1],x$aov.tab$Df[3], sep = ",")
    tab[i ,"F"]  <- x$aov.tab$F[1]
    tab[i ,"R2"] <- x$aov.tab$R2[1]
    tab[i ,"P"]  <- x$aov.tab$P[1]
  }
  return(tab)
}
################################################################################

### Prepare taxonomic data frame for input into metacoder ######################
tax_to_string <- function(tax, root = "Root"){
  tax <- cbind(root = rep(root, nrow(tax)), tax)
  result <-  apply(tax, 1, paste, collapse = ";")
  return(result)
}
################################################################################

### Calculate number of reads per lineage in metacoder's taxmap object #########
reads_per_lineage <- function(obj, var = "obs_reads") {
  vapply(obs(obj, recursive = T, simplify = F), function(x) 
    sum(obs_data(obj, row_subset = x, col_subset = var, drop = F)), numeric(1))
}
################################################################################

### Calculate number of groups per lineage in metacoder's taxmap object ########
get_groups <- function(obj, x, var) {
  g <- obs_data(obj, row_subset = x, col_subset = var, drop = F)
  g[g > 0] <- 1
  return(sum(g))
}
groups_per_lineage <- function(obj, var = "obs_reads") {
  vapply(obs(obj, recursive = T, simplify = F), function(x) 
    get_groups(obj, x, var), numeric(1))
}
################################################################################

### Prepare a metacoder's taxmap with diversity and abundance per group ########
cdm_to_taxmap <- function(tax, cmd, factor) {
  taxmap     <- extract_taxonomy(tax, class_sep = ";")
  taxmap <- mutate_obs(taxmap, obs_reads = colSums(cmd))
  taxmap <- mutate_taxa(taxmap,
    prop_reads = sqrt(reads_per_lineage(taxmap, "obs_reads")))
  
  for(i in levels(factor)) {
    taxmap <- mutate_obs(taxmap,
      reads = colSums(cmd[factor == i, ]))
    names(taxmap$obs_data) <-
     c(names(taxmap$obs_data)[1:(length(names(taxmap$obs_data)) - 1)], i) 

    taxmap <- mutate_taxa(taxmap,
      reads = sqrt(reads_per_lineage(taxmap, i)))
    names(taxmap$taxon_data) <-
     c(names(taxmap$taxon_data)[1:(length(names(taxmap$taxon_data)) - 1)],
     paste(i, "abund", sep = "_"))

    taxmap <- mutate_taxa(taxmap,
      group = groups_per_lineage(taxmap, i))
    names(taxmap$taxon_data) <-
     c(names(taxmap$taxon_data)[1:(length(names(taxmap$taxon_data)) - 1)],
     paste(i, "groups", sep = "_"))
  }
  return(taxmap)
}
################################################################################

### Function to plot PCoAs with factors ########################################
plot_pcoa <- function(pcoa, expl, index = sample) {
  fac <- droplevels(index[rownames(pcoa$points), ])
  date_fill <- col_comp[as.character(fac$compartment)]
  date_fill[fac$date == "February"] <- 0
  plot(scores(pcoa),
    type = "n",
    xlab = paste("PCo1 (", expl[1], "%)", sep = ""),
    ylab = paste("PCo2 (", expl[2], "%)", sep = ""))
  lines(c(0, 0), c(-1, 1), col = gray(0.5), lty = 3)
  lines(c(-1, 1), c(0, 0), col = gray(0.5), lty = 3)
  points(scores(pcoa),
    pch = as.numeric(fac$species) + 20,
    col = col_comp[as.character(fac$compartment)],
    bg = date_fill)
}
################################################################################

### Function to plot CAPs with factors #########################################
plot_cap <- function(cap, expl, index = sample) {
  fac <- droplevels(index[rownames(cap$CCA$u), ])
  date_fill <- col_comp[as.character(fac$compartment)]
  date_fill[fac$date == "February"] <- 0
  plot(scores(cap)$sites,
    type = "n",
    xlab = paste("Axis 1 (", expl[1], "%)", sep = ""),
    ylab = paste("Axis 2 (", expl[2], "%)", sep = ""))
  lines(c(0, 0), c(-10, 10), col = gray(0.5), lty = 3)
  lines(c(-10, 10), c(0, 0), col = gray(0.5), lty = 3)
  points(scores(cap)$sites,
    pch = as.numeric(fac$species) + 20,
    col = col_comp[as.character(fac$compartment)],
    bg = date_fill)
}
################################################################################

### Modify names in taxonomy file ##############################################
# version: 2019-02-20
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
################################################################################

### Boxplot with data points ###################################################
# version: 2019-02-20
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
################################################################################

### Result of multipatt to network #############################################
indssp_to_network <- function(x, p_val = 0.01) {
  require(igraph)
  tab <- na.omit(x$sign[x$sign$p.value < p_val, ])
  conn <- c()
  for(i in colnames(tab)[1:(which(colnames(tab) == "index") - 1)]) {
    for(j in rownames(tab[tab[, i] == 1, ])) {
      conn <- c(conn, j, i) 
    }
  }
  net <- graph(conn, directed = F)
  net$object <- rep("species", length(V(net)))
  net$object[grep("s.", V(net)$name)] <- "sample"
  net$shape <- rep("circle", length(V(net)))
  net$shape[net$object == "sample"] <- "square"
  return(net)
}
################################################################################

### Return table with proportions at selected taxonomic level ##################
tax_prop_table <- function(m, t, s = NULL, n = 5) {
  tab <- t(apply(m, 1, function(x) tapply(x, t, sum)))
  if(!is.null(s)) tab <- apply(tab, 2, function(x) tapply(x, s, sum))
  tab <- tab[, names(sort(colSums(tab), decreasing = T))]
  tab <- tab / rowSums(tab)
  tab <- cbind(tab[, 1:n], others = rowSums(tab[, (n + 1):ncol(tab)]))
  return(tab)
}
################################################################################

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


### End ########################################################################
