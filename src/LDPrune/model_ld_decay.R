source(file.path("src", "R", "file_paths.R"))
library(tidyverse)
library(GGally)
library(pgda)
library(parallel)
install.packages("fda")

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

phys_order <- match(phys_data$snp$id, gen_data$snp$id)

snp_data <- tibble(
  phys_id = phys_data$snp$id, gen_id = gen_data$snp$id,
  chrom = phys_data$snp$chrom, phys_pos = (phys_data$snp$pos) / 1e6,
  gen_pos = (gen_data$snp$pos[phys_order]) / 100
)

library(mgcv)
plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
  xy <- data.frame(x = chrom_data$phys_pos, y = chrom_data$gen_pos)
  k <- 9
  unc <- gam(y ~ s(x, k = k, bs = "cr"), data = xy)
  sm <- smoothCon(s(x, k = k, bs = "cr"), xy, knots = NULL)[[1]]
  F <- mono.con(sm$xp)
  G <- list(
    y = xy$y, w = xy$y * 0 + 1, X = sm$X, C = matrix(0, 0, 0), S = sm$S,
    off = 0, sp = unc$sp, p = sm$xp, Ain = F$A, bin = F$b
  )
  p <- pcls(G)
  newx <- with(xy, data.frame(x = seq(min(x), max(x), length = 100)))
  fv <- Predict.matrix(sm, newx) %*% p
  newxy <- transform(newx, yhat = fv[,1])

  ggplot() +
    geom_point(aes(x = xy$x, y = xy$y), size = 0.5, alpha = 0.3) +
    geom_line(aes(x = newxy$x, y = newxy$yhat), col = "blue")
})

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Position in cM",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c(
    "Comparison of Order by Position between Location\n and Genetic Maps"
  )
)

# plot the matrix
png(
  file.path("chrom_pruned_marker_phys_vs_gen.png"),
  family = "Times New Roman", width = 120, height = 267, pointsize = 5,
  units = "mm", res = 300)
plots_matrix
dev.off()