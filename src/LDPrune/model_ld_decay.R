source(file.path("src", "R", "file_paths.R"))
library(tidyverse)
library(GGally)
library(pgda)
library(parallel)
# install.packages("fda")

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

# snp_data <- tibble(
#   id = phys_data$snp$id, chrom = phys_data$snp$chrom,
#   phys_pos = phys_data$snp$pos,
#   gen_pos = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)]
# )
snp_data <- tibble(
  phys_id = phys_data$snp$id, gen_id = gen_data$snp$id,
  chrom = phys_data$snp$chrom, phys_pos = (phys_data$snp$pos) / 1e6,
  gen_pos = (gen_data$snp$pos) / 100
)
# snp_data <- tibble(
#   phys_id = phys_data$snp$id, gen_id = gen_data$snp$id,
#   chrom = phys_data$snp$chrom, phys_pos = phys_data$snp$pos / 1e6,
#   gen_pos = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)]
# )


marker_densities <- function(pos, n) {
  half_n <- n / 2
  half_n_less_1 <- half_n - 1
  gaps <- diff(pos)
  mclapply(1:length(pos), function (i) {
    nearest <- vector()
    if (i - half_n <= 0 && i + half_n_less_1 < length(gaps)) {
      nearest <- 1:(i + half_n_less_1)
    } else if (i - half_n > 0 && i + half_n_less_1 >= length(gaps)) {
      nearest <- (i - half_n):length(gaps)
    } else if (i - half_n < 0 && i + half_n_less_1 > length(gaps)) {
      nearest <- 1:length(gaps)
    } else {
      nearest <- (i - half_n):(i + half_n_less_1)
    }
    # calc the density of the regions aorund marker i and add it to the
    # densities vector
    gaps[nearest] %>% mean()
  }, mc.cores = detectCores()) %>% unlist()
}

# nrow(snp_data)
start.time <- Sys.time()
n <- 20
densities <- by(snp_data, snp_data$chrom, function (chrom_data) {
  half_n <- n / 20
  marker_densities <- tibble(
    marker = 1:nrow(chrom_data),
    pos = chrom_data$gen_pos,
    density = marker_densities(chrom_data$gen_pos, n)
  )
  while (max(marker_densities$density) > 1) {
    densist <- which.max(marker_densities$density)
    marker_densities <- marker_densities[-densist, ]
    if (densist - n <= 0 && densist + n < nrow(marker_densities)) {
      new_densities <- marker_densities(
        marker_densities$pos[1:(densist + n)], n
      )
      marker_densities$density[1:(densist + half_n)] <- new_densities[
        (length(new_densities) - (half_n + n)):(length(new_densities) - half_n)
      ]
    } else if (densist - n > 0 && densist + n >= nrow(marker_densities)) {
      new_densities <- marker_densities(
        marker_densities$pos[(densist - n):nrow(marker_densities)], n
      )
      marker_densities$density[
        (densist - half_n):nrow(marker_densities)
      ] <- new_densities[(1 + half_n):(1 + n + half_n)]
    } else if (densist - n <= 0 && densist + n >= nrow(marker_densities)) {
       marker_densities$density <- marker_densities(marker_densities$pos, n)
    } else {
      new_densities <- marker_densities(
        marker_densities$pos[(densist - n):(densist + n)], n
      )
      marker_densities$density[
        (densist - half_n):(densist + half_n)
      ] <- new_densities[(1 + half_n):(1 + n + half_n)]
    }
  }
  marker_densities
})
end.time <- Sys.time()
end.time - start.time

densities[[1]]

# which.max(densities)

# all.equal(densities1, densities2)

# library(mgcv)
# plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
#   xy <- data.frame(x = chrom_data$phys_pos, y = chrom_data$gen_pos)
#   k <- 9
#   unc <- gam(y ~ s(x, k = k, bs = "cr"), data = xy)
#   sm <- smoothCon(s(x, k = k, bs = "cr"), xy, knots = NULL)[[1]]
#   F <- mono.con(sm$xp)
#   G <- list(
#     y = xy$y, w = xy$y * 0 + 1, X = sm$X, C = matrix(0, 0, 0), S = sm$S,
#     off = 0, sp = unc$sp, p = sm$xp, Ain = F$A, bin = F$b
#   )
#   p <- pcls(G)
#   newx <- with(xy, data.frame(x = seq(min(x), max(x), length = 100)))
#   fv <- Predict.matrix(sm, newx) %*% p
#   newxy <- transform(newx, yhat = fv[,1])
#   # print(newx)
#   # tibble(x = newx, y = yhat)
#   ggplot() +
#     geom_point(aes(x = xy$x, y = xy$y), size = 0.5, alpha = 0.3) +
#     geom_line(aes(x = newxy$x, y = newxy$yhat), col = "blue")
#   # plot(x ~ y, data = xy, pch = 16)
#   # lines(x ~ yhat, data = newxy, col = "blue")
#   # chrom_data %>% ggplot() +
#   #   geom_point(aes(phys_pos, gen_pos), size = 0.5) +
#   #   geom_smooth(aes(phys_pos, gen_pos), size = 0.5, method = lm, formula = )
#   # break
# })

# plots_matrix <- ggmatrix(
#   plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Position in cM",
#   xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
#   title = str_c(
#     "Comparison of Order by Position between Location\n and Genetic Maps"
#   )
# )

# # plot the matrix
# png(
#   file.path("chrom_pruned_marker_phys_vs_gen.png"),
#   family = "Times New Roman", width = 120, height = 267, pointsize = 5,
#   units = "mm", res = 300)
# plots_matrix
# dev.off()