marker_densities <- function(pos, n) {
  half_n <- floor(n / 2)
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
    # calc the density of the regions around marker i and add it to the
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