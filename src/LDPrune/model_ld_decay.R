source(file.path("src", "R", "file_paths.R"))
library(tidyverse)
library(GGally)
library(pgda)

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

snp_data <- tibble(
  id = phys_data$snp$id, chrom = phys_data$snp$chrom,
  phys_pos = phys_data$snp$pos,
  gen_pos = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)]
)

plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
  pos_ratio <- chrom_data$phys_pos / chrom_data$gen_pos
  kept_index <- vector()
  for (i in 2:length(pos_ratio)) {
    if (
      pos_ratio[i] > (pos_ratio[i - 1] * 0.95) &
      pos_ratio[i] < (pos_ratio[i - 1] * 1.05)
    ) {
      kept_index <- c(kept_index, i)
    }
  }
  # kept_index <- 1:length(pos_ratio)
  chrom <- chrom_data$chrom[1]
  chrom_data[kept_index, ] %>%
    ggplot() +
      ylim(
          0,
          gen_data$max_lengths[[
            ifelse(grepl("1", chrom), "one", 
              ifelse(grepl("2", chrom), "two",
                ifelse(grepl("3", chrom), "three",
                  ifelse(grepl("4", chrom), "four",
                    ifelse(grepl("5", chrom), "five",
                      ifelse(grepl("6", chrom), "six", "seven")
                    )
                  )
                )
              )
            )
          ]] / 100
        ) +
      xlim(
        0,
        phys_data$max_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ] / 1e6
      ) +
      geom_point(aes((phys_pos / 1e6), (gen_pos / 100)), size = 0.5)
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