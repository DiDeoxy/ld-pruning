# import file paths and functions
source(file.path("src", "R", "file_paths.R"))
import::from(dplyr, "as_tibble")
import::from(igraph, "graph_from_edgelist", "max_cliques")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse", "snpgds_sample_subset", "snpgds_snp_subset")
import::from(
  SNPRelate, "snpgdsClose", "snpgdsIBS", "snpgdsLDpruning", "snpgdsOpen",
  "snpgdsSelectSNP"
)
import::from(stringr, "str_c")

# eliminate those individuals that show identity by state
# (IBS, fractional identity) greater than 0.99
wheat_gds <- snpgdsOpen(file.path(gds, "full_phys.gds"))
IBS <- snpgdsIBS(wheat_gds, autosome.only = F)
snpgdsClose(wheat_gds)

pairs <- which(IBS$ibs >= 0.99, arr.ind = T)
pairs <- cbind(IBS$sample.id[pairs[, 1]], IBS$sample.id[pairs[, 2]])

indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i, 1] == pairs[i, 2]) {
    indices <- c(indices, i)
  }
}
pairs <- pairs[-indices, ]

# print out a table of the maximal cliques constructed from the pairs
graph_from_edgelist(pairs) %>%
  max_cliques() %>%
  lapply(names) %>%
  lapply(`length<-`, max(lengths(.))) %>%
  do.call(rbind, .) %>%
  cbind(str_c("Clique ", 1:nrow(.)), .) %>%
  as_tibble() %>%
  write.table(
    file.path("results", "clique_table.csv"), sep = ",",
    row.names = FALSE, quote = FALSE,
    col.names = c("Clique", str_c("Cultivar ", 1:(ncol(.) - 1)))
  )

# used clique_table.csv to identify cultivars to prune from each clique
NILs <- c(
  "PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1", "BW395",
  "PT616", "BW427 1", "BW492", "BW948", "Carberry 1", "CDC Stanley 1",
  "PT754", "SWS349", "Somerset 1", "Stettler 1", "SWS241", "SWS345",
  "AC Reed 1", "SWS87", "SWS390", "SWS408", "SWS410"
)

# mr pruned phys map
wheat_data <- snpgds_parse(file.path(gds, "full_phys.gds"))
sample_index <- match(NILs, wheat_data$sample$id)
# create gds object without the NILs
snpgds_sample_subset(
  wheat_data, file.path(gds, "full_phys_sample_subset.gds"), sample_index
)

# mr pruned gen map
wheat_data <- snpgds_parse(file.path(gds, "full_gen.gds"))
sample_index <- match(NILs, wheat_data$sample$id)
snpgds_sample_subset(
  wheat_data, file.path(gds, "full_gen_sample_subset.gds"), sample_index
)

# identify snps with a maf below 0.05
wheat_gds <- snpgdsOpen(file.path(gds, "full_phys_sample_subset.gds"))
kept_id <- snpgdsSelectSNP(
  wheat_gds, maf = 0.05, missing.rate = 0.10, autosome.only = F
)
snpgdsClose(wheat_gds)

# reomve these markers from the phys map
wheat_data <- snpgds_parse(file.path(gds, "full_phys_sample_subset.gds"))
kept_index <- match(kept_id, wheat_data$snp$id)
snpgds_snp_subset(wheat_data, file.path(phys_gds), kept_index)

# remove these markers from the gen map
wheat_data <- snpgds_parse(file.path(gds, "full_gen_sample_subset.gds"))
kept_index <- match(kept_id, wheat_data$snp$id) %>% sort()
snpgds_snp_subset(wheat_data, file.path(gen_gds), kept_index)