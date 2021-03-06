# raw data
raw_data <- file.path("data", "raw")
cultivars <- file.path(raw_data, "cultivars")
markers <- file.path(raw_data, "markers")

################################################################################
# intermediate data
intermediate <- file.path("data", "intermediate")
ifelse(! dir.exists(intermediate), dir.create(intermediate), FALSE)

inter_markers <- file.path(intermediate, "markers")
ifelse(! dir.exists(inter_markers), dir.create(inter_markers), FALSE)

gds <- file.path(intermediate, "gds")
ifelse(! dir.exists(gds), dir.create(gds), FALSE)

phys_gds <- file.path(gds, "maf_mr_sample_filtered_phys.gds")
gen_gds <- file.path(gds, "maf_mr_sample_filtered_gen.gds")