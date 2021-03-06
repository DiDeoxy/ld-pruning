# import file paths and functions
source(file.path("src", "R", "file_paths.R"))
import::from(ape, "read.gff")
import::from(
  dplyr, "arrange", "do", "filter", "group_by", "left_join", "n", "rename",
  "select", "ungroup"
)
import::from(magrittr, "%>%")
import::from(
  readr, "col_character", "col_double", "col_factor", "col_integer", "read_csv",
  "read_rds", "type_convert", "write_rds"
)
import::from(stringr, "str_c")
import::from(tibble, "as_tibble", "tibble")

# load the filtered gen map
poz_filtered <- read_rds(file.path(inter_markers, "pozniak_filtered_map.rds"))

# create a vector of the chromosomes
chrs <- str_c(
  "chr",
  outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
    t() %>% as.vector()
)

# parse the GFF3 format file of the wheat 90K snp chip physical map positions
alignments <- data.frame()
for (chr in chrs) {
  chr_feats <- read.gff(
    file.path(
      markers, "90K_RefSeqv1_probe_alignments", str_c("Infinium90K-", chr, ".gff3")
    )
  )
  parsed_pos <- apply(chr_feats, 1, function (feat) {
    attributes <- strsplit(feat[9], split = ";")[[1]]
    name <- strsplit(attributes[2], split = "=")[[1]][2]
    ID <- strsplit(attributes[4], split = "=")[[1]][2]
    coverage <- strsplit(attributes[1], split = "=")[[1]][2]
    per_id <- strsplit(attributes[3], split = "=")[[1]][2]
    chrom <- substr(feat[1], 4, nchar(feat[1]))
    pos <- floor( (as.integer(feat[4]) + as.integer(feat[5]) ) / 2)
    return(c(name, ID, chrom, pos, coverage, per_id, use.names = F))
  }) %>% t()
  alignments <- rbind(alignments, parsed_pos, stringsAsFactors = F)
}

## format the alignemts into a tibble with named columns
alignments <- alignments %>%
  as_tibble() %>%
  type_convert(
    col_type = list(
      col_character(), col_character(), col_character(), col_integer(),
      col_double(), col_double()
    )
  ) %>%
  rename(
    marker = V1, ID = V2, chrom = V3, pos = V4, coverage = V5, per_id = V6
  )

# a function for returning the best alignment on the correct chromosome
best_alignment <- function(aligns, poz_filtered) {
  # finds if the same marker is on the filtered gen map
  poz_marker <- filter(poz_filtered, marker == aligns$marker[1])
  if (nrow(poz_marker)) {
    # finds if the marker has a single phys alignment on a chromosome 
    # equivalent to the gen map linkage group
    aligns_chrom <- filter(aligns, chrom == poz_marker$group)
    if (nrow(aligns_chrom) == 1) {
        return(aligns_chrom)
    } else {
        return(tibble())
    }
  } else {
    return(tibble())
  }
}

# create the physical map by filtering based on quality and using the
# best_alignment function
phys_map <- alignments %>%
  filter(coverage >= 90, per_id >= 98) %>%
  group_by(marker) %>%
  do(best_alignment(., poz_filtered)) %>%
  select(marker, chrom, pos) %>%
  arrange(chrom, pos)
maps <- phys_map %>%
  left_join(poz_filtered, by = "marker") %>%
  rename(phys_pos = pos.x, gen_pos = pos.y) %>%
  select(marker, chrom, phys_pos, gen_pos)

# format the genotype data into the proper format for snpgds format
genotypes <- read_csv(
  file.path(markers, "Jan_6_wheat_genotypes_curtis.csv")
) %>%
  select(-X1, -X3, -X4, -X5, -Name) %>%
  .[-1:-2, ] %>%
  rename(marker = X2) %>%
  replace(. == "C1", 0) %>%
  replace(. == "c1", 0) %>%
  replace(. == "C2", 2) %>%
  replace(. == "NC", 3)

# combine the phys map with the genotypes 
maps_genotypes <- maps %>%
  left_join(genotypes) %>%
  type_convert() %>%
  group_by(chrom, phys_pos)
nrow(maps_genotypes) %>%
  str_c("Num markers: ", .) %>%
  print()

# find the groups of markers which have the same position in a chromosome
duplicates <- maps_genotypes %>%
  group_by(chrom, phys_pos) %>%
  filter(n() >= 2)
nrow(duplicates) %>%
  str_c("Num markers with identical postions: ", .) %>%
  print()

# remove at least one marker mapping to the same position in phys map
unique_marker <- function(markers) {
  if (nrow(markers) == 1) {
    return(markers[1, ])
  } else if (nrow(markers) == 2) {
    marker1 <- markers[1, 5:ncol(markers)]
    marker2 <- markers[2, 5:ncol(markers)]
    per_diff <- (sum(abs(marker1 - marker2), na.rm = TRUE) / 2) / 100
    if (per_diff > 0.01) {
        return(tibble())
    } else {
        return(markers[1, ])
    }
  } else {
    return(tibble())
  }
}

# filter the duplicates to identify the number retained
duplicates_filtered <- duplicates %>%
  do(unique_marker(.))
nrow(duplicates_filtered) %>%
  str_c(
    "Of markers with identical postions num retained after pruning: ", .
  ) %>% print()

# de-duplicate the full data set
maps_genotypes_deduplicated <- maps_genotypes %>%
  do(unique_marker(.)) %>%
  ungroup()
nrow(maps_genotypes_deduplicated) %>%
  str_c("Num markers remaining after pruning: ", .) %>% print()

# write out the maps with genotypes
write_rds(
  maps_genotypes_deduplicated, 
  file.path(inter_markers, "maps_genotypes.rds")
)
