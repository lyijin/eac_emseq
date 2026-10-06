#!/usr/bin/env Rscript

"> compile_atlas.R <

This script reads in tsv files with per-CpG meth and cov values, then summarises
those values across tissues.

The output is in the format that should (hopefully) be accepted by CelFiE-ISH as
an atlas.

Authored by Alice McAtamney, cosmetic fixes by Yi Jin Liew.
" -> doc

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(this.path)
})

setwd(this.path::here())

### Input/Output dirs
in_dir  <- "./filtered_beta_files/"
out_tsv <- "./Atlas.U250.per-CpG.hg38.tsv"

# Grouping metadata
meta_xlsx <- file.path(in_dir, "Loyfer_2023_TableS1.xlsx")

# List data files
files <- list.files(in_dir, pattern = "\\.atlas_cpg\\.tsv\\.gz$", full.names = TRUE)

# Clean grouping metadata names
meta <- read_excel(meta_xlsx, sheet = 1, col_names = FALSE)
cols <- meta[3,] %>%
  unlist() %>%
  as.character() %>%
  str_trim()
meta <- meta[-c(1:3), , drop = FALSE]
names(meta) <- cols

# create clean group mapping file
groups <- meta %>%
  transmute(
    sample_name_raw = str_trim(as.character(`Sample name`)),
    # normalise Excel sample names to match actual filenames
    sample_name = sample_name_raw %>%
      str_replace_all("Epithelium", "Epithelial"),
    group = str_trim(as.character(`Group`))
  ) %>%
  filter(!is.na(sample_name), sample_name != "") %>%
  distinct(sample_name, .keep_all = TRUE)


# Function: Infer sample name from .tsv filenames
# assumes that tissue comes after first "_" and before "-Z..."
infer_sample <- function(path) {
  b <- basename(path)
  sample_name <- str_match(b, "^[^_]+_(.+?)\\.atlas_cpg\\.tsv\\.gz$")[,2]
  sample_name <- str_trim(sample_name)
  if (is.na(sample_name) || sample_name == "") stop("Could not infer sample from filename: ", b)
  sample_name
}


# Function: map sample name to group (puts samples with no group into "UNMAPPED" group, which I drop later)
map_group <- function(sample_name, fallback = "UNMAPPED") {
  grp <- groups$group[match(sample_name, groups$sample_name)]
  ifelse(is.na(grp) | grp == "", fallback, grp)
}


# Function: read in .tsv, add inferred tissue column, and give celfie-specific colnames
read_one <- function(path) {
  sample_name <- infer_sample(path)
  grp <- map_group(sample_name)
  
  read_tsv(path, col_names = FALSE, show_col_types = FALSE) %>%
    transmute(
      tissue = grp,
      CHROM = X1,
      START = as.integer(X2),
      END   = as.integer(X3),
      METH  = as.numeric(X4),
      COV   = as.numeric(X5)
    )
}


# Use map_dfr to read in all .tsvs as one big tibble
# then calculate per cpg sum of meth and cov
dat_sum <- files %>% map_dfr(read_one) %>%
  group_by(tissue, CHROM, START, END) %>%
  summarise(
    METH = sum(METH, na.rm = TRUE),
    COV  = sum(COV, na.rm = TRUE),
    .groups = "drop" #ungroup
  )


# drop unmapped files
dat_sum <- dat_sum %>% filter(tissue != "UNMAPPED")

# Check if all groups are represented
not_present <- setdiff(groups$group %>% unique(), dat_sum$tissue %>% unique())
if (length(not_present) == 0) {
  message("All groups are present in dat_sum.")
} else {
  message("Groups not present: ",
          paste(sort(not_present), collapse = ", "))
}

# Pivot wider to get meth and cov column per tissue
atlas <- dat_sum %>% pivot_wider(
  id_cols = c(CHROM, START, END),
  names_from = tissue,
  values_from = c(METH, COV),
  names_glue = "{tissue}_{.value}",
  values_fill = list(METH = 0, COV = 0) # change cpgs with missing data to 0 meth 0 cov
)


# Order columns alphabetically so that tissue meth and cov cols are adjacent
coord_cols <- c("CHROM", "START", "END")
tissues <- sort(unique(dat_sum$tissue))
ordered_measure_cols <- as.vector(rbind(
  paste0(tissues, "_METH"),
  paste0(tissues, "_COV")
)) #keep _METH and _COV in the right order...

atlas <- atlas %>% select(all_of(coord_cols), all_of(ordered_measure_cols))

n_cpg <- atlas %>% distinct(CHROM, START, END) %>% nrow()
n_cols <- ncol(atlas)
n_groups <- (n_cols - 3)/2
message("Number of cpgs: ", n_cpg)
message("Number of groups: ", n_groups)


write_tsv(atlas, out_tsv)
message("Wrote: ", out_tsv)


# for replicability purposes
sessionInfo()
