#!/usr/bin/env Rscript

"> 03_check_overlap_cell_type_classification.R <

Are the differentially hypermethylated promoters that appear in the
inter-cell-type comparisons similar to those from inter-classification
comparisons?

Focus on (all) six cell-type comparisons:
- Columnar vs. Squamous
- Lymphoid vs. Squamous
- Erythroid vs. Squamous
- Lymphoid vs. Columnar
- Erythroid vs. Columnar
- Erythroid vs. Lymphoid

and the six 'relevant' classification comparisons:
- NDBE vs. NSq
- EAC vs. NSq
- LGD vs. NDBE
- HGD/EAC vs. NDBE
- HGD vs. LGD
- EAC vs. HGD
" -> doc

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(this.path)
})

setwd(this.path::here())
CELLTYPE_DMR_FOLDER <- '../15_dmr_by_cell_types/filtered_dmrs/'
CELLTYPE_FILES <- c(
  paste0(CELLTYPE_DMR_FOLDER, 'Columnar_vs_Squamous.hypermeth.tsv.gz'),
  paste0(CELLTYPE_DMR_FOLDER, 'Lymphoid_vs_Squamous.hypermeth.tsv.gz'),
  paste0(CELLTYPE_DMR_FOLDER, 'Erythroid_vs_Squamous.hypermeth.tsv.gz'),
  paste0(CELLTYPE_DMR_FOLDER, 'Lymphoid_vs_Columnar.hypermeth.tsv.gz'),
  paste0(CELLTYPE_DMR_FOLDER, 'Erythroid_vs_Columnar.hypermeth.tsv.gz'),
  paste0(CELLTYPE_DMR_FOLDER, 'Erythroid_vs_Lymphoid.hypermeth.tsv.gz')
)
CLASSIFICATION_DMR_FOLDER <- '../11_analyse_dmr_metilene/filtered_dmrs/'
CLASSIFICATION_FILES <- c(
  paste0(CLASSIFICATION_DMR_FOLDER, 'allB_vs_allN.hypermeth.tsv.gz'),
  paste0(CLASSIFICATION_DMR_FOLDER, 'allC_vs_allN.hypermeth.tsv.gz'),
  paste0(CLASSIFICATION_DMR_FOLDER, 'allL_vs_allB.hypermeth.tsv.gz'),
  paste0(CLASSIFICATION_DMR_FOLDER, 'allHC_vs_allB.hypermeth.tsv.gz'),
  paste0(CLASSIFICATION_DMR_FOLDER, 'allH_vs_allL.hypermeth.tsv.gz'),
  paste0(CLASSIFICATION_DMR_FOLDER, 'allC_vs_allH.hypermeth.tsv.gz')
)

# error out if any of the files do not exist
stopifnot(file.exists(CELLTYPE_FILES))
stopifnot(file.exists(CLASSIFICATION_FILES))


# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# store hypermeth promoters in a list
hmp_list <- list()

# read all 12 relevant files into the list of vectors
cellfs <- CELLTYPE_FILES
classfs <- CLASSIFICATION_FILES
for (cf in c(cellfs, classfs)) {
  # get list of unique promoters and store it in a char array
  # note: array is unsorted to preserve the order in the original files,
  #       i.e., ordered by delta difference
  hmp_list[[cf]] <- read_tsv(cf, show_col_types=FALSE) |>
    filter(promoter == 1) |>
    pull(gene_name_prot) |>
    unique()
}

# create a 6x6 data frame (not tibble, i want rownames) to store overlap
overlap_df <- data.frame(matrix(NA_character_, nrow=6, ncol=6))
rownames(overlap_df) <- cellfs
colnames(overlap_df) <- classfs

# replace the NAs column-by-column
for (classf in classfs) {
  overlap_df[[classf]] <- sapply(cellfs, function(x) intersect(hmp_list[[classf]], hmp_list[[x]]))
}

# prettify row/col names, now that cells are populated
rownames(overlap_df) <- cellfs |>
  str_replace('.*/', '') |>
  str_replace('\\..*$', paste0(' (n = ', sapply(hmp_list[cellfs], length), ')')) |>
  str_replace('_vs_', ' vs. ')
colnames(overlap_df) <- classfs |>
  str_replace('.*/', '') |>
  str_replace('\\..*$', paste0(' (n = ', sapply(hmp_list[classfs], length), ')')) |>
  str_replace('_vs_', ' vs. ') |>
  str_replace('allHC', 'HGD/EAC') |>
  str_replace('allN', 'NSq') |>
  str_replace('allB', 'NDBE') |>
  str_replace('allL', 'LGD') |>
  str_replace('allH', 'HGD') |>
  str_replace('allC', 'EAC')

# view counts shared between the rows/cols
counts_df <- overlap_df
for (classf in colnames(counts_df)) {
  counts_df[[classf]] <- sapply(counts_df[[classf]], length)
}
counts_df
write.table(counts_df, paste0('raw_supptable27a.tsv'), quote=FALSE, sep='\t',
            col.names=NA)

# view top 10 gene names of hypermeth promoters common to both comparisons
top10_df <- overlap_df
for (classf in colnames(top10_df)) {
  for (cellf in cellfs) {
    top10_df[[classf]][[cellf]] <-
      c(top10_df[[classf]][[cellf]][1:10],
        paste('+', length(top10_df[[classf]][[cellf]])-10, 'others')) |>
      str_flatten(collapse=';')
  }
}
top10_df[] <- lapply(top10_df, as.character)
write.table(top10_df, paste0('raw_supptable27c.tsv'), quote=FALSE, sep='\t',
            col.names=NA)


# for replicability purposes
sessionInfo()
