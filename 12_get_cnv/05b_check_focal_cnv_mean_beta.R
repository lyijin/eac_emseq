#!/usr/bin/env Rscript

"> 05b_check_focal_cnv_mean_beta.R <

Whole point of this script is to basically replot the heatmap in
`04_check_focal_cnv.R`, but overlay changes in methylation with up/down/sideways
arrows.
" -> doc

suppressPackageStartupMessages({
  library(Cairo)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(readr)
  library(stringr)
  library(this.path)
  library(tidyr)
})


setwd(this.path::here())
SUPP_FILES <- c('./raw_supptable23a.tsv', './raw_supptable23b.tsv')
FOCAL_CNV_MEAN_BETA_FILE <- './focal_cnv_mean_betas.tsv'

# error out if any of the files do not exist
stopifnot(file.exists(SUPP_FILES))
stopifnot(file.exists(FOCAL_CNV_MEAN_BETA_FILE))

# agreed upon colour scheme, based on Set1:
#   Normal squamous                 NSq    N   #984ea3   purple
#   Non-dysplastic Barrett's (IM)   NDBE   B   #377eb8   blue
#   LGD                             LGD    L   #4daf4a   green
#   HGD                             HGD    H   #a65628   brown
#   Esophageal adenocarcinoma       EAC    C   #e41a1c   red
disease_to_color_vec <- c(
  'N'='#984ea3', 'B'='#377eb8', 'L'='#4daf4a', 'H'='#a65628', 'C'='#e41a1c',
  'NSq'='#984ea3', 'NDBE'='#377eb8', 'LGD'='#4daf4a',
  'HGD'='#a65628', 'EAC'='#e41a1c')
short_to_long_disease_vec <- c(
  'N'='NSq', 'B'='NDBE', 'L'='LGD', 'H'='HGD', 'C'='EAC', HC='HGD/EAC')
disease_progression_vec <- c('NSq', 'NDBE', 'LGD', 'HGD', 'EAC')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# read the supplementary files (to replot the heatmap)
del_tib <- read_tsv(SUPP_FILES[1], show_col_types=FALSE)
del_tib <- del_tib |>
  filter(str_detect(gene, '^OR', negate=TRUE)) |>
  filter(str_detect(gene, '^POTE', negate=TRUE)) |>
  filter(str_detect(gene, '^IFNA', negate=TRUE))

amp_tib <- read_tsv(SUPP_FILES[2], show_col_types=FALSE)
amp_tib <- amp_tib |>
  filter(str_detect(gene, '^OR', negate=TRUE)) |>
  filter(str_detect(gene, '^POTE', negate=TRUE)) |>
  filter(str_detect(gene, '^IFNA', negate=TRUE))

mean_beta_tib <- read_tsv(FOCAL_CNV_MEAN_BETA_FILE, show_col_types=FALSE)
# easier to work with tibble in its transposed form
mean_beta_tib <- mean_beta_tib |>
  pivot_longer(-`...1`, names_to='gene') |>
  pivot_wider(names_from=`...1`, values_from=value)
mean_beta_tib

# i have visually confirmed the genes with NAs REALLY REALLY REALLY had
# no methylation data at all in the 20 kb window where the focal deletion/
# amplification happened
#
# not much i can do with NAs; drop em
mean_beta_tib <- mean_beta_tib |> drop_na() |>
  # remove # CpGs per gene, was there to confirm why NAs were NAs initially
  # (cus there were no CpGs, basically)
  mutate(gene=str_replace(gene, pattern=' \\| .*$', replacement='')) |>
  # drop N samples, they have no CNVs (repurposed python script did calculate
  # mean betas for N, was too lazy to alter the script to not do it)
  select(-starts_with('N')) |> select(-meanN)

# finally, the fun starts. compare per-sample, per-gene mean beta vs. group
# mean. if mean beta is 2% above group mean, replace with up arrow; 2% below,
# down arrow; within 2%, sideways arrow
arrows_tib <- mean_beta_tib %>%
  mutate(across(
    starts_with('B'), ~ case_when(
      .x > meanB + 0.02 ~ '\u2191',
      .x < meanB - 0.02 ~ '\u2193',
      .default = '\u2194'),
    .names = "{.col}")) %>%
  mutate(across(
    starts_with('L'), ~ case_when(
      .x > meanL + 0.02 ~ '\u2191',
      .x < meanL - 0.02 ~ '\u2193',
      .default = '\u2194'),
    .names = "{.col}")) %>%
  mutate(across(
    starts_with('H'), ~ case_when(
      .x > meanH + 0.02 ~ '\u2191',
      .x < meanH - 0.02 ~ '\u2193',
      .default = '\u2194'),
    .names = "{.col}")) %>%
  mutate(across(
    starts_with('C'), ~ case_when(
      .x > meanC + 0.02 ~ '\u2191',
      .x < meanC - 0.02 ~ '\u2193',
      .default = '\u2194'),
    .names = "{.col}")) |>
  select(-starts_with('mean'))

# visual check that code went alright
mean_beta_tib |> select(c(gene, B02, meanB, C02, meanC)) |> head()
arrows_tib |> select(c(gene, B02, C02)) |> head()
# yay


# plot focal deletion heatmap WITH ARROWS YAY
del_df <- del_tib |>
  select(matches('^[A-Z]', ignore.case=FALSE)) |>
  as.matrix()
rownames(del_df) <- del_tib$gene

# create matching df for arrows
del_arrows_tib <- left_join(
  del_tib |> select(gene),
  arrows_tib,
  by='gene') %>%
  replace(is.na(.), '')
del_arrows_df <- del_arrows_tib |>
  select(-gene) |>
  as.data.frame()
rownames(del_arrows_df) <- del_arrows_tib$gene

#+ fig.width=8, fig.height=3
pd <- pheatmap(
  del_df,
  color=c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'),
  border_color=NA,
  breaks=c(0, 0.8, 1.6, 2.4, 3.2, 4),
  cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize=7,
  display_numbers=del_arrows_df,
  gaps_col=c(24,36,47))
ggsave('raw_suppfig4a.pdf', plot=pd$gtable, device=cairo_pdf, width=6, height=2)


# repeat for amps
amp_df <- amp_tib |>
  select(matches('^[A-Z]', ignore.case=FALSE)) |>
  as.matrix()
rownames(amp_df) <- amp_tib$gene

# create matching df for arrows
amp_arrows_tib <- left_join(
  amp_tib |> select(gene),
  arrows_tib,
  by='gene') %>%
  replace(is.na(.), '')
amp_arrows_df <- amp_arrows_tib |>
  select(-gene) |>
  as.data.frame()
rownames(amp_arrows_df) <- amp_arrows_tib$gene

#+ fig.width=8, fig.height=3
pa <- pheatmap(
  amp_df,
  color=c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'),
  border_color=NA,
  breaks=c(0, 0.8, 1.6, 2.4, 3.2, 4),
  cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize=7,
  display_numbers=amp_arrows_df,
  gaps_col=c(24,36,47))
ggsave('raw_suppfig4b.pdf', plot=pa$gtable, device=cairo_pdf, width=5.91, height=1.9)
# manually jiggle width value so that the heatmaps of "pd" and "pa" have nearly
# identical plot areas for the actual heatmap


# for replicability purposes
sessionInfo()
