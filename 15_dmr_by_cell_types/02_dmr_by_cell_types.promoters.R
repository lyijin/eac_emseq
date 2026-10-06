#!/usr/bin/env Rscript

"> 02_dmr_by_cell_types.promoters.R <

Post-CelFiE-ISH deconvolution, check whether the groupings differ in methylation
patterns.

There are quite a number of DMRs linked with disease progression, but analyses
typically group samples by disease classifications. Would shifting to a
cell-type-based approach lead to a different view?
" -> doc

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(GenomicRanges)
  library(pheatmap)
  library(readr)
  library(readxl)
  library(stringr)
  library(this.path)
  library(tibble)
  library(tidyr)
})

setwd(this.path::here())
HEATMAP_ANNOTS_FILE <- '../14_deconv_celfieish/heatmap_annots.tsv'
COMPILED_BETA_FILE <- '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
DMR_FOLDER <- './filtered_dmrs/'
LIT_METH_GENES <- '../data/supp_table_25.xlsx'

# error out if any of the files do not exist
stopifnot(file.exists(HEATMAP_ANNOTS_FILE))
stopifnot(file.exists(COMPILED_BETA_FILE))
stopifnot(file.exists(DMR_FOLDER))
stopifnot(file.exists(LIT_METH_GENES))

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
  'N'='NSq', 'B'='NDBE', 'L'='LGD', 'H'='HGD', 'C'='EAC', 'HC'='HGD/EAC')
disease_progression_vec <- c('NSq', 'NDBE', 'LGD', 'HGD', 'EAC')
comparison_order <- c('NDBE vs. NSq', 'LGD vs. NDBE', 'HGD vs. LGD', 'EAC vs. HGD')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# load the heatmap annot file to get a list of methylation sample short IDs
heatmap_annot_tib <- read_tsv(HEATMAP_ANNOTS_FILE, show_col_types=FALSE) |>
  dplyr::rename(short_id=`...1`) |>
  mutate(cluster_ID=factor(
    cluster_ID, levels=c('Squamous', 'Erythroid', 'Columnar', 'Lymphoid')))
head(heatmap_annot_tib)


# focus on hypermethylated (higher meth in more serious pathology) markers
dmr_files <- c(paste0(DMR_FOLDER, 'Columnar_vs_Squamous.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'Lymphoid_vs_Squamous.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'Erythroid_vs_Squamous.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'Lymphoid_vs_Columnar.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'Erythroid_vs_Columnar.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'Erythroid_vs_Lymphoid.hypermeth.tsv.gz'))
stopifnot(file.exists(dmr_files))

# create a dt to store annots of DMRs
annot_dmr_tib <- tibble(
  seqnames=character(), start=numeric(), end=numeric(), width=numeric(),
  qval=numeric(), delta_beta=numeric(), nCpG=numeric(), mean_g1=numeric(),
  mean_g2=numeric(), gene_id=character(), gene_name=character(),
  promoter=numeric(), genebody=numeric(), intergenic=numeric(),
  comparison=character(), hyper_or_hypo=character())

for (dmr_file in dmr_files) {
  dmr_tib <- read_tsv(dmr_file, show_col_types=FALSE)
  diag_message(dmr_file, ' has ', formatC(nrow(dmr_tib), big.mark=','), ' rows.')
  
  # filter for promoters! only pick DMRs that sit COMPLETELY WITHIN promoter
  # regions (1.0). partial overlaps result in fractional values
  dmr_tib <- dmr_tib |> filter(promoter == 1)
  diag_message(dmr_file, ' has ', formatC(nrow(dmr_tib), big.mark=','), ' promoter rows.')
  
  # parse filename for contextual clues for sample types / risk states
  fname_split <- unlist(strsplit(dmr_file, '/', fixed=TRUE))
  fname_split <- unlist(strsplit(fname_split[length(fname_split)], '.', fixed=TRUE))
  hyper_or_hypo <- fname_split[2]
  compared_states <- fname_split[1]
  compared_states <- unlist(strsplit(compared_states, '_vs_'))
  compared_states <- paste0(compared_states[1], ' vs. ', compared_states[2])
  diag_message(dmr_file, ' are ', hyper_or_hypo, ' markers for ', compared_states)
  
  dmr_tib$comparison <- compared_states
  dmr_tib$hyper_or_hypo <- hyper_or_hypo
  
  # rbind important columns into `annot_dmr_tib`
  annot_dmr_tib <- bind_rows(annot_dmr_tib, dmr_tib[colnames(annot_dmr_tib)])
}
annot_dmr_tib <- annot_dmr_tib |>
  # remove version in "gene_id"
  mutate(gene_id=sub('\\..*', '', gene_id)) |>
  # convert comparison to a factor
  mutate(comparison=as_factor(comparison))
head(annot_dmr_tib)

# some genes have multiple DMRs in its promoter--which is fine, but just
# remember that the analysis is on a per-gene (specifically, ENSEMBL ID) basis,
# and for each sample, the mean beta will be calculated from all CpGs across
# the multiple DMRs
annot_dmr_tib |> filter(comparison == 'Erythroid vs. Squamous') |> 
  filter(hyper_or_hypo == 'hypermeth') |> arrange(gene_id) |> head(n=10)


# subselect interesting DMRs now--letting the script run on ALL DMRs is TOO
# SLOW. but in doing so, blocks exploratory visualisation on genes that aren't
# in this pre-selected subset
#
# CIMP gene list from Krause et al., Carcinogenesis, 2016
cimp_meth_genes <- c('MLH1', 'CDKN2A', 'MGMT', 'CACNA1G', 'IGF2', 'NEUROG1',
                     'RUNX3', 'SOCS1', 'KCNK13', 'SLIT1', 'RAB31', 'FOXL2',
                     'B3GAT2', 'FAM78A', 'MYOCD', 'KCNC1', 'FSTL1', 'SLC6A4')
# other genes from supplementary table
interesting_meth_genes_tib <- read_excel(
  LIT_METH_GENES, sheet=1, range='A13:M133') |>  # unfortunately readxl doesn't allow
  select(-c(3:10)) |>                            # selection of cols A:B & K:M
  drop_na(Gene) |>                        # skips blank lines with no genes
  filter_out(str_detect(Gene, "'")) |>  # remove genes with "'"
  filter_out(str_detect(`B vs. N`, '(control)')) |>  # remove control genes
  filter_out(str_detect(`HC vs. N`, '(control)')) |>
  filter_out(str_detect(`HC vs. B`, '(control)')) |>
  mutate(Gene=str_replace(Gene, '_.*', '')) |>  # remove anything after '_'
  mutate(`HC vs. B`=str_replace(`HC vs. B`, ' \\(.*', '')) |>  # remove parens
  mutate(Publication=str_replace(Publication, ' \\(.*', ''),  # remove parens
         `B vs. N`=as.logical(`B vs. N`),                     # change chr to lgl
         `HC vs. N`=as.logical(`HC vs. N`),
         `HC vs. B`=as.logical(`HC vs. B`))
interesting_meth_genes_tib

# for comparison, spike in 9 genes with the largest delta betas with biomarker
# potential to illustrate that EM-seq gives way more breadth in discovery work.
# i.e., the truth is out there! *cue X-Files music*
p1 <- annot_dmr_tib |>
  filter(comparison == 'Lymphoid vs. Squamous') |>
  slice_max(order_by=delta_beta, n=3) |>
  pull(gene_name)

p2 <- annot_dmr_tib |>
  filter(comparison == 'Lymphoid vs. Columnar') |>
  slice_max(order_by=delta_beta, n=4) |>
  filter(gene_name != 'PROM2') |>  # lower in Squamous than Columnar, rejected
  pull(gene_name)

p3 <- annot_dmr_tib |>
  filter(comparison == 'Erythroid vs. Columnar') |>
  slice_max(order_by=delta_beta, n=6) |>
  filter(!gene_name %in% c('IFFO1', 'UNKL', 'SYNE2')) |>  # meth not lowest in Columnar
  pull(gene_name)

# a bit manual in terms of gene selection, but as the pairwise comparisons are
# not aware of the overall methylation patterns in the other non-compared
# states, sometimes the baseline has higher methylation than a non-compared
# group. kick those genes out manually if they do
proposed_genes <- c(p1, p2, p3)


# define commercial utility on a gene-by-gene basis. e.g., CDKN2A was used in
# Barrett 1999 for B vs. N and HC vs. N, but Eads 2001 also thinks it's good
# for HC vs. B too. this means that for CDKN2A, it has utility in all three
# mentioned comparisons. in programming terms--use any() to determine whether
# TRUE exists anywhere for the gene across publications
commercial_purpose_tib <- interesting_meth_genes_tib |>
  select(-Publication) |>
  group_by(Gene) |>
  summarize(`NDBE vs NSq`=any(`B vs. N`),
            `HGD/EAC vs. NSq`=any(`HC vs. N`),
            `HGD/EAC vs. NDBE`=any(`HC vs. B`))

# for plotting reasons, don't let genes repeat
interesting_meth_genes_tib <- bind_rows(
  tibble(Gene=cimp_meth_genes, Publication='Krause 2016'),
  interesting_meth_genes_tib,
  tibble(Gene=proposed_genes, Publication='Proposed')) |>
  # first publication to mention it gets the honour of having the gene
  # assigned to it 
  distinct(Gene, .keep_all=TRUE) |>
  select(Gene, Publication) |>
  # merge in the commercial utility tibble
  left_join(commercial_purpose_tib, by='Gene')
interesting_meth_genes_tib
interesting_meth_genes <- interesting_meth_genes_tib |> pull(Gene)

# do the DMR subsetting, finally
annot_dmr_tib <- annot_dmr_tib |> filter(gene_name %in% interesting_meth_genes)

# see the annoying thing about meth data is that it's done at a loci-level, not
# like transcripts or proteins where values can be assigned to genes
# straightaway. the loci-based analysis does lead to a gene (e.g., AKAP2, below)
# having a promoter region with hypo- and hypermeth DMRs simultaneously; and
# A vs. Squamous and B vs. Squamous both having DMRs that aren't exactly
# overlapping (raises the question, which to use for plotting purposes, as it's
# easier to assign loci to genes and talk about biological function at the gene
# level)
annot_dmr_tib |> arrange(gene_name) |> head(n=10)
# so! tiebreak rules:
#   1. focus on hypermethylated DMRs (as most assays are designed for that)
#   2. pick the region with the lowest qval. usually strikes a good balance
#      between "width"/"nCpG" and "delta_beta". as in, some regions can have
#      extreme delta_beta values when windows are small, but qval will be higher
#   3. if qvals are tied (especially when q=0), pick the region with higher
#      delta_beta. easier to design assays on loci with larger changes
annot_dmr_tib <- annot_dmr_tib |>
  filter(hyper_or_hypo == 'hypermeth') |>
  group_by(gene_name) |>
  slice_min(order_by=qval, n=1) |>
  slice_min(order_by=-delta_beta, n=1)
annot_dmr_tib |> arrange(gene_name)

# defensive programming: all genes should be unique after tiebreaking
stopifnot(length(unique(annot_dmr_tib$gene_name)) == nrow(annot_dmr_tib))

# create a GRanges object for overlapping purposes later
dmr_gr <- makeGRangesFromDataFrame(
  annot_dmr_tib, keep.extra.columns=TRUE, ignore.strand=TRUE)


# read in the chonky per-CpG beta values table
heatmap_col_order <- heatmap_annot_tib$short_id[order(heatmap_annot_tib$heatmap_col_order)]
beta_tib <- read_tsv(COMPILED_BETA_FILE, show_col_types=FALSE, lazy=TRUE) |>
  # sort the columns according to the heatmap order, but start with chr/start/end
  select(chr, start, end, all_of(heatmap_col_order))


# convert beta tibble to a GRanges object, then to reduce memory usage, retain
# CpG positions that are in at least one DMR
beta_gr <- makeGRangesFromDataFrame(
  beta_tib,
  keep.extra.columns=TRUE,
  ignore.strand=TRUE)
# note: select='first' because i care that a CpG is in at least one DMR (but not
# which DMR exactly). computationally quicker than the default of select='all'
fo_beta_dmr <- findOverlaps(beta_gr, dmr_gr, minoverlap=2L, select='first')
beta_gr <- beta_gr[!is.na(fo_beta_dmr), ]
beta_gr[, 1:5]  # visual check

# nuke the tibble to free memory
rm(beta_tib)

# i'm sure there's an arcane apply() that achieves the same outcome, but...
# for each comparison...
mean_beta_tib <- NULL
for (co in levels(annot_dmr_tib$comparison)) {
  # for each hyper/hypo...
  for (hoh in unique(annot_dmr_tib$hyper_or_hypo)) {
    # for every gene with DMR(s)...
    for (gi in annot_dmr_tib |> filter(comparison == co, hyper_or_hypo == hoh) |> pull(gene_id) |> unique()) {
      gi_dmr_gr <- dmr_gr[dmr_gr$gene_id == gi, ]
      fo_beta_gidmr <- findOverlaps(beta_gr, gi_dmr_gr,
                                    minoverlap=2L, select='first')
      mean_betas <- colMeans(as.matrix(mcols(beta_gr[!is.na(fo_beta_gidmr), ])), na.rm=TRUE)
      
      # aggregate per-sample mean betas by cell type (Squamous, Erythroid,
      # Columnar, Lymphoid; in that specific order)
      celltype_mean_betas <- sapply(
        levels(heatmap_annot_tib$cluster_ID),
        function(x) {mean(mean_betas[
          heatmap_annot_tib |> filter(cluster_ID == x) |> pull(short_id)], na.rm=TRUE)},
        USE.NAMES=TRUE)
      
      # append results into a single-row tibble first
      mean_beta_row <- bind_cols(
        # descriptive stuff & cor results
        bind_rows(
          gene_id=gi,
          comparison=co,
          hyper_or_hypo=hoh,
        ),
        # add in cell type beta values
        bind_rows(celltype_mean_betas),
        # keep beta values
        bind_rows(mean_betas) %>% rename_with(~ paste0('beta_', .))
      )
      
      # then append the row into the overall tibble
      mean_beta_tib <- bind_rows(mean_beta_tib, mean_beta_row)
    }
  }
}

# merge DMR results with per-sample beta values
mean_beta_tib <- left_join(
  annot_dmr_tib |> select(-c(promoter, genebody, intergenic)),  # drop unimpt cols
  mean_beta_tib,
  by=c('gene_id', 'comparison', 'hyper_or_hypo')) |>
  ungroup() |>  # grouped tibbles has unexpected behaviour with as.matrix()
  relocate(gene_name, gene_id)  # move gene identifiers to the front
# then sort by the order of plotted genes
mean_beta_tib <- mean_beta_tib |>
  mutate(gene_name=factor(
    gene_name, levels=interesting_meth_genes[interesting_meth_genes %in% mean_beta_tib$gene_name])) |>
  arrange(gene_name)
head(mean_beta_tib)

# save this table as a supp table
write_tsv(mean_beta_tib, 'raw_supptable26.tsv')


# for pheatmap: convert tibble into matrix for plotting
goi_df <- mean_beta_tib |> select(starts_with('beta_')) |> as.matrix()
colnames(goi_df) <- sub('^beta_', '', colnames(goi_df))
rownames(goi_df) <- mean_beta_tib$gene_name
goi_df[1:6, 1:6]

# for pheatmap: create row annot for published utility of meth markers
annot_df <- interesting_meth_genes_tib |>
  select(-Publication) %>%
  mutate(across(where(is.logical), ~na_if(., FALSE))) |>  # there are 2 FALSEs: convert to NA
  mutate_all(as.character) |>  # ... yeah, pheatmap hates logicals
  column_to_rownames(var='Gene')

# deduce "gaps_row" for pheatmap programmatically, instead of manually figuring
# things out
pub_cumsum <- interesting_meth_genes_tib |> 
  filter(Gene %in% rownames(goi_df)) |>          # whitelist plotted genes
  mutate(Publication=as_factor(Publication)) |>  # prevent count() from sorting alphabetically
  count(Publication, sort=FALSE) |>
  pull(n) |>  # an array: 13  5  1  1  1 ...
  cumsum()    # an array: 13 18 19 20 21 ...
# split proposed genes from a block of 9 to three blocks of 3
pub_cumsum <- sort(c(pub_cumsum, max(pub_cumsum)-6, max(pub_cumsum)-3))
pub_cumsum

#+ fig.width=10, fig.height=12
p_meth <- pheatmap(
  goi_df,
  #scale='row',  # meth data, better to plot betas 0-1 unscaled
  border_color=NA,
  cluster_cols=FALSE,
  cluster_rows=FALSE,
  # these funky rows basically evaluates to c(17, 35, 57), which is where
  # the gaps are inserted to separate the different clusters
  gaps_col=c(sum(heatmap_annot_tib$cluster_ID == 'Squamous'),
             sum(heatmap_annot_tib$cluster_ID %in% c('Squamous', 'Erythroid')),
             sum(heatmap_annot_tib$cluster_ID %in% c('Squamous', 'Erythroid', 'Columnar'))),
  # introduce gaps for proposed genes (9 in total, 3/3/3)
  gaps_row=pub_cumsum,
  annotation_row=annot_df,
  annotation_legend=FALSE,
  fontsize_col=8)
ggsave('raw_fig4b.pdf', plot=p_meth$gt, width=10.73, height=11)


# for replicability purposes
sessionInfo()
