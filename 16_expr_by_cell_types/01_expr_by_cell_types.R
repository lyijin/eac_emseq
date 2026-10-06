#!/usr/bin/env Rscript

"> 01_expr_by_cell_types.R <

Post-CelFiE-ISH deconvolution, check whether the groupings differ in expression
patterns.

Use Z-scaled log10 (TPM+1) for visualisation.
- log10 crushes the long-tailed distribution into something normal-ish
- Z-scaling to see whether expression is higher in some cell type groups than
  others
" -> doc

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(forcats)
  library(GenomicRanges)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(readr)
  library(this.path)
  library(tibble)
  library(tidyr)
})


setwd(this.path::here())
HEATMAP_ANNOTS_FILE <- '../14_deconv_celfieish/heatmap_annots.tsv'
COMPILED_TPM_FILE <- '../data/udumanne_etal_deg/counts_tpm.csv.gz'
GENOMIC_ANNOT_FILE <- '../data/gencode.v43.annotation_2026-02-19.RData'

# error out if any of the files do not exist
stopifnot(file.exists(HEATMAP_ANNOTS_FILE))
stopifnot(file.exists(COMPILED_TPM_FILE))
stopifnot(file.exists(GENOMIC_ANNOT_FILE))

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
  'N'='NSq', 'B'='NDBE', 'L'='LGD', 'H'='HGD', 'C'='EAC')
disease_progression_vec <- c('NSq', 'NDBE', 'LGD', 'HGD', 'EAC')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# load pre-processed genomic annotations and functions to annotate GRanges
load(GENOMIC_ANNOT_FILE)
ensembl_to_genename_tib <-
  # one of the many variables in the RData file has a mapping of ENSEMBL
  # ENSGxxx IDs to gene names
  mcols(gencode_gene)[c('gene_id', 'gene_name')] |>
  # convert from DataFrame to tibble
  as_tibble() |>
  # remove version in "gene_id"
  mutate(gene_id=sub('\\..*', '', gene_id))


# load the heatmap annot file to get a list of methylation sample short IDs
heatmap_annot_tib <- read_tsv(HEATMAP_ANNOTS_FILE, show_col_types=FALSE) |>
  dplyr::rename(short_id=`...1`)
head(heatmap_annot_tib)


# time to read in the gene expression file, noting the sample IDs with RNA-seq data
tpm_tib <- read_csv(COMPILED_TPM_FILE, show_col_types=FALSE) %>%
  # sort by sample ID, except for the gene_id column
  select(gene_id, sort(colnames(.))) |>
  # remove version in "gene_id"
  mutate(gene_id=sub('\\..*', '', gene_id)) %>%
  # remove genes that are too lowly expressed
  filter(rowMeans(.[-1]) > 1) %>%
  # do a base10 log1p [i.e., log10 (x+1)] transform of raw TPM data for
  # correlation with methylation data
  mutate(across(!gene_id, ~ log1p(.x)/log(10)))
head(tpm_tib)
# confirm vast majority of the remaining values have transformed TPMs > 0
table(tpm_tib[-1] > 0)

# add in columns with no expression data into the expression data table
# (... yeah i know it's counterintuitive for now, but hang on)
short_ids_no_expr <- heatmap_annot_tib$short_id[
  !(heatmap_annot_tib$short_id %in% colnames(tpm_tib[-1]))]
tpm_tib <- tpm_tib |>
  # convert EMSEMBL IDs into gene names
  left_join(ensembl_to_genename_tib, by='gene_id') %>%
  # courtesy https://forum.posit.co/t/creating-new-na-columns-with-mutate-at-worked-in-dplyr-0-8-1-but-not-in-0-8-2/34443/6
  "[<-"(short_ids_no_expr, value=NA_real_) |>
  # then select samples with methylation data. those without expression data
  # will now just have NAs down the column
  select(c(gene_name, heatmap_annot_tib$short_id[order(heatmap_annot_tib$heatmap_col_order)]))
tpm_tib[1:6, 1:10]

# this tibble is now ready for rowwise filtering based on names of interesting
# genes, based on biological function
#
# define interesting gene lists!
goi_tib <- tibble()

# from naeini et al, nat commun, 2023
immune_hot_cold_markers <- c('CCL5', 'CD8A', 'NKG7', 'GZMA', 'GZMB', 'IDO1',
                             'SPP1', 'MMP3', 'SERPINB1', 'AQP3')
immune_hot_cold_markers %in% tpm_tib$gene_name    # should be all TRUE
goi_tib <- goi_tib |> bind_rows(
  tpm_tib |>
    filter(gene_name %in% immune_hot_cold_markers) |>
    arrange(factor(gene_name, levels=immune_hot_cold_markers))
)

# udumanne dysplasia signatures
dysplasia_markers <- c('SLC11A1', 'LUCAT1', 'MIR215', 'IL36A', 'RNU6-954P')
dysplasia_markers %in% tpm_tib$gene_name          # should be all TRUE
goi_tib <- goi_tib |> bind_rows(
  tpm_tib |>
    filter(gene_name %in% dysplasia_markers) |>
    arrange(factor(gene_name, levels=dysplasia_markers))
)


# plot expression of these genes in a heatmap. ordering of cols are from
# an upstream script
goi_df <- as.data.frame(goi_tib)
rownames(goi_df) <- goi_df$gene_name
goi_df <- as.matrix(goi_df[-1])

#+ fig.width=10, fig.height=4
p_expr <- pheatmap(
  goi_df,
  scale='row',
  border_color=NA,
  cluster_cols=FALSE,
  cluster_rows=FALSE,
  # these funky rows basically evaluates to c(17, 35, 57), which is where
  # the gaps are inserted to separate the different clusters
  gaps_col=cumsum(c(sum(heatmap_annot_tib$cluster_ID == 'Squamous'),
                    sum(heatmap_annot_tib$cluster_ID == 'Erythroid'),
                    sum(heatmap_annot_tib$cluster_ID == 'Columnar'))),
  # and same idea for the row gaps to produce "10"
  gaps_row=length(immune_hot_cold_markers),
  fontsize_col=8)
ggsave('raw_fig4c.pdf', plot=p_expr$gt, width=9.72, height=3.5)

# plot boxplot equivalents for each of the genes in the previous plot, to better
# eyeball trends
goi_long_tib <- goi_tib |>
  pivot_longer(!gene_name, names_to='short_id', values_to='logtpm') |>
  drop_na() |>
  # factoring to preserve relative order
  mutate(gene_name=fct_inorder(gene_name)) |>
  # add cluster ID annotations
  left_join(heatmap_annot_tib |> select(c(short_id, cluster_ID)), by='short_id') |>
  # factoring again for order. sigh ggplot why do you always sort alphabetically
  # modify the ordering to mimic disease progression, so that left vs. right
  # comparison of the supp fig is more visually similar
  mutate(cluster_ID=factor(
    cluster_ID, levels=c('Squamous', 'Columnar', 'Lymphoid', 'Erythroid'))) |>
  # for cell type vs disease classification eyeballing
  mutate(classification=factor(
    short_to_long_disease_vec[substr(short_id, 1, 1)], levels=disease_progression_vec))
head(goi_long_tib)

# calculate global mean (for dotted horizontal line in plots)
cluster_mean_tib <- goi_long_tib |>
  group_by(gene_name) |> 
  summarize(mean_logtpm=mean(logtpm))

#+ fig.width=10, fig.height=10
# plot expression by cell type cluster
g1 <- ggplot(goi_long_tib, aes(x=cluster_ID, y=logtpm, fill=cluster_ID)) +
  geom_boxplot(alpha=0.7) +
  geom_hline(data=cluster_mean_tib, aes(yintercept=mean_logtpm), linetype=2, alpha=0.5) +
  stat_compare_means(method='anova', label.y=Inf, vjust=1, size=3) +
  stat_compare_means(label='p.signif', method='t.test', ref.group='.all.',
                     hide.ns=TRUE, size=5) + # pairwise t vs. base mean
  scale_fill_manual(values=c('Squamous'='#8da0cb', 'Columnar'='#66c2a5',
                             'Lymphoid'='#e5c494', 'Erythroid'='#fc8d62')) +
  scale_x_discrete(guide=guide_axis(angle=45)) +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) +
  facet_wrap(~ fct_inorder(gene_name), ncol=3, scales='free_y') +
  labs(x='', y='log10 (TPM+1)') +
  theme_minimal(12) +
  theme(legend.position='top')

# repeat for disease classification
g2 <- ggplot(goi_long_tib, aes(x=classification, y=logtpm, fill=classification)) +
  geom_boxplot(alpha=0.7) +
  geom_hline(data=cluster_mean_tib, aes(yintercept=mean_logtpm), linetype=2, alpha=0.5) +
  stat_compare_means(method='anova', label.y=Inf, vjust=1, size=3) +
  stat_compare_means(label='p.signif', method='t.test', ref.group='.all.',
                     hide.ns=TRUE, size=5) + # pairwise t vs. base mean
  scale_fill_manual(values=disease_to_color_vec) +
  scale_x_discrete(guide=guide_axis(angle=45)) +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) +
  facet_wrap(~ fct_inorder(gene_name), ncol=3, scales='free_y') +
  labs(x='', y='log10 (TPM+1)') +
  theme_minimal(12) +
  theme(legend.position='top')

plot_grid(g1, g2, nrow=1, align='h')
ggsave('raw_suppfig9.pdf', width=11, height=11)

# for replicability purposes
sessionInfo()
