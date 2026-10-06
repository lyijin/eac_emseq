#!/usr/bin/env Rscript

"> 01_celfieish_deconvs_unsupervised.R <

Post-CelFiE-ISH deconvolution, do the differences in tissue compositions
predict pathology or progression across the EAC spectrum? Run unsupervised
analyses (PCA, hierarchical clustering) to tease out higher-order patterns.
" -> doc

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(GGally)
  library(ggplot2)
  library(pheatmap)
  library(readr)
  library(this.path)
  library(tibble)
  library(tidyr)
})


setwd(this.path::here())
COMPILED_DECONV_FILE <- './all_86_samples.celfie-ish_deconv.tsv'
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'

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

# function to plot PCA
run_prcomp <- function(df, most_variable_rows='all') {
  # run prcomp to generate plottable pca df
  
  # remove rows that have 0 variance
  df <- df[!apply(df, 1, var, na.rm=TRUE) == 0, ]
  
  # by default, uses all rows in df. if given a numeric value, selects said
  # number of rows in df
  if (most_variable_rows != 'all' && most_variable_rows == as.integer(most_variable_rows)) {
    df <- df[order(apply(df, 1, var, na.rm=TRUE), decreasing=TRUE)[1:most_variable_rows], ]
  }
  
  prcomp_obj <- prcomp(t(df), scale.=FALSE, rank=20)  # note scale is FALSE!
  prcomp_obj
}

plot_pca <- function(prcomp_obj, x='PC1', y='PC2',
                     size=NULL, alpha=0.3, color_palette=NULL, risk_vector=FALSE,
                     hide_labels=FALSE, title=NULL, subtitle=FALSE, hide_legend=FALSE) {
  # extract important stuff from the prcomp object
  pca_coords <- as.data.frame(prcomp_obj$x)
  eigs <- prcomp_obj$sdev ^ 2
  
  x_pct <- round(eigs[as.numeric(substr(x, 3, nchar(x)))] / sum(eigs) * 100, 2)
  y_pct <- round(eigs[as.numeric(substr(y, 3, nchar(y)))] / sum(eigs) * 100, 2)
  disease_class <- unname(short_to_long_disease_vec[substr(rownames(pca_coords), 1, 1)])
  disease_class <- factor(disease_class, levels=disease_progression_vec)
  g <- ggplot(pca_coords, aes(x=!!sym(x), y=!!sym(y)))
  
  if (isFALSE(risk_vector)) {
    g <- g + geom_point(aes(color=disease_class), 
                        size=if (is.null(size)) 4 else size, alpha=alpha)
  } else {
    # risk_vector is assumed to be an array containing four risk labels:
    #   H_adj
    #   H_NA
    #   H_prog
    #   L_NA, but NOT the name of the column that stores risk
    risk_overall <- substr(risk_vector, 1, 1)
    risk_detail <- substr(risk_vector, 3, 100)
    g <- g + 
      geom_point(aes(color=disease_class, shape=risk_detail, alpha=risk_overall),
                 size=if (is.null(size)) 4 else size) +
      scale_shape_manual(values=c(15, 16, 17)) +  # NA gets circle
      scale_alpha_discrete(range=c(0.5, 0.3)) +   # H gets more opaque
      labs(shape='Specific risk', alpha='Overall risk')
  }
  
  if (isFALSE(hide_labels)) {
    g <- g + geom_text(aes(label=rownames(pca_coords)), check_overlap=TRUE, hjust='inward', vjust='inward')
  }
  
  g <- g +
    labs(x=paste0(x, ' (', x_pct, '%)'),
         y=paste0(y, ' (', y_pct, '%)'),
         color='Classification') +
    theme_minimal(12)
  
  if (!is.null(color_palette)) {
    g <- g + scale_color_manual(values=color_palette)
  }
  
  if (!is.null(title)) {
    g <- g + labs(title=title)
  }
  
  if (subtitle) {
    # only requires TRUE/FALSE here, as subtitles just show number of positions
    g <- g + labs(subtitle=paste0('# of tissues: ', format(nrow(prcomp_obj$rotation), big.mark=',')))
  }
  
  if (hide_legend) {
    g <- g + theme(legend.position='none')
  }
  
  g
}


# subset clinical data for these included samples
clin_dt <- read_tsv(CLINDATA_FILE, show_col_types=FALSE) |> filter(emseq_include)
# change some labels to the current standard terminology
clin_dt <- clin_dt |> mutate(classification = case_when(
  classification == 'Normal squamous' ~ 'NSq',
  classification == 'IM' ~ 'NDBE',
  classification == 'Cancer' ~ 'EAC',
  TRUE ~ classification)) |>
  arrange(short_id)
# distinguish high-risk --> adjacent / progression / NA (for cancers)
clin_dt$risk_type <- paste(clin_dt$risk, clin_dt$H_risk_category, sep='_')

# read the compiled file & convert "emseq_sample_no" in `deconv_dt` to "short_id"
short_id_dt <- clin_dt |> select(emseq_sample_no, short_id)
deconv_dt <- read_tsv(COMPILED_DECONV_FILE, show_col_types=FALSE)
deconv_dt <- left_join(deconv_dt, short_id_dt, by='emseq_sample_no') |>
  relocate(short_id) |>
  select(-emseq_sample_no) |>
  arrange(short_id)
deconv_dt 

# sanity check
stopifnot(deconv_dt$short_id == clin_dt$short_id)  # sanity check to confirm equality


# for PCA/heatmaps, samples have to occupy columns (not rows)
# tidyr isn't a fan of transpose, hence the hacky code
#   https://github.com/tidyverse/tidyr/issues/925
#   https://stackoverflow.com/questions/28917076/transposing-data-frames-using-the-tidyverse/28917212
deconv_dt <- deconv_dt |>
  pivot_longer(-1) |>
  pivot_wider(names_from=1, values_from=value)


#+ fig.width=10, fig.height=10
# plots to check clustering across samples based on deconvoluted tissue
# proportions
prcomp_obj <- run_prcomp(deconv_dt |> select(-name))
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of CelFiE-ISH on', ncol(deconv_dt), 'samples'),
         subtitle=TRUE)

#+ fig.width=10, fig.height=10
# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)
# hmm, not really

# layer in risk information to see whether high-risk vs. low-risk separates on 
# any of the principal components
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type,
         title=paste('PCA of CelFiE-ISH on', ncol(deconv_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are risk-related trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)


# plot unsupervised clustering of heatmap
pheatmap(deconv_dt |> column_to_rownames(var='name'),
         clustering_method='ward.D2')
pheatmap(deconv_dt |> column_to_rownames(var='name'), scale='column',
         clustering_method='ward.D2')
# hmm, quite a number of uninteresting tissues that are low abundance/variance
# that is clagging the plot up

heatmap_dt <- deconv_dt %>% mutate(rowmeans=rowMeans(select(., -name))) 
summary(heatmap_dt$rowmeans)  # 50% of tissues have mean composition of <= 0.34%
heatmap_dt <- heatmap_dt |>
  filter(rowmeans > 1) |>     # set rowMeans threshold at 1%
  select(-rowmeans)
heatmap_df <- heatmap_dt |> column_to_rownames(var='name')
heatmap_df

# plot heatmap with sub-1% tissues removed
annot_dt <- clin_dt |>
  select(short_id, classification, gender, age_at_collection, risk) |>
  rename(age=age_at_collection)
annot_df <- annot_dt |> column_to_rownames(var='short_id')
color_list <- list(
  classification=disease_to_color_vec[6:10],
  gender=c('M'='#67a9cf', 'F'='#ef8a62'),  # RdBu color scheme
  risk=c('H'='#e9a3c9', 'L'='#a1d76a'),    # PiYG color scheme
  cluster_ID=c('Columnar'='#66c2a5', 'Squamous'='#8da0cb',
               'Lymphoid'='#e5c494', 'Erythroid'='#fc8d62')) # Set2 color scheme
# plot with unscaled percentage values
#+ fig.width=10, fig.height=6
p_unscaled <- pheatmap(heatmap_df,
                       border_color=NA,
                       clustering_method='ward.D2',
                       cutree_cols=4, cutree_rows=5,
                       annotation_col=annot_df,
                       annotation_colors=color_list,
                       fontsize_col=8)
# the bad thing about the unscaled plot is that head-neck-ep is overwhelmingly
# present in the N samples (bad as in, from a visualisation perspective),
# which stretches the red/blue colour scale out and relegates most values to
# a sea of blue.
#
# whinging aside, note that the overall tissue groupings (split into five
# groups) are consistent with the groupings from the scaled plot. to me,
# this lends more credence that the scaled plot paints a similar picture of the
# data as the unscaled one (but with nicer colours)

# plot using scaled values (z-scores)
p_scaled <- pheatmap(heatmap_df, scale='column',
                     border_color=NA,
                     clustering_method='ward.D2',
                     cutree_cols=4, cutree_rows=5,
                     annotation_col=annot_df,
                     annotation_colors=color_list,
                     fontsize_col=8)
# looks quite interesting & biologically relevant. based on the four-group
# clustering of the samples, propagate this clustering to all subsequent plots.
# to do so, cluster ID has to be saved into "annot_df"
annot_df$cluster_ID <- cutree(p_scaled$tree_col, k=4)

# names of clusters are given post-hoc
annot_df$cluster_ID[annot_df$cluster_ID == 1] <- 'Columnar'
annot_df$cluster_ID[annot_df$cluster_ID == 2] <- 'Squamous'
annot_df$cluster_ID[annot_df$cluster_ID == 3] <- 'Lymphoid'
annot_df$cluster_ID[annot_df$cluster_ID == 4] <- 'Erythroid'
annot_df$cluster_ID <- factor(annot_df$cluster_ID, levels=c('Columnar', 'Squamous', 'Lymphoid', 'Erythroid'))

# replot with cluster ID, and save the publication-grade viz
p_scaled <- pheatmap(heatmap_df, scale='column',
                     border_color=NA,
                     clustering_method='ward.D2',
                     cutree_cols=4, cutree_rows=5,
                     annotation_col=annot_df,
                     annotation_colors=color_list,
                     fontsize_col=8)
ggsave('raw_fig4a.pdf', plot=p_scaled$gt, width=12, height=5)
# lots of biological nuggets here: N-like cluster 2 dominated by head-neck-ep
# (as expected), C-like cluster 4 dominated by erythocyte progenitors (maybe
# expected?) and seemingly has an inverse correlation with T cell levels too.
# clusters 1 and 3 are a bit handwavy. 1 has more (columnar?) epithelia,
# 3 has higher immune presence (B cells / T cells / granulocytes).

# save the clustering and relative ordering for downstream plots to mimic the
# same order--provide qualitative check on whether clusters differ in terms of
# meth/expr patterns (i.e., any interesting biology?)
annot_df$heatmap_col_order <- order(p_scaled$tree_col$order, rownames(annot_df))
write.table(annot_df, file='heatmap_annots.tsv', quote=FALSE, sep='\t', col.names=NA)
# `rownames(annot_df)[order(annot_df$heatmap_col_order)]` recapitulates the
# ordering of the heatmap


# plot a pairplot for the 7 tissues that contributed the most to the clustering
heatmap_df
pp_df <- as.data.frame(t(heatmap_df))
pp_df$classification <- unname(short_to_long_disease_vec[substr(rownames(pp_df), 1, 1)])
#+ fig.width=8, fig.height=8
pp_df |> as_tibble() |>
  select('classification', 'Blood-B', 'Blood-T', 'Gastric-Ep', 'Colon-Ep',
         'Small-Int-Ep', 'Head-Neck-Ep', 'Eryth-prog') %>%
  ggpairs(., columns=2:ncol(.),
          lower=list(continuous=wrap('smooth', alpha=0.3)),
          diag='blank') +
  theme_minimal()
ggsave('raw_suppfig8.pdf', width=7, height=7)


# plot correlation matrix
#+ fig.width=10, fig.height=9
cor_mat <- cor(deconv_dt |> select(-name), method='pearson')
diag(cor_mat) <- NA
p_cor <- pheatmap(cor_mat,
                  color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
                  border_color=NA,
                  clustering_method='ward.D2',
                  cutree_cols=5, cutree_rows=5,
                  annotation_col=annot_df,
                  annotation_colors=color_list,
                  fontsize=7,
                  main=paste0('Clustered heatmap of pairwise Pearson r values (# of tissues: ', 
                              format(nrow(deconv_dt), big.mark=','), ')'))
# this plot visualises the consistency of the clusters from z-score scaling vs.
# clusters from all pairwise pearson correlations. 4 of 5 clusters show
# remarkable consistency; last one is a bit of a "catch-all" bin


# what if we excluded normals as well, how does IM -> cancer samples look on
# a PCA?
#+ fig.width=10, fig.height=10
sliced_deconv_dt <- deconv_dt |> select(-starts_with('N', ignore.case=FALSE))
sliced_clin_dt <- clin_dt[clin_dt$short_id %in% colnames(sliced_deconv_dt), ]
prcomp_obj <- run_prcomp(sliced_deconv_dt |> select(-name))
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of CelFiE-ISH on non-N', ncol(sliced_deconv_dt), 'samples'), subtitle=TRUE) +
  theme(legend.position='inside',
        legend.position.inside=c(0.05, 0.05),
        legend.justification=c(0, 0),
        legend.box.background=element_rect(fill='#ffffff80', color='#33333366'))
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type,
         title=paste('PCA of CelFiE-ISH on non-N', ncol(sliced_deconv_dt), 'samples'), subtitle=TRUE)
# no clear trends, B are a bit more clustered on the RHS of PC1, but otherwise
# not really

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)
# meh.


# similar heatmap for non-N
#+ fig.width=10, fig.height=9
rowname_subset <- rownames(cor_mat)[!grepl('^N', rownames(cor_mat))]
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         border_color=NA,
         clustering_method='ward.D2',
         cutree_cols=5, cutree_rows=5,
         annotation_col=annot_df,
         annotation_colors=color_list,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values for non-N ',
                     '(# of tissues: ', format(nrow(deconv_dt), big.mark=','), ')'))
# pretty similar to the with-N plot; the N-like cluster 2 shrinks to the card-
# carrying B members, other bins echo previous plot


# plot correlation matrices for each contiguous status separately (N & B; B & L;
# L & H; H & C). hopefully see patterns in transitioning from less serious
# classification to more serious classification. also pay attention to the risk
# labelling. are there blocks of L/H within existing clusters?
cor_mat <- cor(deconv_dt |> select(-name), method='pearson')
diag(cor_mat) <- NA

# plot per-category heatmaps
rowname_subset <- grepl('^N|^B', rownames(cor_mat))
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         border_color=NA,
         clustering_method='ward.D2',
         cutree_cols=5, cutree_rows=5,
         annotation_col=annot_df,
         annotation_colors=color_list,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values for N/B ', 
                     '(# of tissues: ', format(nrow(deconv_dt), big.mark=','), ')'))
# hmm. if we go by cluster IDs, looking at B and % high/low risk
# - cluster 2 (N-like): 1/6 high risk (fisher's exact p=0.35)
# - cluster 1 (squamous): 5/13 high risk
# - cluster 3 (high bloods): 6/10 high risk

rowname_subset <- grepl('^B|^L', rownames(cor_mat))
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         border_color=NA,
         clustering_method='ward.D2',
         cutree_cols=4, cutree_rows=4,
         annotation_col=annot_df,
         annotation_colors=color_list,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values for B/L ', 
                     '(# of tissues: ', format(nrow(deconv_dt), big.mark=','), ')'))
# overall L is at 13/19 high risk. cluster 4 Ls (cancer-like) are 5/6 high risk.
# looks significant but fisher's exact p=0.7208. zzz

rowname_subset <- rownames(cor_mat)[grepl('^L|^H', rownames(cor_mat))]
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         border_color=NA,
         clustering_method='ward.D2',
         cutree_cols=3, cutree_rows=3,
         annotation_col=annot_df,
         annotation_colors=color_list,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values for L/H ', 
                     '(# of tissues: ', format(nrow(deconv_dt), big.mark=','), ')'))
# cluster 4 H are mostly low-risk? make this make sense...

rowname_subset <- rownames(cor_mat)[grepl('^H|^C', rownames(cor_mat))]
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         border_color=NA,
         clustering_method='ward.D2',
         cutree_cols=3, cutree_rows=3,
         annotation_col=annot_df,
         annotation_colors=color_list,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values for H/C ', 
                     '(# of tissues: ', format(nrow(deconv_dt), big.mark=','), ')'))
# surprisingly cluster 3 (high bloods) signature clusters more strongly than
# cluster 4 (cancer-like, high erythrocyte progenitors) when H and C are
# grouped together. maybe the BL samples in cluster 3 add more noise; there
# are no B in cluster 4


# get a table of per-cluster, per-classification, per-risk group count
annot_df |> as_tibble() |> count(classification, cluster_ID) |> print(n=999)
annot_df |> as_tibble() |> count(classification, risk) |> print(n=999)
annot_df |> as_tibble() |> count(cluster_ID, risk) |> print(n=999)
annot_df |> as_tibble() |> count(classification, cluster_ID, risk) |> print(n=999)

# carry annot_df (for cluster IDs) to other scripts
saveRDS(annot_df, 'annot_df.rds')

# for replicability purposes
sessionInfo()
