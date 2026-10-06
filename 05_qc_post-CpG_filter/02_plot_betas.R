#!/usr/bin/env Rscript

"> 02_plot_betas.R <

Compares general trends for beta values from 96 datasets, via PCA and
hierarchical clustering of Pearson corr values.

Checks to see whether 'most variable positions' are a good way to sieve out
biologically interesting positions (TL;DR not really).

Also plots figures without N samples (i.e., B is the baseline) to greatly
reduce cell-type-specific methylation patterns (N are mostly squamous; B --> C
are columnar).
" -> doc

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(this.path)
})


setwd(this.path::here())
COMPILED_BETA_FILE <- '../04_filter_cpgs/all.beta.filt.tsv.gz'
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'
MEAN_COV_FILE <- './mean_cov_dt.tsv'
TOP_VAR_POS <- 100000

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
                     size=NULL, alpha=0.3, color_palette=NULL, 
                     title=NULL, subtitle=FALSE, hide_legend=FALSE) {
  # extract important stuff from the prcomp object
  pca_coords <- as.data.frame(prcomp_obj$x)
  eigs <- prcomp_obj$sdev ^ 2
  
  x_pct <- round(eigs[as.numeric(substr(x, 3, nchar(x)))] / sum(eigs) * 100, 2)
  y_pct <- round(eigs[as.numeric(substr(y, 3, nchar(y)))] / sum(eigs) * 100, 2)
  disease_class <- unname(short_to_long_disease_vec[substr(rownames(pca_coords), 1, 1)])
  disease_class <- factor(disease_class, levels=disease_progression_vec)
  g <- 
    ggplot(pca_coords, aes(x=!!sym(x), y=!!sym(y))) +
    geom_point(aes(color=disease_class), alpha=alpha, size=if (is.null(size)) 4 else size) +
    geom_text(aes(label=rownames(pca_coords)), check_overlap=TRUE, hjust='inward', vjust='inward') +
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
    g <- g + labs(subtitle=paste0('# of CpGs: ', format(nrow(prcomp_obj$rotation), big.mark=',')))
  }
  
  if (hide_legend) {
    g <- g + theme(legend.position='none')
  }
  
  g
}

# read compiled cov data, omit first three columns (chr/start/end)
# then sort by column name
beta_dt <- fread(COMPILED_BETA_FILE, sep='\t', header=TRUE, drop=1:3)
setcolorder(beta_dt, c(order(names(beta_dt))))
beta_dt[11:15, 1:48]  # just to see presence of NA

# how many CpGs are there...?
nrow(beta_dt)

# PCAs and heatmaps tend to be allergic to NAs, i.e., need to either drop
# rows containing NA, or do random forest to impute NAs
#
# how many rows have NAs?
table(rowSums(is.na(beta_dt)))

# hmm, doesn't look too bad. discard NAs outright, use remaining rows
beta_dt <- beta_dt[rowSums(is.na(beta_dt)) == 0]

# code has been crashing due to lack of memory, need to slim down `beta_dt`
# slightly (PCA of 13.6 mil x 96 = too large of a matrix)
# slim down by selecting every other row (split in twain)
beta_dt <- beta_dt[rep(c(TRUE, FALSE), length=.N), ]

#+ fig.width=10, fig.height=10
# plots to check similarity of beta values across datasets
# note: remember that PCAs/heatmaps abhor NAs
prcomp_obj <- run_prcomp(beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title='PCA of methylation levels from all 96 samples', subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

cor_mat <- cor(beta_dt, method='pearson')
diag(cor_mat) <- NA
pheatmap(cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

# layer coverage info onto plot to visually check influence of coverage
mean_cov_dt <- fread(MEAN_COV_FILE, sep='\t', header=TRUE)
plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, color_palette=disease_to_color_vec,
         title='PCA of methylation levels from all 96 samples; point sizes ∝ coverage', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info onto plot, for similar check
# mean_betas ranges (0.48, 0.73); accentuate differences by beta x 50 - 25,
# and fixing values < 3 to 3
mean_beta <- colMeans(beta_dt)
mean_beta
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title='PCA of methylation levels from all 96 samples; point sizes ∝ mean beta', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# C03/C07/C10/C11 look unlike all other samples, likely due to low coverages
# (refer to `01_plot_coverages.R`). exclude these samples, so that PCA patterns
# are not driven too much by technical variability
prcomp_obj <- run_prcomp(beta_dt[, .SD, .SDcols=!c('C03', 'C07', 'C10', 'C11')])
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title='PCA of methylation levels from 92 samples (excl. C03/C07/C10/C11)', subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# no need to recompute pairwise correlation matrix--just exclude cols/rows
# for C03/C07/C10/C11
pheatmap(cor_mat[setdiff(rownames(cor_mat), c('C03', 'C07', 'C10', 'C11')), setdiff(colnames(cor_mat), c('C03', 'C07', 'C10', 'C11'))],
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

# what if we excluded normals as well, how does IM -> cancer samples look on
# a PCA? for these plots, maintain exclusion of low coverage samples
sliced_beta_dt <- 
  beta_dt[, .SD, .SDcols=!grepl('^N', colnames(beta_dt))][, .SD, .SDcols=!c('C03', 'C07', 'C10', 'C11')]
prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title='PCA of methylation levels from 82 samples (excl. N & C03/C07/C10/C11)', subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer coverage info onto non-N and non-C03/C07/C10/C11
mean_cov_dt <- mean_cov_dt[mean_cov_dt$short_id %in% colnames(sliced_beta_dt), ]
plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, color_palette=disease_to_color_vec,
         title='PCA of methylation levels from 82 samples (excl. N & C03/C07/C10/C11); point sizes ∝ coverage', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info onto non-N and non-C03/C07/C10/C11
# mean_betas ranges (0.48, 0.73); accentuate differences by beta x 50 - 25,
# and fixing values < 3 to 3
mean_beta <- colMeans(sliced_beta_dt)
mean_beta
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title='PCA of methylation levels from 82 samples (excl. N & C03/C07/C10/C11); point sizes ∝ mean beta', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# and the heatmap for non-N and non-C03/C07/C10/C11
rowname_subset <- rownames(cor_mat)[!grepl('^N', rownames(cor_mat))]
rowname_subset <- setdiff(rowname_subset, c('C03', 'C07', 'C10', 'C11'))
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))


# would analysis be more meaningful if we constrained it to top 100k most
# variable positions? i.e., var() across samples is amongst the highest 100k
#
# start by excluding low-coverage samples first
cpg_var <- beta_dt[, apply(.SD, 1, var), .SDcols=!c('C03', 'C07', 'C10', 'C11')]
sliced_beta_dt <- beta_dt[order(cpg_var, decreasing=TRUE)[1:TOP_VAR_POS], !c('C03', 'C07', 'C10', 'C11')]

prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title='PCA of methylation levels from 92 samples (excl. C03/C07/C10/C11)', subtitle=TRUE)

mean_beta <- colMeans(sliced_beta_dt)   # mean beta of the top 100k positions
mean_beta
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from', ncol(sliced_beta_dt), 'samples; point sizes ∝ mean beta'), subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# then also exclude normals
cpg_var <- beta_dt[, .SD, .SDcols=!grepl('^N', colnames(beta_dt))][, .SD, .SDcols=!c('C03', 'C07', 'C10', 'C11')][, apply(.SD, 1, var)]
sliced_beta_dt <- beta_dt[, .SD, .SDcols=!grepl('^N', colnames(beta_dt))][, .SD, .SDcols=!c('C03', 'C07', 'C10', 'C11')][order(cpg_var, decreasing=TRUE)[1:TOP_VAR_POS], ]

prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from non-low coverage, non-N', ncol(sliced_beta_dt), 'samples'), subtitle=TRUE)

mean_beta <- colMeans(sliced_beta_dt)   # mean beta of the top 100k positions
mean_beta
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from', ncol(sliced_beta_dt), 'samples; point sizes ∝ mean beta'), subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# correlation matrix based on top 100k positions--different from overall view?
cor_mat <- cor(sliced_beta_dt, method='pearson')
diag(cor_mat) <- NA
pheatmap(cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(TOP_VAR_POS, big.mark=','), ')'))
   

# ok change tack--look at distribution of beta values in all datasets
melted_beta_dt <- melt(beta_dt, measure.vars=colnames(beta_dt),
                       variable.name='short_id', value.name='beta')
melted_beta_dt$disease_class <- factor(unname(short_to_long_disease_vec[substr(melted_beta_dt$short_id, 1, 1)]), levels=disease_progression_vec)
head(melted_beta_dt)

# overall view
ggplot(melted_beta_dt, aes(x=beta, color=disease_class, group=short_id)) + 
  # adjust controls kernel density, 1 is default
  geom_line(stat='density', adjust=0.5, position='identity', alpha=0.2) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(title='Distribution of beta, one line per sample, grouped by disease classification',
       x=expression(beta),
       y='Density',
       color='Classification') +
  facet_wrap(~ disease_class, ncol=1) +
  theme_minimal(12) +
  theme(legend.position='none')

#+ fig.width=10, fig.height=6
# superimpose NSq / BO / OAC to check whether cancer drives
# extensive hypomethylation?
# also kick out samples known to have low coverage (C03/C07/C10/C11)
ggplot(melted_beta_dt[!short_id %in% c('C03', 'C07', 'C10', 'C11')]
                     [disease_class %in% c('NSq', 'NDBE', 'EAC')],
       aes(x=beta, color=disease_class, group=short_id)) + 
  # adjust controls kernel density, 1 is default
  geom_line(aes(linetype=disease_class), stat='density', adjust=0.5,
            position='identity', alpha=0.4) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(title='Distribution of beta, coloured by disease classification',
       x=expression(beta),
       y='Density',
       color='Classification',
       linetype='Classification') +
  theme_minimal(12) +
  theme(legend.position='top')

# look at cumulative density, to judge cancer hypo/hypermethylation relative to norm
ggplot(melted_beta_dt[!short_id %in% c('C03', 'C07', 'C10', 'C11')]
                     [disease_class %in% c('NSq', 'NDBE', 'EAC')],
       aes(x=beta, color=disease_class, group=short_id)) + 
  # adjust controls kernel density, 1 is default
  geom_line(aes(linetype=disease_class), stat='ecdf', pad=FALSE,
            position='identity', alpha=0.4) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(title='Cumulative distribution of beta, coloured by disease classification',
       x=expression(beta),
       y='Cumulative density',
       color='Classification',
       linetype='Classification') +
  theme_minimal(12) +
  theme(legend.position='top')
# some cancers have global demeth, some cancers have global remeth?!

# for replicability purposes
sessionInfo()
