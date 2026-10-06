#!/usr/bin/env Rscript

"> 02_plot_betas_with_excl.R <

Compares general trends for beta values from 86 datasets, via PCA and
hierarchical clustering of Pearson corr values.

Also plots figures without N samples (i.e., B is the baseline) to greatly
reduce cell-type-specific methylation patterns (N are mostly squamous; B --> C
are columnar).

After considering exploratory RNA-seq and EM-seq results, some samples look
obviously misclassified. Decision was made to exclude these samples. To reduce
duplication in code, this exclusion was carried out with a Python script
(in `../00_common/`) to produce the `all.beta.filt.samp_excl.tsv.gz` file.
An R script in the same directory handles sample exclusion in R dataframes.

Coverage issues i.e., in C03/C07/C10/C11 are NOT considered in this script.
Sample exclusion here is purely due to misclassification, not technical issues
arising from sequencing.
" -> doc

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(this.path)
})


setwd(this.path::here())
COMPILED_BETA_FILE <- '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
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
beta_dt[90:95, 30:39]  # just to confirm presence of NA

# how many CpGs are there...?
nrow(beta_dt)

# confirm samples have been excluded. n should be < 96 here
ncol(beta_dt)

# subset clinical data for these included samples
clin_dt <- fread(CLINDATA_FILE, sep='\t', header=TRUE)
clin_dt <- clin_dt[order(short_id)]
clin_dt[classification == 'Normal squamous', classification := 'NSq']
clin_dt[classification == 'IM', classification := 'NDBE']
clin_dt[classification == 'Cancer', classification := 'EAC']
clin_dt <- clin_dt[clin_dt$short_id %in% colnames(beta_dt), ]
stopifnot(clin_dt$short_id == colnames(beta_dt))  # sanity check to confirm equality

# distinguish high-risk --> adjacent / progression / NA (for cancers)
clin_dt$risk_type <- paste(clin_dt$risk, clin_dt$H_risk_category, sep='_')

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

# plot per-sample mean beta, categorised by classification
# calculate per-sample mean cov, add sample classification and short ID
mean_beta_dt <- data.table(classification=clin_dt$classification,
                           risk=clin_dt$risk,
                           short_id=clin_dt$short_id,
                           gender=clin_dt$gender,
                           age=clin_dt$age_at_collection,
                           mean_beta=colMeans(beta_dt))
mean_beta_dt$classification <- factor(mean_beta_dt$classification, 
                                      levels=disease_progression_vec)
mean_beta_dt$age <- as.numeric(mean_beta_dt$age)
#+ fig.width=4, fig.height=2.5
ggplot(mean_beta_dt, aes(x=classification, y=mean_beta)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=classification), position=position_jitter(seed=1337), size=3, alpha=0.6) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(x='Disease classification',
       y='Mean beta') +
  theme_minimal(12) +
  theme(legend.position='none')
ggsave('raw_fig1c2.pdf', width=4, height=2.5)

#+ fig.width=10, fig.height=10
# plots to check similarity of beta values across datasets
# note: remember that PCAs/heatmaps abhor NAs
prcomp_obj <- run_prcomp(beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from non-misclassified', ncol(beta_dt), 'samples'), subtitle=TRUE)

#+ fig.width=5, fig.height=5
# cleaner plot for manuscript
pca_coords <- as.data.frame(prcomp_obj$x)
disease_class <- unname(short_to_long_disease_vec[substr(rownames(pca_coords), 1, 1)])
disease_class <- factor(disease_class, levels=disease_progression_vec)
f1d <-
  plot_pca(prcomp_obj, color_palette=disease_to_color_vec, size=3,
           hide_labels=TRUE, title='', subtitle=TRUE) +
  theme_classic(12) +
  theme(legend.position='inside',
        legend.position.inside=c(0.98, 0.02),
        legend.justification=c(1, 0),
        legend.box.background=element_rect(fill='#ffffff80', color='#33333366'))
xdens <-
  axis_canvas(f1d, axis='x') +
  geom_density(data=as.data.frame(prcomp_obj$x),
               aes(x=PC1, fill=disease_class, color=disease_class), alpha=0.3) +
  scale_color_manual(values=disease_to_color_vec) +
  scale_fill_manual(values=disease_to_color_vec)
ydens <-
  axis_canvas(f1d, axis='y', coord_flip=TRUE) +
  geom_density(data=as.data.frame(prcomp_obj$x),
               aes(x=PC2, fill=disease_class, color=disease_class), alpha=0.3) +
  scale_color_manual(values=disease_to_color_vec) +
  scale_fill_manual(values=disease_to_color_vec) +
  coord_flip()
f1d %>%
  insert_xaxis_grob(xdens, grid::unit(0.1, "null"), position='top') %>%
  insert_yaxis_grob(ydens, grid::unit(0.1, "null"), position='right') %>%
  ggdraw()
ggsave('raw_fig1d.pdf', width=5, height=5)

#+ fig.width=10, fig.height=10
# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer in risk information to see whether high-risk vs. low-risk separates on 
# any of the principal components
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type,
         title=paste('PCA of methylation levels from non-misclassified', ncol(beta_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are risk-related trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer coverage info onto plot to visually check influence of coverage
mean_cov_dt <- fread(MEAN_COV_FILE, sep='\t', header=TRUE)
mean_cov_dt <- mean_cov_dt[mean_cov_dt$short_id %in% colnames(beta_dt), ]
stopifnot(mean_cov_dt$short_id == colnames(beta_dt))  # sanity check to confirm equality
plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from', ncol(beta_dt), 'samples; point sizes ∝ coverage'), subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=mean_cov_dt$mean_cov/2, x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info onto plot, for similar check
# mean_betas ranges (0.48, 0.73); accentuate differences by beta x 50 - 25,
# and fixing values < 3 to 3
mean_beta <- mean_beta_dt$mean_beta
mean_beta
plot_pca(prcomp_obj, size=pmax(mean_beta_dt$mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from', ncol(beta_dt), 'samples; point sizes ∝ mean beta'), subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# plot correlation matrix
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

# what if we excluded normals as well, how does IM -> cancer samples look on
# a PCA?
sliced_beta_dt <- beta_dt[, .SD, .SDcols=!patterns('^N')]
sliced_clin_dt <- clin_dt[short_id %in% colnames(sliced_beta_dt), ]
prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from non-misclassified, non-N', ncol(sliced_beta_dt), 'samples'), subtitle=TRUE) +
  theme(legend.position='inside',
        legend.position.inside=c(1, 0),
        legend.justification=c(1, 0),
        legend.box.background=element_rect(fill='#ffffff80', color='#33333366'))
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type,
         title=paste('PCA of methylation levels from non-misclassified, non-N', ncol(sliced_beta_dt), 'samples'), subtitle=TRUE)
# cleaner plot for manuscript, and add density plots along x/y axes
pca_coords <- as.data.frame(prcomp_obj$x)
disease_class <- unname(short_to_long_disease_vec[substr(rownames(pca_coords), 1, 1)])
disease_class <- factor(disease_class, levels=disease_progression_vec)
sf1 <-
  plot_pca(prcomp_obj, color_palette=disease_to_color_vec, size=3,
           hide_labels=TRUE, title='', subtitle=TRUE) +
  theme_classic(12) +
  theme(legend.position='inside',
        legend.position.inside=c(0.97, 0.97),
        legend.justification=c(1, 1),
        legend.box.background=element_rect(fill='#ffffff80', color='#33333366'))
xdens <-
  axis_canvas(sf1, axis='x') +
  geom_density(data=as.data.frame(prcomp_obj$x),
               aes(x=PC1, fill=disease_class, color=disease_class), alpha=0.3) +
  scale_color_manual(values=disease_to_color_vec) +
  scale_fill_manual(values=disease_to_color_vec)
ydens <-
  axis_canvas(sf1, axis='y', coord_flip=TRUE) +
  geom_density(data=as.data.frame(prcomp_obj$x),
               aes(x=PC2, fill=disease_class, color=disease_class), alpha=0.3) +
  scale_color_manual(values=disease_to_color_vec) +
  scale_fill_manual(values=disease_to_color_vec) +
  coord_flip()
sf1 %>%
  insert_xaxis_grob(xdens, grid::unit(0.1, "null"), position='top') %>%
  insert_yaxis_grob(ydens, grid::unit(0.1, "null"), position='right') %>%
  ggdraw()
ggsave('raw_suppfig1.pdf', width=5, height=5)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info
mean_beta <- colMeans(sliced_beta_dt)
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title='PCA of methylation levels from non-misclassified, non-N; point sizes ∝ mean beta', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# and the heatmap for non-N
rowname_subset <- rownames(cor_mat)[!grepl('^N', rownames(cor_mat))]
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
# for non-N?
cpg_var <- beta_dt[, apply(.SD, 1, var), .SDcols=!patterns('^N')]
sliced_beta_dt <- beta_dt[order(cpg_var, decreasing=TRUE)[1:TOP_VAR_POS], .SD, .SDcols=!patterns('^N')]
sliced_clin_dt <- clin_dt[short_id %in% colnames(sliced_beta_dt), ]

prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt[short_id %in% colnames(sliced_beta_dt), ]$risk_type,
         title=paste('PCA of methylation levels from non-misclassified, non-N', ncol(sliced_beta_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info
mean_beta <- colMeans(sliced_beta_dt)
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title='PCA of methylation levels from non-misclassified, non-N; point sizes ∝ mean beta', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# correlation matrix based on top 100k positions--different from overall view?
most_var_cor_mat <- cor(sliced_beta_dt, method='pearson')
diag(most_var_cor_mat) <- NA
pheatmap(most_var_cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(TOP_VAR_POS, big.mark=','), ')'))

# remove C (as it only has high-risk), to see N -> B -> L -> H progression
# with risk potentially explaining a portion of the progression
cpg_var <- beta_dt[, apply(.SD, 1, var), .SDcols=!patterns('^C')]
sliced_beta_dt <- beta_dt[order(cpg_var, decreasing=TRUE)[1:TOP_VAR_POS], .SD, .SDcols=!patterns('^C')]
sliced_clin_dt <- clin_dt[short_id %in% colnames(sliced_beta_dt), ]

prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt[short_id %in% colnames(sliced_beta_dt), ]$risk_type,
         title=paste('PCA of methylation levels from non-misclassified, non-C', ncol(sliced_beta_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info
mean_beta <- colMeans(sliced_beta_dt)
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from non-misclassified, non-C', ncol(sliced_beta_dt), 'samples; point sizes ∝ mean beta'), subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# correlation matrix based on top 100k positions--different from overall view?
most_var_cor_mat <- cor(sliced_beta_dt, method='pearson')
diag(most_var_cor_mat) <- NA
pheatmap(most_var_cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(TOP_VAR_POS, big.mark=','), ')'))

# remove C (as it only has high-risk) and N, to focus on B -> L -> H progression
# with risk potentially explaining a portion of the progression
cpg_var <- beta_dt[, apply(.SD, 1, var), .SDcols=!patterns('^N|C')]
sliced_beta_dt <- beta_dt[order(cpg_var, decreasing=TRUE)[1:TOP_VAR_POS], .SD, .SDcols=!patterns('^N|C')]
sliced_clin_dt <- clin_dt[short_id %in% colnames(sliced_beta_dt), ]

prcomp_obj <- run_prcomp(sliced_beta_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type,
         title=paste('PCA of methylation levels from non-misclassified, non-N/C', ncol(sliced_beta_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC1', y='PC3', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC1', y='PC4', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC1', y='PC5', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC1', y='PC6', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer mean beta info
mean_beta <- colMeans(sliced_beta_dt)
plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), color_palette=disease_to_color_vec,
         title='PCA of methylation levels from non-misclassified, non-N/C; point sizes ∝ mean beta', subtitle=TRUE)
g1 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, size=pmax(mean_beta*50-25, 3), x='PC1', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# correlation matrix based on top 100k positions--different from overall view?
most_var_cor_mat <- cor(sliced_beta_dt, method='pearson')
diag(most_var_cor_mat) <- NA
pheatmap(most_var_cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(TOP_VAR_POS, big.mark=','), ')'))

# plot correlation matrices for each status separately. risk labels are included
# for every sample, to hopefully see low-risk stratifying from high-risk, etc
#
# make sure row/col order of cor_mat is identical to clin_dt$short_id, as the
# risk labels will soon be concatenated into cor_mat
stopifnot(rownames(cor_mat) == clin_dt$short_id)
stopifnot(rownames(cor_mat) == mean_cov_dt$short_id)

mean_beta <- colMeans(beta_dt)
stopifnot(rownames(cor_mat) == names(mean_beta))

cor_mat_with_risk <- cor_mat
rownames(cor_mat_with_risk) <-
  paste(rownames(cor_mat_with_risk), clin_dt$risk_type, clin_dt$gender,
        gsub('one visit', '1V', clin_dt$worst_diagnosis_known),
        signif(mean_cov_dt$mean_cov, 3), signif(mean_beta, 2), sep='_')
colnames(cor_mat_with_risk) <-
  paste(colnames(cor_mat_with_risk), clin_dt$risk_type, clin_dt$gender,
        gsub('one visit', '1V', clin_dt$worst_diagnosis_known), 
        signif(mean_cov_dt$mean_cov, 3), signif(mean_beta, 2), sep='_')
cor_mat_with_risk[1:5, 1:5]

# plot per-category heatmaps
rowname_subset <- rownames(cor_mat_with_risk)[grepl('^N', rownames(cor_mat_with_risk))]
pheatmap(cor_mat_with_risk[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

rowname_subset <- rownames(cor_mat_with_risk)[grepl('^B', rownames(cor_mat_with_risk))]
pheatmap(cor_mat_with_risk[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

rowname_subset <- rownames(cor_mat_with_risk)[grepl('^L', rownames(cor_mat_with_risk))]
pheatmap(cor_mat_with_risk[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

rowname_subset <- rownames(cor_mat_with_risk)[grepl('^H', rownames(cor_mat_with_risk))]
pheatmap(cor_mat_with_risk[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=8,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

# ok change tack--look at distribution of beta values in all datasets
melted_beta_dt <- melt(beta_dt, measure.vars=colnames(beta_dt),
                       variable.name='short_id', value.name='beta')
# inter-join both tables together
melted_beta_dt <- melted_beta_dt[
  clin_dt[, .SD, .SDcols=c('short_id', 'classification', 'risk_type')],
  on=.(short_id),
  nomatch=NULL]
head(melted_beta_dt)

# overall view
ggplot(melted_beta_dt, aes(x=beta, color=classification, group=short_id)) + 
  # adjust controls kernel density, 1 is default
  geom_line(stat='density', adjust=0.5, position='identity', alpha=0.2) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(title='Distribution of beta, one line per sample, grouped by disease classification',
       x=expression(beta),
       y='Density',
       color='Classification') +
  facet_wrap(~ classification, ncol=1) +
  theme_minimal(12) +
  theme(legend.position='none')

#+ fig.width=10, fig.height=6
# superimpose IM / HGD / cancer to check whether cancer drives
# extensive hypomethylation?
ggplot(melted_beta_dt[classification %in% c('NDBE', 'HGD', 'EAC')],
       aes(x=beta, color=classification, group=short_id)) + 
  # adjust controls kernel density, 1 is default
  geom_line(aes(linetype=classification), stat='density', adjust=0.5,
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
ggplot(melted_beta_dt[classification %in% c('NDBE', 'HGD', 'EAC')],
       aes(x=beta, color=classification, group=short_id)) + 
  # adjust controls kernel density, 1 is default
  geom_line(aes(linetype=classification), stat='ecdf', pad=FALSE,
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
