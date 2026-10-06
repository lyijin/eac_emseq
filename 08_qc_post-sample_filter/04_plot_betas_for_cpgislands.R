#!/usr/bin/env Rscript

"> 04_plot_betas_for_cpgislands.R <

Previously, trends were plotted for all 25.3m CpGs in the dataset. What would
trends look like if we focused on the 1.79m CpG islands / top-0.5% most
variable CpG islands instead?
" -> doc

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(pheatmap)
  library(this.path)
})


setwd(this.path::here())
COMPILED_BETA_FILE <- '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'
MEAN_COV_FILE <- './mean_cov_dt.tsv'
CPGISLAND_GR_FILE <- './island_beta_gr.rds'
TOP_VAR_POS <- 10000

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


# read in locations of CpGs in CpG islands
cpgisland_gr <- readRDS(CPGISLAND_GR_FILE)
cpgisland_dt <- data.table(chr=as.character(seqnames(cpgisland_gr)),
                           start=as.numeric(start(cpgisland_gr)))

# read compiled cov data, and include first three columns (chr/start/end)
beta_dt <- fread(COMPILED_BETA_FILE, sep='\t', header=TRUE)
# subset the 1.79m CpG islands with a left-join, then sort by colname.
# keep chr/start/end to facilitate top-X-most-variableslicing later on
beta_dt <- merge(cpgisland_dt, beta_dt, by=c('chr', 'start'), all.x=TRUE)
setcolorder(beta_dt, c(order(names(beta_dt))))
setcolorder(beta_dt, c('chr', 'start', 'end'))

# how many CpGs are there...?
nrow(beta_dt)

# confirm samples have been excluded. should be 86 here
ncol(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')])

# subset clinical data for these included samples
clin_dt <- fread(CLINDATA_FILE, sep='\t', header=TRUE)
clin_dt <- clin_dt[order(short_id)]
clin_dt[classification == 'Normal squamous', classification := 'NSq']
clin_dt[classification == 'IM', classification := 'NDBE']
clin_dt[classification == 'Cancer', classification := 'EAC']
clin_dt <- clin_dt[clin_dt$short_id %in% colnames(beta_dt), ]
stopifnot(clin_dt$short_id == colnames(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')]))  # sanity check to confirm equality

# distinguish high-risk --> adjacent / progression / NA (for cancers)
clin_dt$risk_type <- paste(clin_dt$risk, clin_dt$H_risk_category, sep='_')

# PCAs and heatmaps tend to be allergic to NAs, i.e., need to either drop
# rows containing NA, or do random forest to impute NAs
#
# how many rows have NAs?
table(rowSums(is.na(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')])))

# hmm, doesn't look too bad. discard NAs outright, use remaining rows
beta_dt <- beta_dt[rowSums(is.na(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')])) == 0]

# plot per-sample mean beta, categorised by classification
# calculate per-sample mean cov, add sample classification and short ID
mean_beta_dt <- data.table(classification=clin_dt$classification,
                           risk=clin_dt$risk,
                           short_id=clin_dt$short_id,
                           gender=clin_dt$gender,
                           age=clin_dt$age_at_collection,
                           mean_beta=colMeans(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')]))
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
#ggsave('raw_fig1c2.pdf', width=4, height=2.5)

#+ fig.width=10, fig.height=10
# plots to check similarity of beta values across datasets
# note: remember that PCAs/heatmaps abhor NAs
prcomp_obj <- run_prcomp(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')])
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from non-misclassified', ncol(beta_dt)-3, 'samples'), subtitle=TRUE)

#+ fig.width=5, fig.height=5
# cleaner plot for manuscript
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, size=3,
         hide_labels=TRUE, title='', subtitle=TRUE) +
  theme(legend.position='inside',
        legend.position.inside=c(1, 0),
        legend.justification=c(1, 0),
        legend.box.background=element_rect(fill='#ffffff80', color='#33333366'))
#ggsave('raw_fig1d.pdf', width=5, height=5)

#+ fig.width=10, fig.height=10
# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer in risk information to see whether high-risk vs. low-risk separates on 
# any of the principal components
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type,
         title=paste('PCA of methylation levels from non-misclassified', ncol(beta_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are risk-related trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)


# plot correlation matrix
cor_mat <- cor(beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')], method='pearson')
diag(cor_mat) <- NA
pheatmap(cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=7,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))

# what if we excluded normals as well, how does IM -> cancer samples look on
# a PCA?
sliced_beta_dt <- beta_dt[, .SD, .SDcols=!patterns('^N')]
sliced_beta_dt <- sliced_beta_dt[, .SD, .SDcols=!c('chr', 'start', 'end')]
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

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, risk_vector=sliced_clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# and the heatmap for non-N
rowname_subset <- rownames(cor_mat)[!grepl('^N', rownames(cor_mat))]
pheatmap(cor_mat[rowname_subset, rowname_subset],  # rows and cols are symmetrical
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=7,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(nrow(beta_dt), big.mark=','), ')'))


# would analysis be more meaningful if we constrained it to top 10k most
# variable positions? use variance that excludes Ns, as we are more interested
# in understanding disease progression than Barrett's onset (and
# disproportionately skews overall trends towards CpGs with large variance
# from N --> B)
topvar_gr <- cpgisland_gr[order(-cpgisland_gr$variance_noN)][1:TOP_VAR_POS]
topvar_dt <- data.table(chr=as.character(seqnames(topvar_gr)),
                        start=as.numeric(start(topvar_gr)))
topvar_dt <- merge(topvar_dt, beta_dt, by=c('chr', 'start'), all.x=TRUE)
# no need to retain chr/start/end
topvar_dt <- topvar_dt[, .SD, .SDcols=!c('chr', 'start', 'end')]
# drop any rows containing NAs
na_bool <- rowSums(is.na(topvar_dt)) > 0
topvar_dt <- topvar_dt[!na_bool]
topvar_gr <- topvar_gr[!na_bool]

#+ fig.width=10, fig.height=10
# plots to check similarity of beta values across datasets
# note: remember that PCAs/heatmaps abhor NAs
prcomp_obj <- run_prcomp(topvar_dt)
plot_pca(prcomp_obj, color_palette=disease_to_color_vec,
         title=paste('PCA of methylation levels from non-misclassified', ncol(topvar_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)

# layer in risk information to see whether high-risk vs. low-risk separates on 
# any of the principal components
plot_pca(prcomp_obj, color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type,
         title=paste('PCA of methylation levels from non-misclassified', ncol(beta_dt), 'samples'), subtitle=TRUE)

# also check other PCs to see whether there are risk-related trends there?
g1 <- plot_pca(prcomp_obj, x='PC2', y='PC3', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g2 <- plot_pca(prcomp_obj, x='PC3', y='PC4', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g3 <- plot_pca(prcomp_obj, x='PC4', y='PC5', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
g4 <- plot_pca(prcomp_obj, x='PC5', y='PC6', color_palette=disease_to_color_vec, risk_vector=clin_dt$risk_type, hide_legend=TRUE)
plot_grid(g1, g2, g3, g4, nrow=2)
# risk does not seem to segregate cleanly on any PC


# correlation matrix based on top 10k positions--different from overall view?
topvar_cor_mat <- cor(topvar_dt, method='pearson')
diag(topvar_cor_mat) <- NA
pheatmap(topvar_cor_mat,
         color=colorRampPalette(RColorBrewer::brewer.pal(n=7, name='YlGnBu'))(200)[1:150],
         clustering_method='ward.D2',
         show_colnames=TRUE,
         show_rownames=TRUE,
         fontsize=7,
         main=paste0('Clustered heatmap of pairwise Pearson r values (# of positions: ', 
                     format(TOP_VAR_POS, big.mark=','), ')'))

# perform hierarchical clustering on beta values, instead of cor values
# prepare col annots and row annots
col_annot_df <- readRDS('../14_deconv_celfieish/annot_df.rds')
row_annot_df <- data.frame(promoter=as.character(as.logical(topvar_gr$promoter)),
                           genebody=as.character(as.logical(topvar_gr$genebody)))
color_list <- list(
  classification=disease_to_color_vec[6:10],
  gender=c('M'='#67a9cf', 'F'='#ef8a62'),  # RdBu color scheme
  risk=c('H'='#e9a3c9', 'L'='#a1d76a'),    # PiYG color scheme
  cluster_celltype=c('Columnar'='#66c2a5', 'Squamous'='#8da0cb',
                     'Lymphoid'='#e5c494', 'Erythroid'='#fc8d62'), # Set2 color scheme
  promoter=c('TRUE'='#d8b365', 'FALSE'='#f5f5f5'),
  genebody=c('TRUE'='#5ab4ac', 'FALSE'='#f5f5f5'))  # BrBG color


# do the hierarchical clustering
topvar_mat <- as.matrix(topvar_dt)
rownames(topvar_mat) <- rownames(row_annot_df)
p <- pheatmap(topvar_mat,
              border_color=NA,
              clustering_method='ward.D2',
              cutree_cols=4, cutree_rows=2,
              treeheight_row=10,
              annotation_col=col_annot_df,
              annotation_row=row_annot_df,
              annotation_colors=color_list,
              show_rownames=FALSE,
              fontsize_col=8)
#ggsave('raw_fig3zz.pdf', plot=p$gt, width=12, height=8)


# for replicability purposes
sessionInfo()
