#!/usr/bin/env Rscript

"> 03_plot_contexts.R <

What sort of genomic contexts houses the most variable CpGs?
" -> doc

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(cowplot)
  library(data.table)
  library(forcats)
  library(ggplot2)
  library(GenomicRanges)
  library(this.path)
})


setwd(this.path::here())
COMPILED_BETA_FILE <- '../04_filter_cpgs/all.beta.filt.tsv.gz'
GENOMIC_ANNOT_FILE <- '../data/gencode.v43.annotation_2026-02-19.RData'
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

# load pre-processed genomic annotations and functions to annotate GRanges
load(GENOMIC_ANNOT_FILE)
source('../01_txdb/annotate_ranges_functions.R')

# function to generate ggplot from GRanges object
plot_trends_by_annots <- function (beta_gr) {
  # subfunction to wide-to-long convert a sliced GenomicRanges, then slap a
  # `feature_label`
  genomic_feature_as_long_dt <- function(beta_gr, feature_label) {
    wide_dt <- as.data.table(
      mcols(beta_gr[, grepl('^mean', names(mcols(beta_gr)))]))
    long_dt <- melt(wide_dt, measure.vars=names(wide_dt),
                    variable.name='classification', value.name='beta')
    
    # convert "meanN" to semantically meaningful "Normal squamous", etc.
    disease_class <-
      short_to_long_disease_vec[substr(long_dt$classification, 5, 5)]
    disease_class <- factor(disease_class, levels=disease_progression_vec)
    long_dt$classification <- disease_class
    
    # slap `feature_label` into the dt
    long_dt$genomic_feature <- paste0(
      feature_label,
      '|(n = ', format(nrow(wide_dt), big.mark=','), ')')
    
    long_dt
  }
  
  
  # plot overall beta value distributions based on CpG island / genic annots
  plot_dt <- rbind(
    genomic_feature_as_long_dt(beta_gr[beta_gr$CpGisland == 1], 'CpG islands'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$CpGshores == 1], 'CpG shores'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$nonCpG == 1], 'Non-island/shore'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$promoter == 1], 'Promoters'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$genebody == 1], 'Gene bodies'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$intergenic == 1], 'Intergenic'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$CRE == 1], 'Cis-regulatory elements'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$repeats == 1], 'RepeatMasker repeats'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$simple_repeats == 1], 'Simple tandem repeats')
  )
  
  g <- ggplot(plot_dt, aes(x=beta, y=fct_inorder(genomic_feature))) +
    geom_boxplot(aes(color=classification, fill=classification), alpha=0.5, width=0.8, outlier.shape=NA) +
    scale_color_manual('Classification', values=disease_to_color_vec) +
    scale_fill_manual('Classification', values=disease_to_color_vec) +
    scale_x_continuous(expression(beta)) +
    scale_y_discrete('', labels=function(x) gsub('|', '\n', x, fixed=TRUE), limits=rev) +
    labs(title=paste0('Beta distributions by genomic feature (overall n = ', format(nrow(mcols(beta_gr)), big.mark=','), ')')) +
    theme_classic(12) +
    theme(legend.position='top', 
          axis.line.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=11))
  
  g
}


# read compiled cov data, and then
#   1. remove 'lambda' and 'pUC19'
#   2. sort first three columns (chr/start/end) to the left
#   3. sort by column name for other columns
beta_dt <- fread(COMPILED_BETA_FILE, sep='\t', header=TRUE)
beta_dt <- beta_dt[chr != 'lambda'][chr != 'pUC19']
setcolorder(beta_dt, names(beta_dt)[c(1:3, order(names(beta_dt)[4:ncol(beta_dt)]) + 3)])
beta_dt[1:5, 1:10]

# how many CpGs are there...?
nrow(beta_dt)

# calculate group means for downstream filtering (e.g., where are the
# hypermethylated positions in C preferentially located in the genome?)
beta_dt$variance <- apply(beta_dt[, 4:ncol(beta_dt)], 1, var, na.rm=TRUE)
beta_dt$meanN <- rowMeans(
  beta_dt[, .SD, .SDcols=grepl('^N', colnames(beta_dt))], na.rm=TRUE)
beta_dt$meanB <- rowMeans(
  beta_dt[, .SD, .SDcols=grepl('^B', colnames(beta_dt))], na.rm=TRUE)
beta_dt$meanL <- rowMeans(
  beta_dt[, .SD, .SDcols=grepl('^L', colnames(beta_dt))], na.rm=TRUE)
beta_dt$meanH <- rowMeans(
  beta_dt[, .SD, .SDcols=grepl('^H', colnames(beta_dt))], na.rm=TRUE)
beta_dt$meanC <- rowMeans(
  beta_dt[, .SD, .SDcols=grepl('^C', colnames(beta_dt))], na.rm=TRUE)

# create GenomicRanges
beta_gr <- makeGRangesFromDataFrame(
  beta_dt[, c('chr', 'start', 'end', 'variance',
              'meanN', 'meanB', 'meanL', 'meanH', 'meanC')],
  keep.extra.columns=TRUE,
  ignore.strand=TRUE)
beta_gr

# beta_dt not used from this point on--free up some mem
rm(beta_dt)
gc()

# annotate `beta_gr`
beta_gr <- annotateDMRs(dmrs=beta_gr, tx=gencode_tx, tss_proximal=2000)
head(beta_gr)

# generate plot
#+ fig.width=10, fig.height=8
plot_trends_by_annots(beta_gr)

# let's now look at the top N most variable positions
most_variable_cpgs <- order(beta_gr$variance, decreasing=TRUE)[1:TOP_VAR_POS]

# confirm variable is, indeed, storing the most variable CpGs
beta_gr[most_variable_cpgs[1:5]]
summary(beta_gr$variance)  # max value should be in the GRanges object slice

# do most variable positions tend to be associated with certain genomic loci?
plot_trends_by_annots(beta_gr[most_variable_cpgs])
# hmm, median betas tend to hover around the 0.5 mark. makes mathematical sense
# (betas that are near-0 or near-1 has only one direction for values to deviate
# away from; betas at 0.5 have two directions for deviations, leading to larger
# variances), but does not a good biological story make


# let's define variability by largest delta difference in beta values, between
# groups. N vs. C first
most_variable_cpgs <- order(abs(beta_gr$meanC - beta_gr$meanN), decreasing=TRUE)[1:TOP_VAR_POS]
plot_trends_by_annots(beta_gr[most_variable_cpgs])

# follow same definition of variability, but for B vs. C
most_variable_cpgs <- order(abs(beta_gr$meanC - beta_gr$meanB), decreasing=TRUE)[1:TOP_VAR_POS]
plot_trends_by_annots(beta_gr[most_variable_cpgs])
# visually confirms code works as expected, abs IS getting a mix of diff up &
# diff downs together. unlike N vs. C plot where C was greater than N in all
# categories
#
# n-value-wise, baseline ratio of CpG island:shores:neither is about 1:1.7:17
# the later plots have more "neither" than expected. 1:1.7:27 for B vs. C. hmm.
# similar increase in promoter:genebody:intergenic, baseline was 1:2.2:1.2,
# B vs. C has ratio of 1:2.9:1.7. hints at dysregulation of gene body meth

# let's split apart hyper and hypometh for B vs. C. hyper in C first
most_variable_cpgs <- order(beta_gr$meanC - beta_gr$meanB, decreasing=TRUE)[1:as.integer(TOP_VAR_POS/2)]
plot_trends_by_annots(beta_gr[most_variable_cpgs])

# let's split apart hyper and hypometh for B vs. C. hypo in C next
most_variable_cpgs <- order(beta_gr$meanC - beta_gr$meanB, decreasing=FALSE)[1:as.integer(TOP_VAR_POS/2)]
plot_trends_by_annots(beta_gr[most_variable_cpgs])
# interesting to see step-ladder like behaviour across disease classification
# severity--code does not enforce this behaviour--so it seems hyper/hypometh
# progression is gradual in nature (analogue-ish), not binary (digital-ish)


# for replicability purposes
sessionInfo()
