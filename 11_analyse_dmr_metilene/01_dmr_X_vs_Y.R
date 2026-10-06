#!/usr/bin/env Rscript

"> 01_dmr_X_vs_Y.R <

This script handles pairwise comparisons that do not need additional analyses
e.g., compare against previously published work. Two other scripts in this
same folder ('02_dmr_HC_vs_N.R'/'03_dmr_HC_vs_B.R') handles those special cases.

Nomenclature alert: in terms of disease progression, X is ALWAYS the more
serious state, Y is the less serious one.
" -> doc

suppressPackageStartupMessages({
  library(eulerr)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(this.path)
})


setwd(this.path::here())
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'
GENOMIC_ANNOT_FILE <- '../data/gencode.v43.annotation_2026-02-19.RData'

# error out if any of the files do not exist
stopifnot(file.exists(CLINDATA_FILE))
stopifnot(file.exists(GENOMIC_ANNOT_FILE))

# load pre-processed genomic annotations and functions to annotate GRanges
load(GENOMIC_ANNOT_FILE)
source('../01_txdb/annotate_ranges_functions.R')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# check the directory for all *.tsv.gz files, except HC_vs_N and HC_vs_B
METILENE_FILES <- Sys.glob('../10_call_dmr_metilene/metilene_output/*.tsv.gz')
METILENE_FILES <- METILENE_FILES[!grepl('allHC_vs_all[B|N]', METILENE_FILES, perl=TRUE)]

#+ fig.width=10, fig.height=3
for (METILENE_FILE in METILENE_FILES) {
  diag_message('Processing ', METILENE_FILE, '...')
  
  dmr_dt <- fread(METILENE_FILE, sep='\t', header=FALSE,
                  col.names=c('chr', 'start', 'stop', 'qval', 'delta_beta',
                              'nCpG', 'p_MWU', 'p_KS', 'mean_g1', 'mean_g2'))
  diag_message(METILENE_FILE, ' has ',
               formatC(nrow(dmr_dt), big.mark=','), ' rows.')
  
  # retain rows that are "great" DMRs
  #   1. metilene has a min delta beta of 0.1 -- do nothing here
  #   2. qval < 0.01 (default is 0.05, but this is a fairly light-touch filter)
  dmr_dt <- dmr_dt[qval < 0.01]
  diag_message('Post qval filtering, ', formatC(nrow(dmr_dt), big.mark=','), ' rows remain.')
  
  # quick histogram to check extent of delta_beta
  g <- ggplot(mapping=aes(dmr_dt$delta_beta)) +
    geom_histogram(binwidth=0.01) +
    ggtitle(METILENE_FILE) +
    theme_minimal(14)
  print (g)
  
  # create GRanges from this data.table, and retain other cols as metadata
  # start/stop coords from metilene are already 1-based, confirmed in JBrowse
  dmr_gr <- makeGRangesFromDataFrame(
    dmr_dt[, c('chr', 'start', 'stop', 'qval', 'delta_beta', 'nCpG', 'mean_g1', 'mean_g2')],
    keep.extra.columns=TRUE,
    ignore.strand=TRUE)
  
  # annotate these regions
  dmr_gr <- annotateDMRs(dmrs=dmr_gr, tx=gencode_tx, tss_proximal=2000)
  
  # write all hypermeth regions into external file for eyeballing
  HYPER_OUTPUT_FILE <- sub('dmr.tsv.gz', 'hypermeth.tsv', METILENE_FILE)
  HYPER_OUTPUT_FILE <- sub('10_call_dmr_metilene/metilene_output/', 
                           '11_analyse_dmr_metilene/filtered_dmrs/', HYPER_OUTPUT_FILE)
  tmp_gr <- dmr_gr[order(dmr_gr$delta_beta, decreasing=TRUE), ]
  tmp_gr <- tmp_gr[tmp_gr$delta_beta > 0, ]
  write.table(tmp_gr, file=HYPER_OUTPUT_FILE, quote=FALSE, sep='\t', row.names=FALSE)
  
  # repeat for hypometh regions
  HYPO_OUTPUT_FILE <- sub('dmr.tsv.gz', 'hypometh.tsv', METILENE_FILE)
  HYPO_OUTPUT_FILE <- sub('10_call_dmr_metilene/metilene_output/',
                          '11_analyse_dmr_metilene/filtered_dmrs/', HYPO_OUTPUT_FILE)
  tmp_gr <- dmr_gr[order(dmr_gr$delta_beta, decreasing=FALSE), ]
  tmp_gr <- tmp_gr[tmp_gr$delta_beta < 0, ]
  write.table(tmp_gr, file=HYPO_OUTPUT_FILE, quote=FALSE, sep='\t', row.names=FALSE)
}


# for replicability purposes
sessionInfo()
