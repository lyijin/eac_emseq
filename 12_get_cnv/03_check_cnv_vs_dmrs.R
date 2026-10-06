#!/usr/bin/env Rscript

"> 03_check_cnv_vs_dmrs.R <

Do DMRs tend to sit in CNVs--or not really?
" -> doc

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(pheatmap)
  library(this.path)
  library(tidyr)
})


setwd(this.path::here())
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'

# error out if any of the files do not exist
stopifnot(file.exists(CLINDATA_FILE))

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


# check the DMR directory for the core six metilene comparisons
METILENE_FILES <-
  c('../11_analyse_dmr_metilene/filtered_dmrs/allB_vs_allN.hypermeth.tsv.gz',
    '../11_analyse_dmr_metilene/filtered_dmrs/allC_vs_allN.hypermeth.tsv.gz',
    '../11_analyse_dmr_metilene/filtered_dmrs/allL_vs_allB.hypermeth.tsv.gz',
    '../11_analyse_dmr_metilene/filtered_dmrs/allHC_vs_allB.hypermeth.tsv.gz',
    '../11_analyse_dmr_metilene/filtered_dmrs/allH_vs_allL.hypermeth.tsv.gz',
    '../11_analyse_dmr_metilene/filtered_dmrs/allC_vs_allH.hypermeth.tsv.gz')

# create a list to store all GRanges objects
dmr_gr <- list()
for (mf in METILENE_FILES) {
  # defensive programming
  stopifnot(file.exists(mf))
  
  # parse filename for contextual clues for disease states
  fname_split <- unlist(strsplit(mf, '/', fixed=TRUE))
  fname_split <- unlist(strsplit(fname_split[length(fname_split)], '.', fixed=TRUE))
  hyper_or_hypo <- fname_split[2]
  compared_states <- fname_split[1]
  compared_states <- unlist(strsplit(compared_states, '_vs_'))
  # for this analysis, drop the "all" label for stuff like "allN"
  compared_states <- sub('all', '', compared_states)
  compared_states <- paste0(short_to_long_disease_vec[compared_states[1]], ' vs. ',
                            short_to_long_disease_vec[compared_states[2]])
  diag_message('Reading ', mf, ' (', hyper_or_hypo, ' markers for ',
               compared_states, ').')
  
  dmr_gr[[compared_states]] <- makeGRangesFromDataFrame(
    fread(mf, select=c(1,2,3,7,8,9,10,21)),
    keep.extra.columns=TRUE,
    ignore.strand=TRUE)
  
  # CONSTRAIN ANALYSIS TO AUTOSOMES. HAVING DMRs ON CHROM X GAVE THE IMPRESSION
  # THAT THERE WERE MASSIVE CNV, BUT TURNS OUT BEING MALE WAS THE REASON WHY
  dmr_gr[[compared_states]] <- 
    dmr_gr[[compared_states]][!seqnames(dmr_gr[[compared_states]]) %in% c('chrX', 'chrY', 'chrM')]
  # there weren't any DMRs on Y/M, but to be safe, at the expense of a bit
  # more computational power
}
dmr_gr

# read per-sample inferred CNV values from `cnvkit.py call`
cnvcall_files <- c(Sys.glob('03_cnvkit_call_cnv/B*.call.cns'),
                   Sys.glob('03_cnvkit_call_cnv/L*.call.cns'),
                   Sys.glob('03_cnvkit_call_cnv/H*.call.cns'),
                   Sys.glob('03_cnvkit_call_cnv/C*.call.cns'))
cnvcall_gr <- list()
for (cf in cnvcall_files) {
  # defensive programming
  stopifnot(file.exists(cf))
  
  # get the three-char short_id from the filename
  cf_id <- gsub('.*/', '', gsub('\\..*$', '', cf))
  
  # read the CNV file, convert into GenomicRanges object for `findOverlaps()`
  cnvcall_gr[[cf_id]] <- makeGRangesFromDataFrame(
    fread(cf, select=c(1,2,3,8)),
    keep.extra.columns=TRUE,
    ignore.strand=TRUE)
}
head(cnvcall_gr, 3)  # to illustrate that every sample has different ranges

# carry out the per-comparison, per-sample tallying of # DMRs are in CNVs,
# and how many aren't
dmr_cn_tib <- NULL
all_dmr_cn_tib <- NULL
for (cf_id in names(cnvcall_gr)) {
  # check whether any of the six comparisons have DMRs overlapping CNVs
  for (compared_states in names(dmr_gr)) {
    fo_dmr_cnvcall <- findOverlaps(dmr_gr[[compared_states]], cnvcall_gr[[cf_id]], type='within')
    
    # create tibble to store CNVs for top 50 DMRs (this will be visualised with
    # pheatmap later)
    top_n <- 50
    
    # holy spaghetti code. this bit rbinds the previous incarnation of the
    # DMR copy number tibble ("dmr_cn_tib") with the next one...
    dmr_cn_tib <- bind_rows(
      dmr_cn_tib,
      tibble(
        # ... chrom/start was used instead of rank as rank order can differ
        # if a loci has no copy number info in some samples (findOverlaps
        # not being sensitive to no-overlaps)
        chrom=as.character(seqnames(dmr_gr[[compared_states]])[queryHits(fo_dmr_cnvcall)[1:top_n]]),
        start=start(dmr_gr[[compared_states]])[queryHits(fo_dmr_cnvcall)[1:top_n]],
        comparison=as.factor(compared_states),
        short_id=as.factor(cf_id),
        delta_beta=mcols(dmr_gr[[compared_states]])$delta_beta[queryHits(fo_dmr_cnvcall)[1:top_n]],
        gene=mcols(dmr_gr[[compared_states]])$gene_name_prot[queryHits(fo_dmr_cnvcall)[1:top_n]],
        cn=mcols(cnvcall_gr[[cf_id]])$cn[subjectHits(fo_dmr_cnvcall)[1:top_n]]) |>
        # this bit combines four columns together with unite() to create a
        # index column. this is for the pheatmap later, as it wants a single
        # index column and many value columns
        unite('chrom_start_comparison_gene', c(chrom, start, comparison, gene)))
    
    all_dmr_cn_tib <- bind_rows(
      all_dmr_cn_tib,
      tibble(
        # repeat for all DMRs to check whether observations of top 50 is
        # representative of the entire set
        chrom=as.character(seqnames(dmr_gr[[compared_states]])[queryHits(fo_dmr_cnvcall)]),
        start=start(dmr_gr[[compared_states]])[queryHits(fo_dmr_cnvcall)],
        comparison=as.factor(compared_states),
        short_id=as.factor(cf_id),
        delta_beta=mcols(dmr_gr[[compared_states]])$delta_beta[queryHits(fo_dmr_cnvcall)],
        gene=mcols(dmr_gr[[compared_states]])$gene_name_prot[queryHits(fo_dmr_cnvcall)],
        cn=mcols(cnvcall_gr[[cf_id]])$cn[subjectHits(fo_dmr_cnvcall)]) |>
        unite('chrom_start_comparison_gene', c(chrom, start, comparison, gene)))
  }
}
head(dmr_cn_tib)

# long-to-wide transformation for pheatmap
plot_tib <- pivot_wider(
  dmr_cn_tib,
  id_cols=chrom_start_comparison_gene,
  names_from=short_id,
  values_from=cn)

# it's pretty rare but there ARE rows in "plot_tib" that are majority NA.
# it's because some of the cancer samples have poor coverage, i.e., more
# genomic regions without copy number calls. those samples would have rank 51
# and/or rank 52 creeping up into the top 50 because `findOverlaps()` doesn't
# report no-overlaps -_____________-".
#
# one could assume NA calls could default to mean "copy number = 2"; but for
# cancer samples, i don't think that's a good assumption (for NDBE, sure--but
# not cancer). i rather nuke rows with majority NAs, then remaining rare NAs are
# shaded in grey in the heatmap later.
na_count <- rowSums(is.na(plot_tib[-1]))
na_count
plot_tib <- plot_tib[na_count < ncol(plot_tib) * .9, ]
head(plot_tib)
# 300 rows is correct! 6 comparisons x top 50

# convert to a matrix to feed into `pheatmap`
plot_matrix <- as.matrix(plot_tib[-1])
rownames(plot_matrix) <- plot_tib$chrom_start_comparison_gene
head(plot_matrix)
table(plot_matrix)  # confirms there's no CN == 0 (no dark red in plot)

# store metadata for comparison type
metadata_dmr <- data.frame(
  sub('.*_(.*? vs. .*?)_.*', '\\1', plot_tib$chrom_start_comparison_gene),
  row.names=rownames(plot_matrix))
colnames(metadata_dmr) <- c('comparison')

# highlight gene names with consistent CNV across all samples
# define "interesting" as CNV change (any direction) in > 20% of all samples
interesting_dmr <- rowSums(plot_matrix != 2, na.rm=TRUE) > 0.2 * ncol(plot_matrix)
plot_matrix[interesting_dmr, ]
apply(plot_matrix[interesting_dmr, ], 1, table)
dmr_labels <- sub('^.*_', '', rownames(plot_matrix))
dmr_labels[!interesting_dmr] <- ''

# ... and to defend against reviewer 2's, "why 20%?!", do a bit of stats across
# all DMRs to see overall % of samples having DMRs in CNVs
#
# long-to-wide transformation to make life easier
all_dmr_cn_tib <- pivot_wider(
  all_dmr_cn_tib,
  id_cols=chrom_start_comparison_gene,
  names_from=short_id,
  values_from=cn)
frac_dmr_cnv <- rowSums(all_dmr_cn_tib[-1] != 2, na.rm=TRUE) / rowSums(!is.na(all_dmr_cn_tib))
summary(frac_dmr_cnv)                                      # five number summary
quantile(frac_dmr_cnv, c(0.9, 0.95, 0.99, 0.999, 0.9999))  # the nines
# oh nice, the intuitive "> 20% samples means interesting" corresponds to the
# top 1% percentile of all DMRs

# okay, plot things out
#+ fig.width=10, fig.height=8
p <- pheatmap(plot_matrix,
              color=c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'),
              breaks=c(0, 0.8, 1.6, 2.4, 3.2, 4),
              cluster_rows=FALSE, cluster_cols=FALSE,
              annotation_row=metadata_dmr,
              fontsize=8,
              labels_row=dmr_labels)
ggsave('raw_fig3b.pdf', plot=p$gtable, width=8, height=6)

# for replicability purposes
sessionInfo()
