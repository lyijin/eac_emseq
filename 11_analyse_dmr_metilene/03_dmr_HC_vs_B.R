#!/usr/bin/env Rscript

"> 03_dmr_HC_vs_B.R <

Carry out a specific 2-way comparison for [HGD + EAC] vs. NDBO, so that results
can be compared to previous efforts e.g., Melissa Thomas' thesis.
" -> doc

suppressPackageStartupMessages({
  library(data.table)
  library(eulerr)
  library(GenomicRanges)
  library(ggplot2)
  library(this.path)
})


setwd(this.path::here())
METILENE_FILE <- '../10_call_dmr_metilene/metilene_output/allHC_vs_allB.dmr.tsv.gz'
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'
GENOMIC_ANNOT_FILE <- '../data/gencode.v43.annotation_2026-02-19.RData'

# downloaded from https://zwdzwd.github.io/InfiniumAnnotation, as an alternative
# to importing massive libraries...
HM450_HG38_FILE <- './data/HM450.hg38.manifest.gencode.v36.tsv.gz'  # coords are 0-based

# for comparison vs. thomas' thesis
THOMAS_HYPER_FILE <- './data/thomas.B_vs_HC.hyper.top100probes.txt'
THOMAS_HYPO_FILE <- './data/thomas.B_vs_HC.hypo.top100probes.txt'
THOMAS_AMPLICONS_FILE <- './data/thomas.B_vs_HC.amplicons.txt'

# error out if any of the files do not exist
stopifnot(file.exists(METILENE_FILE))
stopifnot(file.exists(CLINDATA_FILE))
stopifnot(file.exists(GENOMIC_ANNOT_FILE))
stopifnot(file.exists(HM450_HG38_FILE))
stopifnot(file.exists(THOMAS_HYPER_FILE))
stopifnot(file.exists(THOMAS_HYPO_FILE))
stopifnot(file.exists(THOMAS_AMPLICONS_FILE))

# load pre-processed genomic annotations and functions to annotate GRanges
load(GENOMIC_ANNOT_FILE)
source('../01_txdb/annotate_ranges_functions.R')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# read the metilene results file
dmr_dt <- fread(METILENE_FILE, sep='\t', header=FALSE,
                col.names=c('chr', 'start', 'stop', 'qval', 'delta_beta',
                            'nCpG', 'p_MWU', 'p_KS', 'mean_g1', 'mean_g2'))
diag_message(METILENE_FILE, ' has ',
             formatC(nrow(dmr_dt), big.mark=','), ' rows.')

# retain rows that are "great" DMRs
#   1. metilene has a min delta beta of 0.1 -- do nothing here
#      cutoff to select DMRs with larger changes in methylation levels
#   2. qval < 0.01 (while it removes heaps of predictions; qval 0.05 (147k) -->
#      0.01 (39k), but better DMRs with larger effect sizes (> 0.2) have
#      qval < 0.01 anyway
dmr_dt <- dmr_dt[qval < 0.01]
diag_message('Post qval filtering, ', formatC(nrow(dmr_dt), big.mark=','), ' rows remain.')

# quick scatterplot to check extent of delta_beta between B vs. HC
#+ fig.width=10, fig.height=3
ggplot(mapping=aes(dmr_dt$delta_beta)) +
  geom_histogram(binwidth=0.01) +
  ggtitle(METILENE_FILE) +
  theme_minimal(14)
# can't go too stringent, will end up with almost nothing if abs(delta) > 0.2
# proceed for now

# create GRanges from this data.table, and retain other cols as metadata
# start/stop coords from metilene are already 1-based, confirmed in JBrowse
dmr_gr <- makeGRangesFromDataFrame(
  dmr_dt[, c('chr', 'start', 'stop', 'qval', 'delta_beta', 'nCpG', 'mean_g1', 'mean_g2')],
  keep.extra.columns=TRUE,
  ignore.strand=TRUE)

# annotate these regions, and see top 10 most hypermeth regions
dmr_gr <- annotateDMRs(dmrs=dmr_gr, tx=gencode_tx, tss_proximal=2000)
dmr_gr[order(dmr_gr$delta_beta, decreasing=TRUE)[1:10], ]

# write all hypermeth regions into external file for eyeballing
tmp_gr <- dmr_gr[order(dmr_gr$delta_beta, decreasing=TRUE), ]
tmp_gr <- tmp_gr[tmp_gr$delta_beta > 0, ]
write.table(tmp_gr,
            file='./filtered_dmrs/allHC_vs_allB.hypermeth.tsv',
            quote=FALSE, sep='\t', row.names=FALSE)

# repeat for hypometh regions
tmp_gr <- dmr_gr[order(dmr_gr$delta_beta, decreasing=FALSE), ]
tmp_gr <- tmp_gr[tmp_gr$delta_beta < 0, ]
write.table(tmp_gr,
            file='./filtered_dmrs/allHC_vs_allB.hypometh.tsv',
            quote=FALSE, sep='\t', row.names=FALSE)


# prediction checking #1: check prediction overlap with thomas' thesis
#                         (100 hyper + 100 hypo)
hm450_dt <- fread(HM450_HG38_FILE, sep='\t', header=TRUE)
hm450_gr <- makeGRangesFromDataFrame(
  hm450_dt[, c('CpG_chrm', 'CpG_beg', 'CpG_end', 'probeID', 'genesUniq')],
  keep.extra.columns=TRUE,
  seqnames.field='CpG_chrm',
  start.field='CpG_beg',
  end.field='CpG_end',
  starts.in.df.are.0based=TRUE,
  na.rm=TRUE)
thomas_hyp_probes <- c(scan(THOMAS_HYPER_FILE, what=character()),
                       scan(THOMAS_HYPO_FILE, what=character()))
thomas_hyp_gr <- hm450_gr[hm450_gr$probeID %in% thomas_hyp_probes, ]
dmrs_with_hm450_obj <- findOverlaps(hm450_gr, dmr_gr)
dmrs_overlap_thomas_hyp_obj <- findOverlaps(thomas_hyp_gr, dmr_gr)

thomas_hm450_venn <- euler(
  c('Thomas'            = 200 - length(dmrs_overlap_thomas_hyp_obj),
    'This study'        = length(dmrs_with_hm450_obj) - length(dmrs_overlap_thomas_hyp_obj),
    'Thomas&This study' = length(dmrs_overlap_thomas_hyp_obj)))
plot(thomas_hm450_venn,
     labels=list(cex=1.5),
     quantities=list(type=c('counts'), cex=1.5),
     fill=c("lightgrey", "darkgrey"))
# while overlap doesn't look very extensive, consider that this study's DMRs are
# ~2% of HM450 probes, i.e., the expected overlap between this study and thomas'
# is ~2% = 4/200 probes. we're getting 25ish%, i.e., defo statistically
# significant if we did Fisher's exact

# force printing out of DMRs that overlaps with thomas'
as.data.frame(subsetByOverlaps(dmr_gr, thomas_hyp_gr))

# prediction checking #2: checking vs. thomas' amplicons
thomas_ampl_dt <- fread(THOMAS_AMPLICONS_FILE, sep='\t', header=FALSE,
                        col.names=c('genename', 'chrom_region', 'chr',
                                    'start', 'end', 'size', 'nCpG'))
thomas_ampl_gr <- makeGRangesFromDataFrame(
  thomas_ampl_dt, keep.extra.columns=TRUE, ignore.strand=TRUE)
thomas_ampl_gr
dmrs_overlap_thomas_ampl_obj <- findOverlaps(thomas_ampl_gr, dmr_gr)
dmrs_overlap_thomas_ampl_obj
# hmmmmmm. only one overlap. upon further inspection, thomas' amplicon regions
# do not contain any methylated CpGs in our dataset. was i too stringent in
# culling positions initially...?

# check #3: text search against reported biomarkers (from nieto 2018 review &
# thejaani)
REPORTED_BIOMARKERS <-
  unique(c(
    # barrett 1999
    'CDKN2A',
    # eads 2001
    'CDKN2A', 'ESR1', 'MYOD1', 'CALCA', 'MGMT', 'TIMP3',
    # schulmann 2005
    'RUNX3', 'HPP1', 'CDKN2A',
    # clement 2006
    'APC', 'TIMP3', 'TERT',
    # clement 2008
    'WIF1',
    # smith 2008
    'CDKN2A', 'RUNX3',
    # jin 2009
    'HPP1', 'CDKN2A', 'RUNX3', 'CDH13', 'TAC1', 'NELL1', 'AKAP12', 'SST',
    # wang 2009
    'APC', 'CDKN2A',
    # agarwal 2012
    'BCL11B', 'CCDC57', 'F10', 'CKB', 'ELAVL3', 'ZNF653', 'GPR177', 'HOXB7',
    'IGF1R', 'ITGB8', 'JARID', 'JUND', 'LAMA5', 'LOC55565', 'MGC35308', 'MMD2',
    'TAF10', 'TLX3', 'UBP1', 'WNK4', 'WWC1', 'ZBTB7B', 'ZNF358',
    # krause 2016 (supp table 5, top 100 sorted by e val
    'TRIB3', 'CAPN10', 'TRPM4', 'CAPN10', 'EPHB3', 'SND1', 'BTBD9', 'CAPN10',
    'UBR4', 'HCCA2', 'ASPSCR1', 'KLHL33', 'C21orf67', 'SORBS1', 'SRFBP1',
    'SNORD103A', 'ARHGEF10L', 'DHX32', 'CAPN10', 'KLHL33', 'AP2A2', 'FHAD1',
    'MFSD9', 'VGLL4', 'ELOVL6', 'SNORD103A', 'HDLBP', 'RPS15AP10', 'APBB2',
    'BAIAP2', 'LRRC8D', 'JDP2', 'AP2A2', 'ABHD6', 'SEPT9', 'NUMA1', 'C7orf50',
    'AP2A2', 'HMG20A', 'SLC25A12', 'TRPM4', 'SEPT9', 'SEPW1', 'FAM118A',
    'LRRFIP2', 'BAIAP2', 'RPS15AP10', 'KREMEN1', 'NUMA1', 'SH3TC1', 'C7orf29',
    'IFT140', 'JDP2', 'SNORD103A', 'XRCC3', 'MAN1C1', 'ZFYVE21', 'TEAD1',
    'CLDN14', 'GTF2H5', 'DIP2A', 'RNF144B', 'TSC2', 'DAG1', 'CDC42BPB', 'CASZ1',
    'NUMA1', 'UNC84A', 'RPH3AL', 'FAM60A', 'NFIX', 'TPPP', 'GNAZ', 'JDP2',
    'ARL8A', 'CNKSR3', 'BAIAP2', 'C15orf42', 'ACVR2B', 'CSK', 'CEP72', 'LRP5',
    'TOP2B', 'LRIG1',
    # thomas thesis 2018
    'AKAP12', 'APC', 'CDH13', 'DAPK1', 'GPX', 'GST', 'HPP1', 'MGMT', 'NELL1',
    'p16', 'CDKN2A', 'RPRM', 'PRDM2', 'RUNX3', 'SFRP', 'SOCS', 'SST', 'TAC1',
    'TERT', 'TIMP3', 'WIF1',
    # salta 2020
    'COL14A1', 'ZNF569',
    # pinto 2022
    'APC', 'CDKN2A', 'MGMT', 'TIMP3'))
dmr_match_reported_gr <- dmr_gr[dmr_gr$gene_name_prot %in% REPORTED_BIOMARKERS, ]
dmr_match_reported_gr
sort(unique(dmr_match_reported_gr$gene_name_prot))
plot(dmr_match_reported_gr$delta_beta)
# very few hypermeth markers, most are hypo
dmr_match_reported_gr[dmr_match_reported_gr$delta_beta > 0]


# for replicability purposes
sessionInfo()
