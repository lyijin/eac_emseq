#!/usr/bin/env Rscript

"> 03_plot_contexts_with_excl.R <

What sort of genomic contexts houses the most variable CpGs?

After considering exploratory RNA-seq and EM-seq results, some samples look
obviously misclassified. Decision was made to exclude these samples. To reduce
duplication in code, this exclusion is carried out using a function residing
in `00_common/` instead of in individual scripts--and this reduces re-editing
of scripts if the excluded samples were to change.

Do trends hold higher biological significances if we restricted the analyses
to the top-X most variable positions (a la the 'top-5000-most-variable-probes')
hierarchical clustering & analyses done on many array-based papers (e.g.,
Yu, Gut, 2019; Jammula, Gastroenterology, 2020)?

In our case, as we have 25m CpGs, we opted for the most variable 100,000
CpGs, representing ~0.4% of the dataset. The % magnitude is in the same order
as subsetting the top 5k most variable probes in Illumina 450k arrays
(5k/450k = 1.1%) or MethylationEPIC arrays (5k/850k = 0.59%).

Similarly, do trends have higher biological significances if we restricted the
analyses to CpG islands only (like what Krause, Carcinogenesis, 2016 did on
ARRAY data)?
" -> doc

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(cowplot)
  library(data.table)
  library(forcats)
  library(GenomicRanges)
  library(ggplot2)
  library(this.path)
})


setwd(this.path::here())
COMPILED_BETA_FILE <- '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
GENOMIC_ANNOT_FILE <- '../data/gencode.v43.annotation_2026-02-19.RData'
SNP_ANNOT_FILE <- '../data/snp151Common_hg38.txt.gz'
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


# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}

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
    genomic_feature_as_long_dt(beta_gr[beta_gr$CpGisland == 1], 'CpG island'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$CpGshores == 1], 'CpG shores'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$nonCpG == 1], 'Non-island/shore'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$promoter == 1], 'Promoter'),
    genomic_feature_as_long_dt(beta_gr[beta_gr$genebody == 1], 'Gene body'),
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
beta_dt$variance_noN <-
  apply(beta_dt[, .SD, .SDcols=grepl('^B|^L|^H|^C', colnames(beta_dt))], 1, var, na.rm=TRUE)
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
  beta_dt[, c('chr', 'start', 'end', 'variance', 'variance_noN',
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
#+ fig.width=10, fig.height=7
plot_trends_by_annots(beta_gr)

# save in more compact form for manuscript
ggsave('raw_fig1e.pdf', width=5, height=5)


#+ fig.width=10, fig.height=7
# let's now look at the top N most variable positions
most_variable_cpgs <- order(beta_gr$variance, decreasing=TRUE)[1:TOP_VAR_POS]

# confirm variable is, indeed, storing the most variable CpGs
beta_gr[most_variable_cpgs[1:5]]
summary(beta_gr$variance)  # max value should be in the GRanges object slice

# check whether most variable CpGs are disproportionately in some chromosomes
# (guard against imprinting artefacts)
table(seqnames(beta_gr))
table(seqnames(beta_gr[most_variable_cpgs]))
# hmm, X has ~3x overrepresentation relative to expectation based on chromosomal
# lengths, but not to an extent where it'd distort overall trends. but anyway,
# to be safe, remove XYM
beta_gr <- beta_gr[!seqnames(beta_gr) %in% c('chrX', 'chrY', 'chrM')]
table(seqnames(beta_gr))

# check whether most variable CpGs disproportionately coincide with SNPs
# read SNP file, and remove intersection with "beta_gr"
snp_dt <- fread(SNP_ANNOT_FILE, select=1:4, sep='\t', header=FALSE,
                col.names=c('chrom', 'start', 'end', 'rsID'))
snp_gr <- makeGRangesFromDataFrame(snp_dt, keep.extra.columns=TRUE)
# drop alt chromosomes
seqlevels(snp_gr, pruning.mode='coarse') <- seqlevels(beta_gr)

# get overlaps, then remove CpGs that overlap *any* SNPs
ov <- findOverlaps(beta_gr, snp_gr, type='any', select='all')
diag_message('There are ', length(unique(queryHits(ov))),
             ' CpGs overlapping common SNPs (dbSNP151).')
beta_gr <- beta_gr[-unique(queryHits(ov))]

# do most variable positions tend to be associated with certain genomic loci?
most_variable_cpgs <- order(beta_gr$variance, decreasing=TRUE)[1:TOP_VAR_POS]
plot_trends_by_annots(beta_gr[most_variable_cpgs])
# hmm, median betas tend to hover around the 0.5 mark. makes mathematical sense
# (betas that are near-0 or near-1 has only one direction for values to deviate
# away from; betas at 0.5 have two directions for deviations, leading to larger
# variances), but this is not biologically interesting

# are trends similar if we focused on non-N variance? (B --> C)
most_variable_noN_cpgs <- order(beta_gr$variance_noN, decreasing=TRUE)[1:TOP_VAR_POS]
plot_trends_by_annots(beta_gr[most_variable_noN_cpgs])
# ... hmm, very extensive overlap


# if we shift focus to looking at the n values of CpG islands vs. shores vs.
# non-island/shores--the plain jane "top 100k most variable" positions had
# < 1% CpG islands; whilst the redefinitions of most variable as the largest
# differences between N/B vs. C led to higher % of positions in CpG islands.
# hints at CpG islands being biologically more interesting and/or relevant
# to disease progression!
#
# speculation: one reason why top-X-most-variable isn't working as well in our
# methylation dataset could be because we have five states (NBLHC); whilst
# vast majority of other EAC papers only have two (NC or BC). i think it's
# much easier to get probes that go low/high across two states than across five.
#
# there are 1,789,643 CpGs in CpG islands. applying the ballpark top ~0.5% idea,
# this results in ~10k positions (10000/1789643 = 0.59%)
TOP_VAR_POS <- 10000

# note that "beta_gr" is getting subsetted so it can be saved later
#+ fig.width=10, fig.height=6
beta_gr <- beta_gr[beta_gr$CpGisland == 1]
plot_trends_by_annots(beta_gr)

most_variable_cpgs <- order(beta_gr$variance, decreasing=TRUE)[1:TOP_VAR_POS]
plot_trends_by_annots(beta_gr[most_variable_cpgs])

most_variable_noN_cpgs <- order(beta_gr$variance_noN, decreasing=TRUE)[1:TOP_VAR_POS]
plot_trends_by_annots(beta_gr[most_variable_noN_cpgs])
# huh, surprisingly similar!

# check overlap of most variable CpGs
table(most_variable_cpgs %in% most_variable_noN_cpgs)
# both sets are equally sized (10k), FALSEs mean # CpGs unique to either
# selection criteria. overlap _is_ pretty extensive

# save "beta_gr" for PCA/hierarchical clustering of island-specific CpGs
saveRDS(beta_gr, 'island_beta_gr.rds')


# for replicability purposes
sessionInfo()
