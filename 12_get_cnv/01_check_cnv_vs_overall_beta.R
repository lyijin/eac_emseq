#!/usr/bin/env Rscript

"> 01_check_cnv_vs_overall_beta.R <

Do CNVs alter beta values in NSq-to-EAC progression?
" -> doc

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(forcats)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
  library(this.path)
  library(tidyr)
})


setwd(this.path::here())
COMPILED_BETA_FILE <- '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
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


# read compiled cov data, omit first three columns (chr/start/end)
# then sort by column name
beta_dt <- fread(COMPILED_BETA_FILE, sep='\t', header=TRUE)
# N samples do not, by definition, have CNVs (coverages from N samples form
# the baseline where all other samples are compared against)
beta_dt[, grep('^N', colnames(beta_dt)) := NULL]
# sort columns alphabetically, keep chr/start/end at the start
setcolorder(beta_dt, order(names(beta_dt)))
setcolorder(beta_dt, c('chr', 'start', 'end'))
beta_dt[90:95, 30:39]  # just to confirm presence of NA

# how many CpGs are there...?
nrow(beta_dt)

# confirm samples have been excluded. n should be 76 + 3 (chr/start/end) here
ncol(beta_dt)

# convert to GenomicRanges
beta_gr <- makeGRangesFromDataFrame(
  beta_dt,
  keep.extra.columns=TRUE,
  ignore.strand=TRUE)
beta_gr


# read per-sample inferred CNV values from `cnvkit.py call`
cnvcall_files <- c(Sys.glob('03_cnvkit_call_cnv/B*.call.cns'),
                   Sys.glob('03_cnvkit_call_cnv/L*.call.cns'),
                   Sys.glob('03_cnvkit_call_cnv/H*.call.cns'),
                   Sys.glob('03_cnvkit_call_cnv/C*.call.cns'))
summary_tib <- NULL
for (cf in cnvcall_files) {
  # get the three-char short_id from the filename
  cf_id <- gsub('.*/', '', gsub('\\..*$', '', cf))
  
  # read the CNV file, convert into GenomicRanges object for `findOverlaps()`
  cnvcall_gr <- makeGRangesFromDataFrame(
    fread(cf, select=c(1,2,3,8)),
    keep.extra.columns=TRUE,
    ignore.strand=TRUE)
  cnvcall_gr
  
  # when searching for overlaps, enforce complete overlap (prevents edge case
  # where CpG is straddles the boundaries of two 20 kb CNV buckets
  fo_beta_cnvcall <- findOverlaps(beta_gr, cnvcall_gr, minoverlap=2L)
  beta_cn_tib <- tibble(
    beta=mcols(beta_gr)[[cf_id]][queryHits(fo_beta_cnvcall)],
    cn=mcols(cnvcall_gr)$cn[subjectHits(fo_beta_cnvcall)],
    short_id=as.factor(cf_id)) |>
    drop_na() |> # some beta values can be NA
    mutate(cn_type=cut(cn, c(-1,2,3,999), right=FALSE,
                       labels=c('CN < 2', 'CN = 2', 'CN > 2')))
  summary_tib <- bind_rows(
    summary_tib,
    beta_cn_tib |>
      group_by(short_id, cn_type) |>
      summarize(mean=mean(beta),
                q1=quantile(beta, 0.25),
                median=median(beta),
                q3=quantile(beta, 0.75),
                .groups='drop'))
}

# add in another row that calculates mean across all samples, to see overall
# trends
summary_tib <- bind_rows(
  summary_tib |>
    group_by(cn_type) |>
    summarize(mean=mean(mean),
              q1=mean(q1, 0.25),
              median=mean(median),
              q3=mean(q3, 0.75),
              .groups='drop') |>
    mutate(short_id=as.factor('Overall mean')),
  summary_tib)

# have a tibble to store values for dotted lines / solid lines later
line_minmax_tib <- summary_tib |>
  group_by(short_id) |>
  summarize(start=min(mean), end=max(mean))

#+ fig.width=8, fig.height=10
ggplot(summary_tib, aes(x=mean, y=fct_rev(fct_inorder(short_id)), color=cn_type)) +
  geom_segment(data=line_minmax_tib,
               aes(x=0.50, xend=start, y=fct_rev(fct_inorder(short_id)), yend=fct_rev(fct_inorder(short_id))),
               color='#999999', linetype='dotted', linewidth=0.5, alpha=0.5) +
  geom_segment(data=line_minmax_tib,
               aes(x=start, xend=end, y=fct_rev(fct_inorder(short_id)), yend=fct_rev(fct_inorder(short_id))),
               color='#999999', linewidth=1, alpha=0.5) +
  geom_point(size=4, alpha=1) +
  geom_text_repel(aes(label=sprintf('%.2f', mean))) +
  scale_color_manual(values=c(`CN < 2`='#0571b0', `CN = 2`='#666666', `CN > 2`='#ca0020')) +
  labs(color='Copy number',
       x='Mean \u03b2',
       y='') +
  theme_minimal(12) +
  theme(legend.position='top',
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank())
ggsave('raw_suppfig5.pdf', width=8, height=10)


# for replicability purposes
sessionInfo()
