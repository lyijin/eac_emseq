#!/usr/bin/env Rscript

"> 01_overall_meth_vs_expr.promoters.R <

Are DMRs disproportionally co-located in the promoter of genes that are also
differentially expressed?

Hypothesis is: HYPERmethylated DMRs should be enriched in promoters of
differentially UNDERexpressed genes; and vice versa, HYPOmeth-OVERexpressed

This analysis works on overall bulk data, there are some samples where we only
have methylation/expression data, but for now let's just compare meth/expr means
across stages (as it's much easier to do than analysing individual samples).

To simplify plots and the narrative, only four pairwise comparisons is included
in this script, i.e., the progression-to-next-worst-pathology

  NDBE vs. NSq
  LGD vs. NDBE
  HGD vs. LGD
  EAC vs. HGD
" -> doc

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(ggpp)
  library(GenomicRanges)
  library(readr)
  library(tidyr)
  library(this.path)
})

setwd(this.path::here())
DMR_FOLDER <- '../11_analyse_dmr_metilene/filtered_dmrs/'
DEG_FOLDER <- '../data/udumanne_etal_deg/'

# error out if any of the files do not exist
stopifnot(file.exists(DMR_FOLDER))
stopifnot(file.exists(DEG_FOLDER))


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
  'N'='NSq', 'B'='NDBE', 'L'='LGD', 'H'='HGD', 'C'='EAC', 'HC'='HGD/EAC')
disease_progression_vec <- c('NSq', 'NDBE', 'LGD', 'HGD', 'EAC')
comparison_order <- c('NDBE vs. NSq', 'LGD vs. NDBE', 'HGD vs. LGD', 'EAC vs. HGD')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# focus on DMRs from pairwise comparisons in the N B L H C sequence
dmr_files <- c(paste0(DMR_FOLDER, 'allB_vs_allN.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allB_vs_allN.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allL_vs_allB.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allL_vs_allB.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allH_vs_allL.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allH_vs_allL.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allC_vs_allH.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allC_vs_allH.hypometh.tsv.gz'))
stopifnot(file.exists(dmr_files))

# create a dt to store annots of DMRs
annot_dmr_tib <- tibble(
  seqnames=character(), start=numeric(), end=numeric(), width=numeric(),
  qval=numeric(), delta_beta=numeric(), nCpG=numeric(), mean_g1=numeric(),
  mean_g2=numeric(), gene_id=character(), gene_name=character(),
  promoter=numeric(), genebody=numeric(), intergenic=numeric(),
  comparison=character(), hyper_or_hypo=character())

for (dmr_file in dmr_files) {
  dmr_tib <- read_tsv(dmr_file, show_col_types=FALSE)
  diag_message(dmr_file, ' has ', formatC(nrow(dmr_tib), big.mark=','), ' rows.')
  
  # filter for promoters! only pick DMRs that sit COMPLETELY WITHIN promoter
  # regions (1.0). partial overlaps result in fractional values
  dmr_tib <- dmr_tib |> filter(promoter == 1)
  diag_message(dmr_file, ' has ', formatC(nrow(dmr_tib), big.mark=','), ' promoter rows.')
  
  # parse filename for contextual clues for sample types / risk states
  fname_split <- unlist(strsplit(dmr_file, '/', fixed=TRUE))
  fname_split <- unlist(strsplit(fname_split[length(fname_split)], '.', fixed=TRUE))
  hyper_or_hypo <- fname_split[2]
  compared_states <- fname_split[1]
  compared_states <- unlist(strsplit(compared_states, '_vs_'))
  diag_message(dmr_file, ' are ', hyper_or_hypo, ' markers for ',
               compared_states[1], ' vs ', compared_states[2])
  
  # for this plot, drop the "all" label for stuff like "allN"
  compared_states <- sub('all', '', compared_states)
  compared_states <- paste0(short_to_long_disease_vec[compared_states[1]], ' vs. ',
                            short_to_long_disease_vec[compared_states[2]])
  
  dmr_tib$comparison <- compared_states
  dmr_tib$hyper_or_hypo <- hyper_or_hypo
  
  # rbind important columns into `annot_dmr_tib`
  annot_dmr_tib <- bind_rows(annot_dmr_tib, dmr_tib[colnames(annot_dmr_tib)])
}
annot_dmr_tib <- annot_dmr_tib |>
  # remove version in "gene_id"
  mutate(gene_id=sub('\\..*', '', gene_id)) |>
  # tidy table
  arrange(match(comparison, comparison_order)) |> 
  # convert comparison to a factor
  mutate(comparison=as_factor(comparison))
head(annot_dmr_tib)

# some genes have multiple DMRs in its promoter
annot_dmr_tib |> filter(comparison == 'NDBE vs. NSq') |> 
  filter(hyper_or_hypo == 'hypermeth') |> arrange(gene_id) |> head(n=10)
# calculate the WEIGHTED MEAN of delta beta, weighted by # of CpGs in the DMRs
wm_dmr_tib <- annot_dmr_tib |>
  summarize(wm_delta_beta=weighted.mean(delta_beta, nCpG),
            .by=c(comparison, hyper_or_hypo, gene_id)) |>
  # left-join to pull important info from "annot_dmr_tib"
  left_join(annot_dmr_tib |> select(gene_id, gene_name) |> distinct(), by='gene_id')
head(wm_dmr_tib)
table(wm_dmr_tib$comparison, wm_dmr_tib$hyper_or_hypo)  # of diff meth genes


# read the pairwise comparisons for DEGs. note that to qualify as a DEG, the
# gene has to have a padj < 0.05 and |log2FC| > 1 in BOTH files. the one without
# a 'cov_' prefix is for a simple model; the one with a 'cov_' prefix is
# covariate-adjusted (covariates were immune and stromal scores of the samples;
# Udumanne, personal communication)
deg_files <- c(paste0(DEG_FOLDER, 'im_nsq_deg.csv.gz'),
               paste0(DEG_FOLDER, 'cov_im_nsq_deg.csv.gz'),
               paste0(DEG_FOLDER, 'lgd_im_deg.csv.gz'),
               paste0(DEG_FOLDER, 'cov_lgd_im_deg.csv.gz'),
               paste0(DEG_FOLDER, 'hgd_lgd_deg.csv.gz'),
               paste0(DEG_FOLDER, 'cov_hgd_lgd_deg.csv.gz'),
               paste0(DEG_FOLDER, 'cancer_hgd_deg.csv.gz'),
               paste0(DEG_FOLDER, 'cov_cancer_hgd_deg.csv.gz'))
stopifnot(file.exists(deg_files))
annot_deg_tib <- tibble()
# argh, too complicated to write fancy code for this. make sure the DEG files
# are ordered in the same N B L H C pairwise comparison, then rip the
# comparison labels from "annot_dmr_tib"
compared_states <- levels(annot_dmr_tib$comparison)

for (n in seq(1, 7, 2)) {
  simple_deg_file <- deg_files[n]
  covadj_deg_file <- deg_files[n+1]
  
  simple_deg_tib <- read_csv(simple_deg_file, show_col_types=FALSE) |>
    dplyr::rename(gene_id=`...1`)
  covadj_deg_tib <- read_csv(covadj_deg_file, show_col_types=FALSE) |>
    dplyr::rename(gene_id=`...1`)
  
  # row orders & nrows of both tibbles should be identical, but to be safe, use
  # "inner_join()" to make sure only rows in both tibbles are combined 
  deg_tib <- inner_join(simple_deg_tib, covadj_deg_tib,
                        by='gene_id', suffix=c('', '_cov')) |>
    # remove rows with any NA. NAs in padj are from DESeq2's filtering, driven
    # by low or outlier counts
    drop_na() |>
    # remove version numbers from ENSGs (e.g., "ENSGx.11" --> "ENSGx")
    mutate(gene_id=sub('\\..*', '', gene_id)) |>
    # then apply cutoffs to find DEGs
    mutate(deg=abs(log2FoldChange) > 1 & abs(log2FoldChange_cov) > 1 &
             padj < 0.05 & padj_cov < 0.05) |>
    # and slap the correct comparison label
    mutate(comparison=compared_states[(n+1)/2])
  head(deg_tib)
  
  # select only important columns to keep for downstream analysis
  deg_tib <- deg_tib |> 
    select(c('gene_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'padj', 'deg', 'comparison'))
  
  diag_message('There are ', sum(deg_tib$deg), ' DEGs (out of ', nrow(deg_tib),
               ') for ', compared_states[(n+1)/2], '.')
  
  # rbind important columns into `annot_deg_tib`
  annot_deg_tib <- bind_rows(annot_deg_tib, deg_tib)
}
head(annot_deg_tib)
table(annot_deg_tib$comparison, annot_deg_tib$deg)  # show # of diff expr genes


# inner-join DMR info with DEG fold change (i.e., DMRs in promoters with no
# gene expression values will be dropped--i mean, not like NAs can be plotted
# right, why keep em?)
dmr_deg_tib <- inner_join(
  wm_dmr_tib |> select(c(gene_id, gene_name, comparison, wm_delta_beta)),
  annot_deg_tib |> select(c(gene_id, comparison, log2FoldChange, deg)),
  by=c('gene_id', 'comparison'))

dmr_deg_tib <- dmr_deg_tib |>
  # drop non-DEG
  filter(deg == TRUE) |>
  # convert comparison into a factor (suppresses alphabetical sorting later)
  mutate(across(comparison, as_factor))
head(dmr_deg_tib)

#+ fig.width=6, fig.height=6
# scatterplot of delta beta of DMR vs. log2 fold change of associated DEG
ggplot(dmr_deg_tib, aes(x=wm_delta_beta, y=log2FoldChange)) +
  # scatterplot
  geom_point(size=1.5, alpha=0.1) +
  # count # of points in every quadrant
  geom_quadrant_lines(color='#de2d26') +
  stat_quadrant_counts(color='#de2d26') +
  facet_wrap(. ~ comparison, nrow=2, ncol=2) +
  labs(title='Promoters: differential methylation vs. gene expression',
       x='Methylation: delta beta',
       y='Expression: log2 FC') +
  theme_minimal(12)
ggsave('raw_suppfig6.pdf', width=6, height=6)
# yes, the quadrants fitting the standard hypothesis (top-left + bottom-right)
# has more datapoints than the other two combined--across all four pairwise
# comparisons, but the numbers aren't spectacularly different really. hmm


# for replicability purposes
sessionInfo()
