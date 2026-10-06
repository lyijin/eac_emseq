#!/usr/bin/env Rscript

"> 03_specific_meth_vs_expr.promoters.R <

Are DMRs disproportionally co-located in the promoter of genes that are also
differentially expressed?

Hypothesis is: HYPERmethylated DMRs should be enriched in promoters of
differentially UNDERexpressed genes; and vice versa, HYPOmeth-OVERexpressed

Constrain analysis to samples that had had RNA-seq and EM-seq data produced
(~50% of all samples). On a per-gene basis, check whether log expression
inversely correlates with untransformed methylation levels.

Only do this for promoters; gene bodies seem to be far weaker at correlating
with expression data (see script #02 in this folder).

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
  library(readxl)
  library(stringr)
  library(tidyr)
  library(this.path)
})

setwd(this.path::here())
COMPILED_BETA_FILE <- '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
COMPILED_TPM_FILE <- '../data/udumanne_etal_deg/counts_tpm.csv.gz'
DMR_FOLDER <- '../11_analyse_dmr_metilene/filtered_dmrs/'
DEG_FOLDER <- '../data/udumanne_etal_deg/'
LIT_METH_GENES <- '../data/supp_table_25.xlsx'

# error out if any of the files do not exist
stopifnot(file.exists(COMPILED_BETA_FILE))
stopifnot(file.exists(COMPILED_TPM_FILE))
stopifnot(file.exists(DMR_FOLDER))
stopifnot(file.exists(DEG_FOLDER))
stopifnot(file.exists(LIT_METH_GENES))

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

# some genes have multiple DMRs in its promoter--which is fine, but just
# remember that the analysis is on a per-gene (specifically, ENSEMBL ID) basis,
# and for each sample, the mean beta will be calculated from all CpGs across
# the multiple DMRs
annot_dmr_tib |> filter(comparison == 'NDBE vs. NSq') |> 
  filter(hyper_or_hypo == 'hypermeth') |> arrange(gene_id) |> head(n=10)
dmr_gr <- makeGRangesFromDataFrame(annot_dmr_tib, ignore.strand=TRUE)


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
annot_deg_tib <- annot_deg_tib |> filter(deg == TRUE)  # discard non-DEGs


# time to read in the smaller gene expression file, noting the sample IDs with
# RNA-seq data
tpm_tib <- read_csv(COMPILED_TPM_FILE, show_col_types=FALSE) %>%
  # sort by sample ID, except for the gene_id column
  select(gene_id, sort(colnames(.))) |>
  # remove version in "gene_id"
  mutate(gene_id=sub('\\..*', '', gene_id)) %>%
  # remove genes that are too lowly expressed
  filter(rowMeans(.[-1]) > 1) |>
  # filter for genes that are differentially expressed in at least one comparison
  filter(gene_id %in% annot_deg_tib$gene_id) %>%
  # do a base10 log1p [i.e., log10 (x+1)] transform of raw TPM data for
  # correlation with methylation data
  mutate(across(!gene_id, ~ log1p(.x)/log(10)))
head(tpm_tib)
# confirm vast majority of the remaining values have transformed TPMs > 0
table(tpm_tib[-1] > 0)


# then read in the chonky per-CpG beta values table
beta_tib <- read_tsv(COMPILED_BETA_FILE, show_col_types=FALSE, lazy=TRUE) %>%
  select(chr, start, end, sort(intersect(colnames(.), colnames(tpm_tib))))

# subselect "tpm_tib" by common samples
tpm_tib <- tpm_tib |>
  select(gene_id, sort(intersect(colnames(beta_tib), colnames(tpm_tib))))

# convert beta tibble to a GRanges object, then to reduce memory usage, retain
# CpG positions that are in at least one DMR
beta_gr <- makeGRangesFromDataFrame(
  beta_tib,
  keep.extra.columns=TRUE,
  ignore.strand=TRUE)
# note: select='first' because i care that a CpG is in at least one DMR (but not
# which DMR exactly). computationally quicker than the default of select='all'
fo_beta_dmr <- findOverlaps(beta_gr, dmr_gr, minoverlap=2L, select='first')
beta_gr <- beta_gr[!is.na(fo_beta_dmr), ]
beta_gr[, 1:5]  # visual check

# nuke the tibble to free memory
rm(beta_tib)


# i'm sure there's an arcane apply() that achieves the same outcome, but...
# for each comparison...
cor_tib <- NULL
for (co in comparison_order) {
  # for each hyper/hypo...
  for (hoh in unique(annot_dmr_tib$hyper_or_hypo)) {
    # get ENSEMBL IDs for genes with DMR which are also DEG 
    deg_co <- annot_deg_tib |> filter(comparison == co) |> pull(gene_id)
    comp_hyp_dmr_tib <- annot_dmr_tib |>
      filter(comparison == co) |>
      filter(hyper_or_hypo == hoh) |>
      filter(gene_id %in% deg_co) |>
      arrange(gene_id)
    
    # for every gene with DMR(s) and with expression data...
    for (gi in unique(intersect(comp_hyp_dmr_tib$gene_id, tpm_tib$gene_id))) {
      comp_hyp_dmr_gr <- makeGRangesFromDataFrame(
        comp_hyp_dmr_tib |> filter(gene_id == gi))
      fo_beta_chdt <- findOverlaps(beta_gr, comp_hyp_dmr_gr,
                                   minoverlap=2L, select='first')
      mean_betas <- colMeans(as.matrix(mcols(beta_gr[!is.na(fo_beta_chdt), ])))
      mean_exprs <- unlist(tpm_tib[tpm_tib$gene_id == gi, -1])
      
      # carry out correlation (pearson) between betas and exprs
      stopifnot(names(mean_betas) == names(mean_exprs))  # confirm samples share same order
      
      # "pearson" & "na.omit" are default params for `cor.test()`, but just be
      # explicit here
      ct <- cor.test(mean_betas, mean_exprs, method='pearson', na.action='na.omit')
      
      # append results into a single-row tibble first
      cor_row <- bind_cols(
        # descriptive stuff & cor results
        bind_rows(
          gene_id=gi,
          comparison=co,
          hyper_or_hypo=hoh,
          cor_pearson=unname(ct$estimate),
          cor_pval=ct$p.value,
        ),
        # keep beta values
        bind_rows(mean_betas) %>% rename_with(~ paste0('beta_', .)),
        # and also keep the tpm values
        bind_rows(mean_exprs) %>% rename_with(~ paste0('logtpm_', .))
      )
      
      # then append the row into the overall tibble
      cor_tib <- bind_rows(cor_tib, cor_row)
    }
  }
}

# add gene annotation
cor_tib <- left_join(
  cor_tib, 
  annot_dmr_tib |> select(gene_id, gene_name) |> distinct(),
  by='gene_id')

# perform BH multiple testing correction, tidy up table
cor_tib <- cor_tib |>
  arrange(cor_pval) |>
  mutate(cor_qval=p.adjust(cor_pval, method='BH')) |>
  relocate(cor_qval, .after=cor_pval) |>
  relocate(gene_name, .after=gene_id)
head(cor_tib)

# look at how the q values are distributed
nrow(cor_tib)
table(cut(cor_tib$cor_qval, breaks=c(0, 0.001, 0.01, 0.05, 1)))
# hmm, a more stringent cutoff of 0.01 is still okay, as it picks out 1.3k genes
# of 5.7k genes
cor_tib <- cor_tib |> 
  filter(cor_qval < 0.01) |>
  # for downstream plots
  mutate(comparison=factor(comparison, levels=comparison_order)) |>
  arrange(comparison, hyper_or_hypo)

# get overall stats for hypermethylated DMRs & anticorrelated meth-expr
cor_tib |> filter(hyper_or_hypo == 'hypermeth') |> filter(cor_pearson < 0) |>
  count(comparison)

# export "cor_tib" as a supplementary table
write_tsv(cor_tib, 'raw_supptable24.tsv')


# create plots to visually appraise genes with the strongest anticorrelations
plot_tib <- cor_tib |>
  # focus on genes that exhibit hypermethylation, and are -vely correlated
  # with expression (the expected biological relationship)
  filter(hyper_or_hypo == 'hypermeth') |>
  filter(cor_pearson < 0) |>
  # focus on genes with the lowest corrected p values
  group_by(comparison) |>
  slice_min(cor_qval, n=10) |>
  # don't really p/q values
  select(!ends_with('val')) |>
  # convert wide-to-long
  pivot_longer(
    cols=starts_with(c('beta_', 'logtpm_')),
    names_to=c('.value', 'short_id'),
    names_pattern='(.*?)_(.*$)') |>
  drop_na() |>
  # for nicer labeling: add disease classifications
  mutate(classification=factor(short_to_long_disease_vec[substr(short_id, 1, 1)],
                               levels=disease_progression_vec))
head(plot_tib)

#+ fig.width=10, fig.height=7
pearsonr_text_tib <- plot_tib |> 
  select(comparison, gene_name, cor_pearson) |> 
  distinct() |>
  mutate(cor_pearson=sprintf('%.3f', round(cor_pearson,3)))
ggplot(plot_tib, aes(x=beta, y=logtpm)) +
  geom_point(aes(color=classification), alpha=0.5) +
  geom_text(data=pearsonr_text_tib, 
            aes(x=Inf, y=Inf, label=paste0('r=', cor_pearson)), hjust=1, vjust=1) +
  scale_color_manual(values=disease_to_color_vec) +
  scale_x_continuous(labels=scales::number_format(accuracy=0.1), n.breaks=4) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), n.breaks=4) +
  facet_wrap(vars(comparison, gene_name), nrow=4, scales='free') +
  theme_minimal(12) +
  theme(legend.position='top',
        panel.grid.minor=element_blank())
ggsave('raw_suppfig7a.pdf', width=12, height=8)


# plot similar scatterplot for commercially interesting genes
interesting_meth_genes <- read_excel(
  LIT_METH_GENES, sheet=1, range='A13:A133') |>  # this is from a supp table
  drop_na() |>
  distinct(Gene) |>                         # get unique members
  filter_out(str_detect(Gene, "'|_")) |>    # remove genes with "'" or "_"
  pull(Gene)                                # converts tibble to vector
interesting_meth_genes

plot_tib <- cor_tib |>
  # focus on genes that exhibit hypermethylation, and are -vely correlated
  # with expression (the expected biological relationship)
  filter(hyper_or_hypo == 'hypermeth') |>
  filter(cor_pearson < 0) |>
  # focus on genes previously associated with NDBE/EAC
  filter(gene_name %in% interesting_meth_genes) |>
  # don't really p/q values
  select(!ends_with('val')) |>
  # convert wide-to-long
  pivot_longer(
    cols=starts_with(c('beta_', 'logtpm_')),
    names_to=c('.value', 'short_id'),
    names_pattern='(.*?)_(.*$)') |>
  drop_na() |>
  # for nicer labeling: add disease classifications
  mutate(classification=factor(short_to_long_disease_vec[substr(short_id, 1, 1)],
                               levels=disease_progression_vec))
head(plot_tib)

#+ fig.width=10, fig.height=2.5
pearsonr_text_tib <- plot_tib |>
  select(comparison, gene_name, cor_pearson) |>
  distinct() |>
  mutate(cor_pearson=sprintf('%.3f', round(cor_pearson,3)))
ggplot(plot_tib, aes(x=beta, y=logtpm)) +
  geom_point(aes(color=classification), alpha=0.5) +
  geom_text(data=pearsonr_text_tib, 
            aes(x=Inf, y=Inf, label=paste0('r=', cor_pearson)), hjust=1, vjust=1) +
  scale_color_manual(values=disease_to_color_vec) +
  scale_x_continuous(labels=scales::number_format(accuracy=0.1), n.breaks=4) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), n.breaks=4) +
  facet_wrap(vars(comparison, gene_name), nrow=1, scales='free') +
  theme_minimal(12) +
  theme(legend.position='top',
        panel.grid.minor=element_blank())
ggsave('raw_suppfig7b.pdf', width=12, height=2.7)


# for replicability purposes
sessionInfo()
