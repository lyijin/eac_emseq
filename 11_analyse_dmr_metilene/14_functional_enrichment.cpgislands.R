#!/usr/bin/env Rscript

"> 14_functional_enrichment.cpgislands.R <

Find out what the DMRs mean biologically; but restrict this analysis to CpG
islands as they are potentially(?) more interesting functionally.

To simplify plots and the narrative, only four pairwise comparisons is included
in this script, i.e., the progression-to-next-worst-pathology

  NDBE vs. NSq
  LGD vs. NDBE
  HGD vs. LGD
  EAC vs. HGD
" -> doc

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(dplyr)
  library(forcats)
  library(GenomicRanges)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(this.path)
})

setwd(this.path::here())
DMR_FOLDER <- './filtered_dmrs/'
LIT_METH_GENES <- '../data/supp_table_25.xlsx'
CPGISLAND_GR_FILE <- '../08_qc_post-sample_filter/island_beta_gr.rds'

# error out if any of the files do not exist
stopifnot(file.exists(DMR_FOLDER))

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


# read in locations of CpGs in CpG islands
cpgisland_gr <- readRDS(CPGISLAND_GR_FILE)

# focus on the hypermeth DMRs first, in the N B L H C sequence + HCvB + CvN
# (last two to fit what commercial tests are doing)
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
annot_tib <- tibble(
  seqnames=character(), start=numeric(), end=numeric(), width=numeric(),
  qval=numeric(), delta_beta=numeric(), nCpG=numeric(), mean_g1=numeric(),
  mean_g2=numeric(), gene_id_prot=character(), gene_name_prot=character(),
  promoter=numeric(), genebody=numeric(), intergenic=numeric(),
  comparison=character(), hyper_or_hypo=character(), rank=numeric())

for (dmr_file in dmr_files) {
  dmr_tib <- read_tsv(dmr_file, show_col_types=FALSE)
  diag_message(dmr_file, ' has ', formatC(nrow(dmr_tib), big.mark=','), ' rows.')
  
  # convert "dmr_tib" to a GRanges object so we can check overlaps
  dmr_gr <- makeGRangesFromDataFrame(
    dmr_tib,
    keep.extra.columns=FALSE,
    ignore.strand=TRUE)
  
  # "dmr_gr" are DMRs of width 200+ bp; "cpgisland_gr" are CpGs of width 2 bp.
  # be generous; whitelist a DMR as being in a CpG island if it overlaps ANY
  # known CpGs in CpG islands
  fo_dmr_cpgislands <- findOverlaps(dmr_gr, cpgisland_gr)
  dmr_tib <- dmr_tib |> slice(unique(queryHits(fo_dmr_cpgislands)))  # chop chop
  
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
  
  # add rank column, preserve how metilene ranks DMRs (by absolute delta_beta)
  dmr_tib$rank <- 1:nrow(dmr_tib)
  
  # rbind important columns into `annot_tib`
  annot_tib <- bind_rows(annot_tib, dmr_tib[colnames(annot_tib)])
}

# visual check that tibble's fine
annot_tib <- annot_tib |> 
  arrange(match(comparison, comparison_order)) |> 
  mutate(comparison=as_factor(comparison))
head(annot_tib)

# `clusterProfiler` wants ENTREZ gene IDs, need to use `AnnotationDbi` to
# map IDs across. gsub(blah) converts "ENSG12345678901.123" --> "ENSG12345678901"
annot_tib$gene_entrez_prot <- mapIds(org.Hs.eg.db,
                                     keys=gsub("\\.\\d+$", "", annot_tib$gene_id_prot, perl=TRUE),
                                     keytype='ENSEMBL', column='ENTREZID')
diag_message('There are ', sum(is.na(annot_tib$gene_entrez_prot)),
             ' DMRs where the ENSGs did not have an ENTREZ ID, out of ',
             nrow(annot_tib), ' DMRs.')
# hmm that's < 1%, can live with that. onwards!


# before we delve into making dot plots / cnetplots (network plots), note that
# cnetplot() has a "foldChange" argument to only plot a subset of interesting
# genes. designed more for transcriptomic data--allows viz of a subset of genes
# with diff expression greater than a certain threshold. this obviously does not
# work well on meth data (beta has no "fold change"). 
# 
# but i can artifically define "foldChange" and focus the network plot on two
# sets of relevant genes. either way, interesting genes have a foldChange
# value of 99999; while uninteresting genes are set to 1
#
#   1. novel DMR rank <= 50 in any pairwise comparison
dummy_hyper_fc_tib <- annot_tib |> 
  filter(hyper_or_hypo == 'hypermeth') |>
  select(rank, gene_name_prot) |>
  group_by(gene_name_prot) |>
  summarise_at(vars(rank), list(min_rank=min)) |>
  mutate(dummy_fc=if_else(min_rank <= 50, 99999, 1))
table(dummy_hyper_fc_tib$dummy_fc)

#   2. previously identified as a DMR in NDBE/EAC studies. subset of these
#      are commercially interesting too
interesting_meth_genes <- read_excel(
  LIT_METH_GENES, sheet=1, range='A13:A133') |>  # this is from a supp table
  drop_na() |>
  distinct(Gene) |>                         # get unique members
  filter_out(str_detect(Gene, "'|_")) |>    # remove genes with "'" or "_"
  pull(Gene)                                # converts tibble to vector
interesting_meth_genes
dummy_hyper_fc_tib[
  dummy_hyper_fc_tib$gene_name_prot %in% interesting_meth_genes, ]$dummy_fc <- 99999
table(dummy_hyper_fc_tib$dummy_fc)

# coax this into a variable that can be fed into the "foldChange" arg later
dummy_hyper_fc <- dummy_hyper_fc_tib$dummy_fc
names(dummy_hyper_fc) <- dummy_hyper_fc_tib$gene_name_prot

# repeat same exercise for hypomethylated genes
dummy_hypo_fc_tib <- annot_tib |> 
  filter(hyper_or_hypo == 'hypometh') |>
  select(rank, gene_name_prot) |>
  group_by(gene_name_prot) |>
  summarise_at(vars(rank), list(min_rank=min)) |>
  mutate(dummy_fc=if_else(min_rank <= 50, 99999, 1))
table(dummy_hypo_fc_tib$dummy_fc)

dummy_hypo_fc_tib[
  dummy_hypo_fc_tib$gene_name_prot %in% interesting_meth_genes, ]$dummy_fc <- 99999
table(dummy_hypo_fc_tib$dummy_fc)

dummy_hypo_fc <- dummy_hypo_fc_tib$dummy_fc
names(dummy_hypo_fc) <- dummy_hypo_fc_tib$gene_name_prot


# use clusterProfiler to carry out GO term & pathway enrichment
# code inspired by https://bioinformatics.ccr.cancer.gov/docs/btep-coding-club/CC2023/FunctionalEnrich_clusterProfiler/
# these steps take forever to run (~2 hours)
hyper_go_cclust <- compareCluster(
  gene_entrez_prot ~ comparison,
  data=annot_tib |> filter(hyper_or_hypo == 'hypermeth'),
  fun=enrichGO, OrgDb=org.Hs.eg.db, ont='BP')
simplified_hyper_go_cclust <- simplify(
  hyper_go_cclust, cutoff=0.7, by='p.adjust', select_fun=min, measure='Wang')
simplified_hyper_go_cclust <- setReadable(simplified_hyper_go_cclust, 'org.Hs.eg.db', 'ENTREZID')

#+ fig.width=10, fig.height=8
dotplot(simplified_hyper_go_cclust, showCategory=8)
cnetplot(simplified_hyper_go_cclust, showCategory=4, size_category=1.5, size_item=1,
         foldChange=dummy_hyper_fc, fc_threshold=2) +
  scale_fill_manual(values=c('#377eb8', '#4daf4a', '#a65628', '#e41a1c'))


# repeat process for GO terms for hypomethylated DMRs
hypo_go_cclust <- compareCluster(
  gene_entrez_prot ~ comparison,
  data=annot_tib |> filter(hyper_or_hypo == 'hypometh'),
  fun=enrichGO, OrgDb=org.Hs.eg.db, ont='BP')
simplified_hypo_go_cclust <- simplify(
  hypo_go_cclust, cutoff=0.7, by='p.adjust', select_fun=min, measure='Wang')
simplified_hypo_go_cclust <- setReadable(simplified_hypo_go_cclust, 'org.Hs.eg.db', 'ENTREZID')

#+ fig.width=10, fig.height=8
dotplot(simplified_hypo_go_cclust, showCategory=8)
cnetplot(simplified_hypo_go_cclust, showCategory=4, size_category=1.5, size_item=1,
         foldChange=dummy_hypo_fc, fc_threshold=2) +
  scale_fill_manual(values=c('#377eb8', '#4daf4a', '#a65628', '#e41a1c'))


# find KEGG pathways for hypermethylated genes. restrict to genes in the
# vicinity of the top-100 hypermethylated DMRs per category to reduce visual
# clutter
hyper_kegg_cclust <- compareCluster(
  gene_entrez_prot ~ comparison,
  data=annot_tib |> filter(hyper_or_hypo == 'hypermeth'),
  fun=enrichKEGG)
hyper_kegg_cclust <- setReadable(hyper_kegg_cclust, 'org.Hs.eg.db', 'ENTREZID')

#+ fig.width=10, fig.height=8
dotplot(hyper_kegg_cclust, showCategory=8) + 
  ggtitle('Enriched KEGG pathways for hypermethylated genes')

#+ fig.width=10, fig.height=10
cnetplot(hyper_kegg_cclust, showCategory=4, size_category=2, size_item=1,
         foldChange=dummy_hyper_fc, fc_threshold=2) +
  scale_fill_manual(values=c('#377eb8', '#4daf4a', '#a65628', '#e41a1c'))


# repeat KEGG network analysis for hypomethylated DMRs
hypo_kegg_cclust <- compareCluster(
  gene_entrez_prot ~ comparison,
  data=annot_tib |> filter(hyper_or_hypo == 'hypometh'),
  fun=enrichKEGG)
hypo_kegg_cclust <- setReadable(hypo_kegg_cclust, 'org.Hs.eg.db', 'ENTREZID')

#+ fig.width=10, fig.height=8
dotplot(hypo_kegg_cclust, showCategory=8) + 
  ggtitle('Enriched KEGG pathways for hypomethylated genes')
hypo_kegg_cclust  # confirm that HGD vs. LGD and EAC vs. HGD do not have
# enriched networks

#+ fig.width=10, fig.height=10
cnetplot(hypo_kegg_cclust, showCategory=4, size_category=1.5, size_item=1,
         foldChange=dummy_hypo_fc, fc_threshold=2) +
  scale_fill_manual(values=c('#377eb8', '#4daf4a', '#a65628', '#e41a1c'))


# for replicability purposes
sessionInfo()
