#!/usr/bin/env Rscript

"> 04_check_focal_cnv.R <

Are there focal CNVs that regularly occur in particular disease classifications?

Copy numbers of segments (per-sample) are in `03_cnvkit_call_cnv/*.call.cns`
files. Contents of these files are not directly comparable against each other,
as the coords/linecounts are not identical across all the files. Despite this,
these segments are made up of more granular segments of ~20 kb, coords
of which are in the `01_run_cnvkit/*.cnr` files.

NOTE: all cnr files have the same segments coordinates, i.e., first three
columns are completely identical across all files and line counts are identical.
  $ wc -l *.cnr | head -5
     143588 B02.cnr
     143588 B03.cnr
     143588 B04.cnr
     143588 B06.cnr
     143588 B08.cnr
  $ for a in *.cnr; do diff <(cut -f 1-3 ${a}) <(cut -f 1-3 B02.cnr); done
     (null output)

Script works by initialising a GenomicRanges of 143,588 segments corresponding
to the coords of a single `*.cnr` file, then iterates over the `*.call.cns`
files and deduces the copy number of each individual segment.

REMEMBER: cnvkit coords are 0-based!
  https://cnvkit.readthedocs.io/en/stable/fileformats.html#
" -> doc

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(ggplot2)
  library(pheatmap)
  library(readr)
  library(stringr)
  library(this.path)
  library(tidyr)
})


setwd(this.path::here())
CNR_FILE <- './01_run_cnvkit/B02.cnr'
CNS_FILES <- c(Sys.glob('./03_cnvkit_call_cnv/B*.call.cns'),
               Sys.glob('./03_cnvkit_call_cnv/L*.call.cns'),
               Sys.glob('./03_cnvkit_call_cnv/H*.call.cns'),
               Sys.glob('./03_cnvkit_call_cnv/C*.call.cns'))

# error out if any of the files do not exist
stopifnot(file.exists(CNR_FILE))
stopifnot(file.exists(CNS_FILES))

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


# read the cnr file to get the coords for the 143,588 segments
cnr_gr <- makeGRangesFromDataFrame(
  read_tsv(CNR_FILE, col_select=c(1,2,3,4), show_col_types=FALSE),
  keep.extra.columns=TRUE,
  ignore.strand=TRUE,
  starts.in.df.are.0based=TRUE)
# create a master tibble to store the copy number data, much easier to work on
# tibbles than mdata(cnr_gr)
cn_tib <- as_tibble(cnr_gr)

# start populating "cnv_gr" with copy number parsed from cns files
for (cf in CNS_FILES) {
  # parse filename for short ID
  fname_split <- unlist(strsplit(cf, '/', fixed=TRUE))
  fname_split <- unlist(strsplit(fname_split[length(fname_split)], '.', fixed=TRUE))
  short_id <- fname_split[1]
  diag_message('Reading ', cf, '...')
  
  # read the individual cns file for cnvkit-called copy number
  cns_gr <- makeGRangesFromDataFrame(
    read_tsv(cf, col_select=c(1,2,3,8), show_col_types=FALSE),
    keep.extra.columns=TRUE,
    ignore.strand=TRUE,
    starts.in.df.are.0based=TRUE)
  
  # check overlaps between master cnr and per-sample cns
  fo_cnr_cns <- findOverlaps(cnr_gr, cns_gr, type='within')
  
  # create null array of nrows matching the master tibble
  sample_cn_array <- rep(NA, nrow(cn_tib))
  
  # populate array with copy number from subjectHits
  sample_cn_array[queryHits(fo_cnr_cns)] <- mcols(cns_gr)$cn[subjectHits(fo_cnr_cns)]
  
  # sanity check, how many NAs are there? fill NAs with 2
  diag_message(short_id, ' has ', sum(is.na(sample_cn_array)),
               ' segments with no copy number data. Treat them as "2" (no CNV).')
  sample_cn_array[is.na(sample_cn_array)] <- 2
  
  # shove it into the master tibble
  cn_tib[short_id] <- sample_cn_array
}

# eyeball tibble
head(cn_tib)

# light-touch filtering
cn_tib <- cn_tib |>
  # remove strand info
  select(-strand) |>
  # ignore CNVs in X/Y, dataset is overwhelmingly male
  filter(seqnames != 'chrX') |>
  filter(seqnames != 'chrY') |>
  # noticed that some tiny windows are sneaking in (e.g., width=23; SERF1B
  # chr5:70037687..70037709). set a minimum 10 kb window threshold
  filter(width > 10000)

# tally # of samples with CN = 0, CN = 1, ..., CN = 4
n_samples <- cn_tib |> select(matches('^[A-Z]', ignore.case=FALSE)) |> ncol()
cn_tib <- bind_cols(
  cn_tib,
  cn_tib |> select(matches('^[A-Z]', ignore.case=FALSE)) %>%
    mutate(cn0 = rowSums(. == 0),
           cn1 = rowSums(. == 1),
           cn2 = rowSums(. == 2),
           cn3 = rowSums(. == 3),
           cn4 = rowSums(. == 4)) |>
    select(matches('cn\\d')))
# add columns to tally number of samples with deletions and amplifications
cn_tib <- cn_tib |>
  mutate(cndel=cn0+cn1,
         cnamp=cn3+cn4)

# check # windows with consistent del/amp & some quick maffs
consistent_del_windows <- cn_tib |> filter(cndel > 0.3 * n_samples) |> nrow()
diag_message('There are ', consistent_del_windows, ' consistent deleted windows ',
             '(20 kb), out of a universe of ', nrow(cn_tib), ' windows (',
             round(consistent_del_windows/nrow(cn_tib)*100, 2), '%).')
consistent_amp_windows <- cn_tib |> filter(cnamp > 0.3 * n_samples) |> nrow()
diag_message('There are ', consistent_amp_windows, ' consistent amplified windows ',
             '(20 kb), out of a universe of ', nrow(cn_tib), ' windows (',
             round(consistent_amp_windows/nrow(cn_tib)*100, 2), '%).')


# exploratory plot to find # of 20 kb windows with consistent focal deletions;
# define "consistent" as > 30% (i.e., n >= 17 out of 54 samples)
ggplot(cn_tib, aes(x=cndel)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=0.3*n_samples, color='#e41a1c') +
  annotate('label', label=paste0(
    '30% consistency cutoff\n', consistent_del_windows, ' windows (',
    round(consistent_del_windows/nrow(cn_tib)*100, 2), '% all windows)'),
           x=0.3*n_samples, y=Inf, hjust=0, vjust=1, color='#e41a1c') +
  coord_cartesian(xlim=c(0, n_samples)) +
  labs(title='Histogram of genomic windows with consistent focal deletions',
       x='# samples with focal deletion at same genomic window',
       y='Count of 20 kb genomic windows') +
  theme_minimal(12)

ggplot(cn_tib, aes(x=cnamp)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=0.3*n_samples, color='#e41a1c') +
  annotate('label', label=paste0(
    '30% consistency cutoff\n', consistent_amp_windows, ' windows (',
    round(consistent_amp_windows/nrow(cn_tib)*100, 2), '% all windows)'),
           x=0.3*n_samples, y=Inf, hjust=0, vjust=1, color='#e41a1c') +
  coord_cartesian(xlim=c(0, n_samples)) +
  labs(title='Histogram of genomic windows with consistent focal amplifications',
       x='# samples with focal amplification at same genomic window',
       y='Count of 20 kb genomic windows') +
  theme_minimal(12)


# what are the genes that are in focal del windows
del_tib <- cn_tib |>
  filter(cndel > 0.3 * n_samples) |>
  # no gene name, no play
  filter(gene != '-') |>
  group_by(gene) |>
  filter(cndel == max(cndel)) %>%
  # only select the first 20 kb window that matches the max() condition, as
  # multiple rows could potentially match the same condition (uninteresting)
  slice_head(n = 1) |>
  ungroup() |>
  arrange(-cndel)

# save this as a supp table before filtering/plotting is carried out
write_tsv(del_tib, 'raw_supptable23a.tsv')

# eyeball gene list
del_tib$gene
# there are some gene families known for (natural) copy number variation in the
# list: use a quick regex to remove
#   - "^OR": olfactory receptor genes
#   - "^POTE": POTE genes (Prostate, Ovary, Testis-Expressed Protein)
#   - "^IFNA": interferon alpha family
del_tib <- del_tib |>
  filter(str_detect(gene, '^OR', negate=TRUE)) |>
  filter(str_detect(gene, '^POTE', negate=TRUE)) |>
  filter(str_detect(gene, '^IFNA', negate=TRUE))

# plot this as a heatmap
del_df <- del_tib |>
  select(matches('^[A-Z]', ignore.case=FALSE)) |>
  as.matrix()
rownames(del_df) <- del_tib$gene

#+ fig.width=8, fig.height=3
pd <- pheatmap(
  del_df,
  color=c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'),
  border_color=NA,
  breaks=c(0, 0.8, 1.6, 2.4, 3.2, 4),
  cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize=7,
  gaps_col=c(24,36,47))


# then do the same for focal amps
amp_tib <- cn_tib |>
  filter(cnamp > 0.3 * n_samples) |>
  # no gene name, no play
  filter(gene != '-') |>
  group_by(gene) |>
  filter(cnamp == max(cnamp)) %>%
  # only select the first 20 kb window that matches the max() condition, as
  # multiple rows could potentially match the same condition (uninteresting)
  slice_head(n = 1) |>
  ungroup() |>
  arrange(-cnamp)

# save this as a supp table before filtering/plotting is carried out
write_tsv(amp_tib, 'raw_supptable23b.tsv')

# similar process as genes involved in focal deletion
amp_tib$gene
amp_tib <- amp_tib |>
  filter(str_detect(gene, '^OR', negate=TRUE)) |>
  filter(str_detect(gene, '^POTE', negate=TRUE)) |>
  filter(str_detect(gene, '^IFNA', negate=TRUE))

amp_df <- amp_tib |>
  select(matches('^[A-Z]', ignore.case=FALSE)) |>
  as.matrix()
rownames(amp_df) <- amp_tib$gene

#+ fig.width=8, fig.height=3
pa <- pheatmap(
  amp_df,
  color=c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'),
  border_color=NA,
  breaks=c(0, 0.8, 1.6, 2.4, 3.2, 4),
  cluster_rows=FALSE, cluster_cols=FALSE,
  fontsize=7,
  gaps_col=c(24,36,47))


# for replicability purposes
sessionInfo()
