#!/usr/bin/env Rscript

"> 10_overall_dmr_survey.progression.R <

Find and plot trends across the numerous DMRs identified from the pairwise
comparisons of samples across disease classifications along the EAC disease
spectrum.
" -> doc

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(forcats)
  library(GenomicRanges)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(this.path)
})

setwd(this.path::here())
DMR_FOLDER <- './filtered_dmrs/'

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
comparison_order <- c('NDBE vs. NSq', 'EAC vs. NSq', 'LGD vs. NDBE',
                      'HGD/EAC vs. NDBE', 'HGD vs. LGD', 'EAC vs. HGD')

# function to pretty-print diagnostic messages
diag_message <- function(...) {
  message('[', format(Sys.time(), "%H:%M:%S"), '] ', ...)
}


# focus on the hypermeth DMRs first, in the N B L H C sequence + HCvB + CvN
# (last two to fit what commercial tests are doing)
dmr_files <- c(paste0(DMR_FOLDER, 'allB_vs_allN.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allB_vs_allN.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allL_vs_allB.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allL_vs_allB.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allH_vs_allL.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allH_vs_allL.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allC_vs_allH.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allC_vs_allH.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allHC_vs_allB.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allHC_vs_allB.hypometh.tsv.gz'),
               paste0(DMR_FOLDER, 'allC_vs_allN.hypermeth.tsv.gz'),
               paste0(DMR_FOLDER, 'allC_vs_allN.hypometh.tsv.gz'))
stopifnot(file.exists(dmr_files))

# create a dt to store plottable information
plot_tib <- tibble(
  seqnames=character(), start=numeric(), end=numeric(), width=numeric(),
  qval=numeric(), delta_beta=numeric(), nCpG=numeric(), mean_g1=numeric(),
  mean_g2=numeric(), gene_name_prot=character(),
  CpGisland=numeric(), CpGshores=numeric(), nonCpG=numeric(),
  comparison=character(), hyper_or_hypo=character())

for (dmr_file in dmr_files) {
  dmr_tib <- read_tsv(dmr_file, show_col_types=FALSE)
  diag_message(dmr_file, ' has ', formatC(nrow(dmr_tib), big.mark=','), ' rows.')
  
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
  
  # rbind important columns into `plot_tib`
  plot_tib <- bind_rows(plot_tib, dmr_tib[colnames(plot_tib)])
}

# visual check that tibble's fine
head(plot_tib)

# first plot: quick survey of how many differentially methylated regions (DMRs)
# and differentially methylated genes (DMGs) there are, for each comparison
n_tib <- plot_tib |>
  group_by(comparison, hyper_or_hypo) |>
  count(name='count') |>
  arrange(match(comparison, comparison_order)) |>
  mutate(type='DMR')
n_tib

n2_tib <- plot_tib |>
  group_by(comparison, hyper_or_hypo) |>
  summarize(count=n_distinct(gene_name_prot)) |>
  arrange(match(comparison, comparison_order)) |>
  mutate(type='DMG')
n2_tib

n_tib <- bind_rows(n_tib, n2_tib)

maxn_2pc <- round(max(n_tib$count) * 0.02)
g1 <- ggplot(n_tib |> filter(hyper_or_hypo == 'hypermeth'), aes(x=count, y=fct_rev(fct_inorder(comparison)), fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=scales::comma(count), x=count+maxn_2pc), position=position_dodge(width=0.9), hjust=0) +
  scale_x_continuous(labels=scales::comma, limits=c(0, 65 * maxn_2pc)) +
  scale_fill_manual(values=c('#666666', '#aaaaaa')) +
  labs(title='# hypermethylated DMRs') +
  theme_void(12) +
  theme(legend.position='top',
        axis.text.y=element_text(hjust=0.5, face='bold'),
        plot.title=element_text(hjust=0.5))
g2 <- ggplot(n_tib |> filter(hyper_or_hypo == 'hypometh'), aes(x=count, y=fct_rev(fct_inorder(comparison)), fill=type)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=scales::comma(count), x=count+maxn_2pc), position=position_dodge(width=0.9), hjust=1) +
  scale_x_reverse(labels=scales::comma, limits=c(0, 65 * maxn_2pc)) +
  scale_fill_manual(values=c('#666666', '#aaaaaa')) +
  labs(title='# hypomethylated DMRs') +
  theme_void(12) +
  theme(legend.position='top',
        axis.text.y=element_text(hjust=0.5, face='bold'),
        plot.title=element_text(hjust=0.5))

#+ fig.width=7, fig.height=3
plot_grid(g2, g1, rel_widths=c(0.5, 0.5))
ggsave('raw_fig2b.pdf', width=7, height=3)


# second plot: check to see whether the DMRs tend to be in CpG islands
island_tib <- plot_tib |>
  select(CpGisland, CpGshores, nonCpG, comparison, hyper_or_hypo) |>
  group_by(comparison, hyper_or_hypo) |>
  summarize(CpGisland=sum(CpGisland)/sum(CpGisland + CpGshores + nonCpG),
            CpGshores=sum(CpGshores)/sum(CpGisland + CpGshores + nonCpG),
            nonCpG=sum(nonCpG)/sum(CpGisland + CpGshores + nonCpG), .groups='keep')
# add in universal proportions of CpG island / shores / neither to provide
# context for these numbers
island_tib <- island_tib |> as_tibble(.) |>
  add_row(comparison='universe', hyper_or_hypo='hypermeth',
          CpGisland=1789643/25348703, CpGshores=1864104/25348703, nonCpG=21694956/25348703) |>
  add_row(comparison='universe', hyper_or_hypo='hypometh',
          CpGisland=1789643/25348703, CpGshores=1864104/25348703, nonCpG=21694956/25348703) |>
  arrange(match(comparison, c('universe', comparison_order))) |>
  pivot_longer(!c(comparison, hyper_or_hypo), names_to='context', values_to='pct') |>
  mutate(pct=pct*100)
island_tib

g3 <- ggplot(island_tib |> filter(hyper_or_hypo == 'hypermeth'), aes(x=pct, y=fct_rev(fct_inorder(comparison)), fill=fct_rev(fct_inorder(context)))) +
  geom_bar(stat='identity') +
  geom_text(aes(label=paste0(round(pct, 0), '%')), position=position_stack(vjust=0.5), size=3) +
  scale_fill_manual(values=c('#a6cee3', '#fe9f38', '#33a02c')) +
  labs(title='% hypermethylated DMRs in CpG islands', fill='context') +
  theme_void(12) +
  theme(legend.position='top',
        axis.text.y=element_text(hjust=0.5, face='bold'),
        plot.title=element_text(hjust=0.5))

g4 <- ggplot(island_tib |> filter(hyper_or_hypo == 'hypometh'), aes(x=pct, y=fct_rev(fct_inorder(comparison)), fill=fct_rev(fct_inorder(context)))) +
  geom_bar(stat='identity') +
  geom_text(aes(label=paste0(round(pct, 0), '%')), position=position_stack(vjust=0.5), size=3) +
  scale_fill_manual(values=c('#a6cee3', '#fe9f38', '#33a02c')) +
  scale_x_reverse() +
  labs(title='% hypomethylated DMRs in CpG islands', fill='context') +
  theme_void(12) +
  theme(legend.position='top',
        axis.text.y=element_text(hjust=0.5, face='bold'),
        plot.title=element_text(hjust=0.5))

#+ fig.width=7, fig.height=2.5
plot_grid(g4, g3, rel_widths=c(0.5, 0.5))
ggsave('raw_fig2d.pdf', width=7, height=2.5)


# third plot: plot delta beta for top 100 DMRs for each comparison and hyper/hypo
# abs() is needed to get the "most negative" hypometh delta betas
top100_tib <- plot_tib |> 
  group_by(comparison, hyper_or_hypo) |> 
  slice_max(abs(delta_beta), n=100)

# to plot a dumbbell plot, only the values of the rank #1 and #100 are required
top100_tib <- bind_rows(
  top100_tib |> group_by(comparison, hyper_or_hypo) |> slice_max(abs(delta_beta), n=1) |> mutate(rank='rank 1'),
  top100_tib |> group_by(comparison, hyper_or_hypo) |> slice_min(abs(delta_beta), n=1) |> mutate(rank='rank 100'))

# arrange comparison column by disease progression (custom sort)
top100_tib <- top100_tib |> arrange(match(comparison, comparison_order))
top100_tib

g5 <- ggplot(top100_tib |> filter(delta_beta > 0), aes(x=delta_beta, y=fct_rev(fct_inorder(comparison)))) +
  geom_line() +
  geom_point(aes(color=rank), size=3) +
  scale_color_manual(values=c('#d8b365', '#5ab4ac')) +
  xlim(0, 0.8) +
  labs(x='delta \u03B2', y=NULL, color='hypermethylated DMR rank') +
  theme_minimal(12) +
  theme(legend.position='top',
        axis.text.y=element_text(hjust=0.5, face='bold'),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank())

g6 <- ggplot(top100_tib |> filter(delta_beta < 0), aes(x=delta_beta, y=fct_rev(fct_inorder(comparison)))) +
  geom_line() +
  geom_point(aes(color=rank), size=3) +
  scale_color_manual(values=c('#d8b365', '#5ab4ac')) +
  xlim(-0.8, 0) +
  labs(x='delta \u03B2', y=NULL, color='hypomethylated DMR rank') +
  theme_minimal(12) +
  theme(legend.position='top',
        axis.text.y=element_text(hjust=0.5, face='bold'),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank())

#+ fig.width=6, fig.height=2.5
plot_grid(g6, g5, rel_widths=c(0.5, 0.5))
ggsave('raw_fig2c.pdf', width=6.8, height=2.5)


# for replicability purposes
sessionInfo()
