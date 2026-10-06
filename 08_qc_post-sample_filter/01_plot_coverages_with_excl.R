#!/usr/bin/env Rscript

"> 01_plot_coverages_with_excl.R <

Compares the post-filtered (not raw!) coverages of the 96 datasets.

Aims:

  1. check whether filtering was overly aggressive.

  2. plot coverages vs. clinical variables like gender, age, and risk class to
     check (and exclude) coverage as a potential reason driving downstream beta
     differences in certain subgroups.

After considering exploratory RNA-seq and EM-seq results, some samples look
obviously misclassified. Decision was made to exclude these samples. To reduce
duplication in code, this exclusion was carried out with a Python script
(in `../00_common/`) to produce the `all.beta.filt.samp_excl.tsv.gz` file.
An R script in the same directory handles sample exclusion in R dataframes.
" -> doc

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(this.path)
})


setwd(this.path::here())
COMPILED_COV_FILE <- '../04_filter_cpgs/all.cov.filt.tsv.gz'
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

# read clindata file as data.table, discard NAs in `emseq_sample_no` column,
# then order by `short_id`
clin_dt <- fread(CLINDATA_FILE, sep='\t', colClasses='character')
clin_dt <- clin_dt[!is.na(emseq_sample_no)][order(short_id)]
clin_dt[classification == 'Normal squamous', classification := 'NSq']
clin_dt[classification == 'IM', classification := 'NDBE']
clin_dt[classification == 'Cancer', classification := 'EAC']
rownames(clin_dt) <- clin_dt$short_id
clin_dt[1:10, ]

# read compiled cov data, omit second and third columns (start/end)
cov_dt <- fread(COMPILED_COV_FILE, sep='\t', header=TRUE, drop=2:3)
# remove lambda and pUC19 coverages
cov_dt <- cov_dt[chr %like% '^chr']
# then remove the first column
cov_dt$chr <- NULL

# get per-CpG coverage by adding odd columns with even columns
cov_dt <-
  cov_dt[, .SD, .SDcols=seq(1,ncol(cov_dt),2)] + 
  cov_dt[, .SD, .SDcols=seq(2,ncol(cov_dt),2)]

# sort cov_dt by column names (short IDs)
setcolorder(cov_dt, c(order(names(cov_dt))))
cov_dt[1:5, 1:10]

# drop excluded samples in datatable (detailed in `exclude_misclassified_samples.R`)
source('../00_common/exclude_misclassified_samples.R')
cov_dt <- exclude_samples_by_col(cov_dt)
clin_dt <- exclude_samples_by_row(clin_dt)

# defensive programming: make sure sample number order is identical, as
# subsequent lines assume this to be true
stopifnot('`clin_dt` should have 86 rows' = nrow(clin_dt) == 86)
stopifnot('`cov_dt` should have 86 cols' = ncol(cov_dt) == 86)
stopifnot('ordering of `clin_dt` and `cov_dt` must be identical' =
            clin_dt$short_id == colnames(cov_dt))

# plot overall coverage 
df <- rbind(
  data.frame(x = rowSums(cov_dt >= 5), mincov='\u2265 5'),
  data.frame(x = rowSums(cov_dt >= 10), mincov='\u2265 10'),
  data.frame(x = rowSums(cov_dt >= 20), mincov='\u2265 20'),
  data.frame(x = rowSums(cov_dt >= 30), mincov='\u2265 30'))
df$mincov <- as.factor(df$mincov)


# create a def to deal with the relabeling of the y-axis (flipped to x)
prettify_num <- function (x) {
  # ecdf produces a range of numbers [0-1], but the log10 transformation
  # changes these values to e.g., -4, -3, -2, -1, 0
  actual_num <- 10^x * nrow(df) / 4
  actual_num <- 
    ifelse(actual_num > 1e6,
           paste(round(actual_num/1e6, 1), 'M'), 
           ifelse(actual_num > 1e3,
                  paste(round(actual_num/1e3, 1), 'K'),
                  round(actual_num, 0)))
  actual_num
}

#+ fig.width=4, fig.height=4
# trial-and-error using this as baseline
# https://stackoverflow.com/questions/77846391/log-scale-on-reverse-cumulative-distribution-plot-in-ggplot2
ggplot(df, aes(x, color=fct_inorder(mincov))) +
  stat_ecdf(aes(y=log10(1-after_stat(y))), geom='line', pad=TRUE) +
  scale_y_continuous('# of CpGs', labels=prettify_num) +
  coord_flip() +
  labs(title='',
       x=expression('Samples with # CpGs \u2265'~italic(t)),
       y='# of CpGs',
       color=expression('Coverage\nthresholds'~italic(t))) +
  theme_minimal(12) +
  theme(legend.position='inside',
        legend.position.inside=c(0.05, 0.05),
        legend.justification=c('left', 'bottom'),
        panel.grid.minor=element_blank())
ggsave('raw_fig1b.pdf', width=4, height=4)
rm(df)

# calculate per-sample mean cov, add sample classification and short ID
mean_cov_dt <- data.table(classification=clin_dt$classification,
                          risk=clin_dt$risk,
                          short_id=clin_dt$short_id,
                          gender=clin_dt$gender,
                          age=clin_dt$age_at_collection,
                          mean_cov=colMeans(cov_dt))
mean_cov_dt$classification <- factor(mean_cov_dt$classification, 
                                     levels=disease_progression_vec)
mean_cov_dt$age <- as.numeric(mean_cov_dt$age)
# bankers' rounding to nearest 10
mean_cov_dt$rounded_age <- round(mean_cov_dt$age / 10) * 10

# do coverages differ...

# per-classification? i.e., are cancer samples sequenced at lower depth?
mean_cov_dt[, .(.N, mean(mean_cov)), by=.(classification)]

# by risk? i.e., are high-risk samples sequenced at lower depth?
mean_cov_dt[, .(.N, mean(mean_cov)), by=.(risk)]

# by gender? i.e., are male samples sequenced at lower depth?
mean_cov_dt[, .(.N, mean(mean_cov)), by=.(gender)]

# by age rounded to nearest tens?
mean_cov_dt[, .(.N, mean(mean_cov)), by=.(rounded_age)]

# when aggregated, coverage differences between groups aren't obvious. good

#+ fig.width=6, fig.height=6
# check which datasets are lowly covered, keep those labels in mind as low
# coverages can affect PCA/hierarchical clustering (more "digital" betas)
ggplot(mean_cov_dt, aes(x=classification, y=mean_cov)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=classification), position=position_jitter(seed=1337), size=4, alpha=0.6) +
  geom_text(aes(label=short_id), position=position_jitter(seed=1337), hjust=0.5, vjust=1) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(title='Per-sample mean coverage (for filtered positions)',
       x='Disease classification',
       y='Mean coverage') +
  theme_minimal(12) +
  theme(legend.position='none')

# plot cleaner version for manuscript (& with bigger axes text)
#+ fig.width=4, fig.height=2.5
ggplot(mean_cov_dt, aes(x=classification, y=mean_cov)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=classification), position=position_jitter(seed=1337), size=3, alpha=0.6) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(x='Disease classification',
       y='Mean coverage') +
  theme_minimal(12) +
  theme(legend.position='none')
ggsave('raw_fig1c1.pdf', width=4, height=2.5)

#+ fig.width=6, fig.height=6
# any converage-related trends wrt risk? null expectation is NO trend
ggplot(mean_cov_dt, aes(x=classification, y=mean_cov)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=classification, shape=risk), position=position_jitter(seed=1337), size=4, alpha=0.6) +
  geom_text(aes(label=short_id), position=position_jitter(seed=1337), hjust=0.5, vjust=1) +
  scale_color_manual(values=disease_to_color_vec) +
  labs(title='Per-sample mean coverage (for filtered positions)',
       x='Disease classification',
       y='Mean coverage') +
  theme_minimal(12)

# any converage-related trends wrt gender & age at collection? null expectation
# is, again, no trend
ggplot(mean_cov_dt, aes(x=classification, y=mean_cov)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=age, shape=gender), position=position_jitter(seed=1337), size=4, alpha=0.8) +
  geom_text(aes(label=short_id), position=position_jitter(seed=1337), hjust=0.5, vjust=1) +
  scale_color_distiller(palette='RdYlBu', limits=c(45,85)) +
  labs(title='Per-sample mean coverage (for filtered positions)',
       x='Disease classification',
       y='Mean coverage') +
  theme_minimal(12)

# also use proper stats to check (also, for manuscript)
anova(lm(mean_cov ~ classification, mean_cov_dt))
summary(lm(mean_cov ~ risk, mean_cov_dt))
summary(lm(mean_cov ~ age, mean_cov_dt))
summary(lm(mean_cov ~ gender, mean_cov_dt))

# write `mean_cov_dt` to file, so that it can be used in downstream scripts
fwrite(mean_cov_dt, file='mean_cov_dt.tsv', sep='\t')

# for replicability purposes
sessionInfo()
