#!/usr/bin/env Rscript

"> 01_plot_coverages.R <

Compares the post-filtered (not raw!) coverages of the 96 datasets.

Aims:

  1. check whether filtering was overly aggressive.

  2. plot coverages vs. clinical variables like gender, age, and risk class to
     check (and exclude) coverage as a potential reason driving downstream beta
     differences in certain subgroups.
" -> doc

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

# make sure sample number order is identical, as subsequent lines assume this
# to be true
stopifnot(clin_dt$short_id == colnames(cov_dt))

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

# also use proper stats to check
summary(lm(mean_cov ~ risk, mean_cov_dt))
summary(lm(mean_cov ~ age, mean_cov_dt))
summary(lm(mean_cov ~ gender, mean_cov_dt))

# write `mean_cov_dt` to file, so that it can be used in downstream scripts
fwrite(mean_cov_dt, file='mean_cov_dt.tsv', sep='\t')

# for replicability purposes
sessionInfo()
