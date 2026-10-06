#!/usr/bin/env Rscript

"> 05_create_clin_characteristics_table.R <

Does what it says on the tin, with `furniture`.
" -> doc

suppressPackageStartupMessages({
  library(data.table)
  library(furniture)
  library(this.path)
})


setwd(this.path::here())
CLINDATA_FILE <- '../00_common/emseq-rnaseq_clin_details.240716.tsv'

# subset clinical data for these included samples
clin_dt <- fread(CLINDATA_FILE, sep='\t', header=TRUE)
clin_dt <- clin_dt[order(short_id)]
clin_dt[classification == 'Normal squamous', classification := 'NSq']
clin_dt[classification == 'IM', classification := 'NDBE']
clin_dt[classification == 'Cancer', classification := 'EAC']

# subselect samples that were used in the analysis (86 samples)
clin_dt <- clin_dt[emseq_include == TRUE]
stopifnot(nrow(clin_dt) == 86)  # sanity check to confirm magic 86 figure

# subselect columns going into clinical characteristics table
clin_dt <- clin_dt[, c('classification', 'risk', 'gender', 'age_at_collection')]
clin_dt$classification <- factor(clin_dt$classification, levels=c('NSq', 'NDBE', 'LGD', 'HGD', 'EAC'))
clin_dt$risk <- factor(clin_dt$risk, levels=c('L', 'H'))
clin_dt$gender <- factor(clin_dt$gender, levels=c('M', 'F'))

# okay, generate table 1. the majority of the comparisons in this manuscript
# is across classifications, hence this table aims to provide descriptive
# statistics of risk/age/sex across classification
table1(clin_dt,
       risk, gender, age_at_collection,
       splitby='classification',
       test=TRUE,
       var_names=c('Risk', 'Sex', 'Age'),
       export='table1')


# for replicability purposes
sessionInfo()
