#!/bin/bash

# > filter_samples.sh <
#
# Post-deciding the removal of 10 samples due to uncertain/incorrect
# classification, implement the removal as Python/R scripts parked in
# `../00_common/` so that it can be used across the project. Also include
# use of scripts to facilitate downstream PCR amplicon design.

### REMOVE MISCLASSFIED SAMPLES ###
python3 ../00_common/exclude_misclassified_samples.py -v ../04_filter_cpgs/all.beta.filt.tsv.gz > tmp

# and properly sort chr/pos in ascending order 
sort -k1,1 -k2,2n tmp | gzip > all.beta.filt.samp_excl.tsv.gz
rm -f tmp
###


### TO ENABLE HIGH-THROUGHPUT AMPLICON DESIGNS ###
# pick out 100 bp windows ("hectobp" windows) containing sufficient # of CpGs
# to facilitate easier amplicon design--i.e., can't design methylation-specific
# amplicons in genomic regions with too few CpGs
python3 calc_per-sample_hectobp_mean_beta.py | gzip > all.beta.filt.samp_excl.hectobp.tsv.gz
# note: this script has hardcoded input files (read the *.py script)

# calculate per-group mean CpGs, dump into plaintext file for next script
# ("calc_per-group_loci_mean_beta.py") to calculate per-loci, per-group mean CpGs
# e.g., to figure out whether existing commercial primers are able to properly
# distinguish OAC vs NDBO based on their amplified region on our data
python3 calc_per-group_cpg_mean_beta.py | gzip > all.beta.filt.samp_excl.per-group.tsv.gz
###
