#!/usr/bin/env python3

docstring = """
> 05a_calc_focal_cnv_mean_beta.py <

Pre-calculate the mean betas of consistent focal CNVs, so that it can be
overlaid onto the plot created by the 05b script.

For each BLHC sample: calculate mean beta across the entire gene (no N
samples because they're the baseline for CNV coverage calculations, so they're
assumed to be, by default, to not have any CNV).

For every disease state (NBLHC): calculate a mean-of-means across all
samples within that state, so that betas of each sample can be compared to the
mean of each state.

Why? This helps visualise whether CNV deletions/amplifications alter the
methylation of the gene of interest.
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd

BETA_FILE = Path(__file__).resolve(strict=True).parent
BETA_FILE /= '../07_filter_samples/all.beta.filt.samp_excl.tsv.gz'
assert BETA_FILE.is_file(), 'BETA_FILE does not exist!'

# read table with consistent focal deletions
cnv_coords = {}
tsv_reader = csv.reader(open('./raw_supp_tableX.tsv'), delimiter='\t')
next(tsv_reader)    # skip header
for line in tsv_reader:
    # NOTE: coords are 1-based, subsequent code logic assumes 1-base
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[4]
    
    # store details in a tuple
    cnv_coords[gene] = (chrom, start, end)

# repeat for focal amplifications; repeating code is ugly but i'm lazier
tsv_reader = csv.reader(open('./raw_supp_tableY.tsv'), delimiter='\t')
next(tsv_reader)    # skip header
for line in tsv_reader:
    # NOTE: coords are 1-based, subsequent code logic assumes 1-base
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[4]
    
    # store details in a tuple
    cnv_coords[gene] = (chrom, start, end)


# filter for samples with overall coverage > 12 (CNV in lower coverage samples
# tend to be very noisy)
included_columns = {'chr': 'string', 'start': 'UInt32', 'end': 'UInt32'}
tsv_reader = csv.reader(open('./heatmap_inclusion.tsv'), delimiter='\t')
next(tsv_reader)    # skip header
for line in tsv_reader:
    if not line: continue
    
    short_id = line[0]
    mean_cov = float(line[3])
    
    # filter out low-coverage samples
    if mean_cov <= 12: continue
    included_columns[short_id] = 'Float32'


# read the chonky table into memory. keep an eye on mem usage
df = pd.read_table(BETA_FILE, sep='\t', usecols=included_columns.keys(),
                   dtype=included_columns)
df = df[included_columns.keys()]

# create dict to store calcs from each cnv loci
cnv_stats = pd.DataFrame()
for gene in cnv_coords:
    chrom, start, end = cnv_coords[gene]
    
    # slice beta df. all coords used are 1-based
    cnv_df = df[(df['chr'] == chrom) & (df['start'] >= start) & (df['end'] <= end)]
    
    # do calcs and append into dict
    temp = cnv_df.iloc[:, 4:].mean()
    temp.name = f'{gene} | {len(cnv_df)} CpGs'    # also include # CpGs
    for x in 'NBLHC':
        temp[f'mean{x}'] = temp[temp.index.str.startswith(x)].mean()
    
    cnv_stats = pd.concat([cnv_stats, temp], axis=1)

# save the stats dataframe for downstream script
cnv_stats.to_csv('focal_cnv_mean_betas.tsv', sep='\t', na_rep='NA')
