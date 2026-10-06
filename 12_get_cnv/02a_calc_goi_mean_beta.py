#!/usr/bin/env python3

docstring = """
> 02a_calc_goi_mean_beta.py <

Pre-calculate the mean betas of the genes of interest, so that it can be
overlaid onto the plot created by the 02b script.

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

# read gene-of-interest table
goi_coords = {}
tsv_reader = csv.reader(open('./genes_of_interest.tsv'), delimiter='\t')
for line in tsv_reader:
    # provide more context to the genes by including stage at which
    # amp/del is expected (+/- respectively), length of amplicon and # CpGs
    coords = line[2]
    chrom, startend = coords.split(':')
    start, end = [int(x) for x in startend.split('-')]
    goi_len = end - start + 1
    
    goi = f'{line[0]}\n({line[1]})\n{goi_len} bp'
    goi_coords[goi] = coords

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

# create dict to store calcs from each goi loci
goi_stats = pd.DataFrame()
for goi in goi_coords:
    # convert goi string e.g., 'chr9:21967752-21995324' into chrom/start/end
    chrom, startend = goi_coords[goi].split(':')
    start, end = [int(x) for x in startend.split('-')]
    
    # slice beta df. all coords used are 1-based
    goi_df = df[(df['chr'] == chrom) & (df['start'] >= start) & (df['end'] <= end)]
    
    # do calcs and append into dict
    temp = goi_df.iloc[:, 4:].mean()
    temp.name = f'{goi}\n{len(goi_df)} CpGs'    # also include # CpGs
    for x in 'NBLHC':
        temp[f'mean{x}'] = temp[temp.index.str.startswith(x)].mean()
    
    goi_stats = pd.concat([goi_stats, temp], axis=1)

# save the stats dataframe for downstream script
goi_stats.to_csv('goi_mean_betas.tsv', sep='\t', na_rep='NA')
