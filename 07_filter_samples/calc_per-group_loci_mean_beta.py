#!/usr/bin/env python3

docstring = """
> calc_per-group_loci_mean_beta.py <

When given chr/start/end coords, calculate the per-group beta means across
the entire loci.

Expects input file to be one-loci-per-line, chr \t start \t end. Coords are
1-based. Empty lines are fine (and ignored).

Uses the hardcoded output file produced by `calc_per-group_cpg_mean_beta.py`,
"all.beta.filt.samp_excl.per-group.tsv.gz".
""".strip()

import argparse
import csv
from pathlib import Path

import pandas as pd

BETA_FILE = Path(__file__).resolve(strict=True).parent / 'all.beta.filt.samp_excl.per-group.tsv.gz'

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('loci_file', metavar='tsv_file', type=Path,
                    help='Coords for loci, one per line, in tsv format.')
args = parser.parse_args()

# use pandas to read the table
df = pd.read_table(BETA_FILE)

# print the header row out
print (*list(df.columns), sep='\t')

# then read the input table row-by-row
tsv_reader = csv.reader(open(args.loci_file, 'r'), delimiter='\t')
for row in tsv_reader:
    # lines are silently ignored by printing it out again with no modification
    if len(row) < 3:
        print (*row, sep='\t')
        continue
    
    try:
        loci_start = int(row[1])
        loci_end = int(row[2])
    except ValueError:
        print (*row, sep='\t')
        continue
    
    # accommodate annots with _OT or _CTOB in them
    chr_annot = row[0].split('_')[0]
    
    # get rows matching the criteria
    sliced_df = df.loc[(df['chr'] == chr_annot) & \
                       (df['start'] >= loci_start) & \
                       (df['end'] <= loci_end)]
    group_means = list(sliced_df.iloc[:, 3:].mean(axis=0))
    group_means = [f'{x:.5f}' for x in group_means]
    
    print (*row[:3], *group_means, sep='\t')
