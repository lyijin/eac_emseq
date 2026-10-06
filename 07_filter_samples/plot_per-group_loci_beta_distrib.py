#!/usr/bin/env python3

docstring = """
> calc_per-group_loci_mean_beta.py <

When given chr/start/end coords, plot the per-group beta distributions across
the entire loci.

Expects input coords to conform to the pattern "chrX:123,456..789,101". Colons
and double dots are crucial; thousand separators are removed automatically
by the script.

Uses the hardcoded output file produced by `filter_CpGs_by_coverage.py`,
"all.beta.filt.samp_excl.tsv.gz".
""".strip()

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

BETA_FILE = Path(__file__).resolve(strict=True).parent / 'all.beta.filt.samp_excl.tsv.gz'
VALID_CHROM_FILE = Path(__file__).resolve(strict=True).parent / '../data/grch38p13_lambda_puc/grch38p13_lambda_puc.seqlen.tsv'
DISEASE_STATES = ['NSq', 'NDBE', 'LGD', 'HGD', 'EAC']
DISEASE_TO_COLOR = {'NSq': '#984ea3', 'NDBE': '#377eb8', 'LGD': '#4daf4a',
                    'HGD': '#a65628', 'EAC': '#e41a1c'}
SHORT_TO_LONG_DISEASE = {'N': 'NSq', 'B': 'NDBE', 'L': 'LGD', 'H': 'HGD', 'C': 'EAC'}

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('coord_string', metavar='coord_string', type=str,
                    help='Coords for loci, in the pattern of "chrX:123,456..789,101".')
parser.add_argument('--gene', '-g', metavar='gene_name', type=str,
                    help='Optional gene name for graph title')
args = parser.parse_args()

# sanity check `coord_string`--reading the file and erroring out later is
# computationally expensive
tsv_reader = csv.reader(open(VALID_CHROM_FILE, 'r'), delimiter='\t')
chrom_lengths = {}
for row in tsv_reader:
    chrom_lengths[row[0]] = int(row[1])

scaf = args.coord_string.split(':')[0]
assert scaf in chrom_lengths, \
    f'{scaf} not a valid scaffold name. Valid chroms are {list(chrom_lengths.keys())}'

coords = args.coord_string.split(':')[1]

start_pos = int(coords.split('..')[0].replace('.', '').replace(',', ''))
assert 0 < start_pos <= chrom_lengths[scaf], 'start position {start_pos} not valid.'

end_pos = int(coords.split('..')[1].replace('.', '').replace(',', ''))
assert 0 < end_pos <= chrom_lengths[scaf], 'end position {end_pos} not valid.'

assert start_pos <= end_pos, f'{start_pos} cannot be larger than {end_pos}.'

# table is large, only read specific rows
# store relevant rows as list-of-tuples, then feed into pd.DataFrame.from_records
valid_rows = []

tsv_reader = csv.reader(gzip.open(BETA_FILE, 'rt'), delimiter='\t')
header = next(tsv_reader)

# to speed things up, assume valid rows are contiguous. exit the loop when
# the appending stops
rows_appended_bool = False
for row in tsv_reader:
    if row[0] != scaf: continue
    
    if start_pos > int(row[1]):
        if rows_appended_bool: break
        continue
    if end_pos < int(row[2]):
        if rows_appended_bool: break
        continue
    
    valid_rows.append(row)
    rows_appended_bool = True

assert valid_rows, 'No CpGs with beta values in the query window, script aborted.'
df = pd.DataFrame.from_records(valid_rows, columns=header)

# lop off chr/start/end columns, then calculate per-sample across-loci mean beta
df.drop(df.iloc[:, 0:3], inplace=True, axis=1)
df = df.replace('NA', np.nan).astype(float)    # there should only be beta values

mean_df = pd.DataFrame({'Beta': df.mean(axis=0)})
mean_df['Classification'] = [SHORT_TO_LONG_DISEASE[x[0]] for x in mean_df.index]

# use seaborn to plot violin plot
sns.set_theme(rc={'figure.figsize':(4, 6)})
sns.set_style('whitegrid')

fig = sns.boxplot(mean_df, x='Classification', y='Beta', hue='Classification',
                  order=DISEASE_STATES, palette=DISEASE_TO_COLOR,
                  width=0.5, legend=False, showcaps=False)
fig.set(ylim=(0, 1))
for patch in fig.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .8))
sns.despine(left=True)

if args.gene:
    fig.set_title(f'{args.gene} ({args.coord_string})', loc='left', weight='bold')
    fig.get_figure().savefig(f'{args.gene}.pdf')
else:
    fig.set_title(f'{args.coord_string}', loc='left', weight='bold')
    fig.get_figure().savefig('out.pdf')


