#!/usr/bin/env python3

docstring = """
> calc_per-group_cpg_mean_beta.py <

Based on the N B L H C disease continuum, calculate per-group mean betas on a 
per-position basis.

Also incorporates risk stratification--some disease states e.g., L, is further
subdivided into high-risk and low-risk. Some disease states (C), however, only
has high-risk. Which is perhaps why pooling all low vs. high might not be
meaningful as the differences might be driven by high-risk having more late-
stage samples than low-risk.

"High-risk" are mostly classed as "adj" (adjacent to tissue with more serious
disease pathologies) than "prog" (samples from patients that progressed later
to more serious pathologies). As there were only three samples in the latter,
they're excluded from "high-risk" averaging. i.e., high-risk refers exclusively
to "adj", which potentially indicates sampling error.

Have defaulted to calculating per-group, per-risk stratification i.e.,
  N_mean      mean of all N samples
  N_L_mean    mean of all N & low-risk samples
  N_H_mean    mean of all N & high-risk samples
  B_mean      mean of all B samples
  [...]

Risk info is read from the clinical details table.
""".strip()

import csv
import gzip
from pathlib import Path
import statistics

import pandas as pd

BETA_FILE = Path('all.beta.filt.samp_excl.tsv.gz')
CLINICAL_DETAILS_FILE = Path('../00_common/emseq-rnaseq_clin_details.240716.tsv')
DISEASE_STATES = ['N', 'B', 'L', 'H', 'C']
RISK_STATES = ['L', 'H']

# parse clinical details file for risk details
clin_df = pd.read_table(CLINICAL_DETAILS_FILE)
clin_df = clin_df[clin_df['emseq_sample_no'].notna()]  # select EM-seq samples
# sample exclusions do not need to be handled here, the BETA_FILE table would
# have those samples excluded already. but take note that `clin_df` might
# contain more sample IDs than present in the header of BETA_FILE

# handle potential compression
if BETA_FILE.suffix == '.gz':
    # assume gzip-compressed file
    tsv_reader = csv.reader(gzip.open(BETA_FILE, 'rt'), delimiter='\t')
else:
    # assume it's plaintext
    tsv_reader = csv.reader(open(BETA_FILE, 'r'), delimiter='\t')

# header has short IDs, cross-ref `clin_df` to get disease/risk states
header = next(tsv_reader)

# meaningful groups are (for now)
#   a) N_mean B_mean L_mean H_mean C_mean
#   b) LR_mean HR_mean
#   c) N_LR_mean B_LR_mean L_LR_mean H_LR_mean
#      N_HR_mean B_HR_mean L_HR_mean H_HR_mean
# this dict contains arrays that contain COLUMN INDICES (integers), to make
# slicing easier later on
meaningful_groups = {}

for ds in DISEASE_STATES:  # (a)
    sample_ids = clin_df[clin_df['short_id'].str.startswith(ds)]['short_id'].values
    
    header_indices = [header.index(x) for x in sample_ids if x in header]
    meaningful_groups[ds] = sorted(header_indices)

for rs in RISK_STATES:  # (b)
    # "high-risk" == "high-risk" && "adj"
    risk_df = clin_df[clin_df['risk'] == rs]
    if rs == 'H':
        risk_df = risk_df[risk_df['H_risk_category'] == 'adj']
    sample_ids = risk_df['short_id'].values
    
    header_indices = [header.index(x) for x in sample_ids if x in header]
    meaningful_groups[f'{rs}R'] = sorted(header_indices)

for rs in RISK_STATES:  # (c)
    # "high-risk" == "high-risk" && "adj"
    risk_df = clin_df[clin_df['risk'] == rs]
    if rs == 'H':
        risk_df = risk_df[risk_df['H_risk_category'] == 'adj']
    
    # cancer is excluded here as there are no high-risk cancer adjacent to
    # anything more serious, as cancer is THE most serious state
    for ds in DISEASE_STATES[:-1]:
        sample_ids = risk_df[risk_df['short_id'].str.startswith(ds)]['short_id'].values
    
        header_indices = [header.index(x) for x in sample_ids if x in header]
        meaningful_groups[f'{ds}_{rs}R'] = sorted(header_indices)

# start reading file and printing group means out
print ('chr', 'start', 'end', *[f'{x}_mean' for x in meaningful_groups], sep='\t')
for row in tsv_reader:
    # print per-disease state mean beta
    mean_betas = []
    for mg in meaningful_groups:
        sliced_betas = [float(row[x]) for x in meaningful_groups[mg] if row[x] != 'NA']
        if sliced_betas:
            mean_betas.append(f'{statistics.mean(sliced_betas):.5f}')
        else:
            mean_betas.append('NA')
    
    # keep chr/start/end from original row
    print (*row[:3], *mean_betas, sep='\t')
