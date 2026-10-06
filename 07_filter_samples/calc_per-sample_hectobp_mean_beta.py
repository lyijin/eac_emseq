#!/usr/bin/env python3

docstring = """
> calc_per-sample_hectobp_mean_beta.py <

When given a tab-separated compiled table of beta values with the columns
    chr \t start \t end \t sample1 \t sample2 \t ..., calculate the beta means
in 100 bp windows (hectobasepair). Expects coords to be 1-based.

Crushing positions into hectobasepair windows is to facilitate downstream
translation on qPCR/dPCR. Amplicon design demands having sufficient CpG
dinucleotides in the window, so there is really no point in calculating mean
betas for windows with too few CpGs.

Assuming we have a diagnostic that gives a positive readout for samples with
amplicon beta >= 0.3 and negative for beta < 0.3, calculates fraction of samples
with positive readouts on a per-disease classification (N B L H C) and per-risk
status (all, high-risk adjacent, low-risk). Would be nice to e.g., create a test
that is low in N B L; then high in H C, while analysing B vs. C results.

Program logic goes:
    1. Read lines belonging to the same window. e.g., chr1 pos12345, 
       chr1 pos12348 and chr1 pos12399 are in the same window
    2. Fewer than 5 CpGs in the window? Discard window, not worth considering.
    3. Pass the min CpG # criteria? Great. Calculate mean beta for the window.
    4. Also calculate fraction of samples with positive readouts on a
       per-disease classification, per-status basis

Assumptions that I made to speed up computation:
    1. Most CpGs are present in the table, avoiding the need to parse the
       actual genome, segmenting it into 100 bp windows, then tallying # CpGs
       in each window
    2. The beta table is a post-processed one where positions with low coverages
       have been removed, so the beta values are "trustable"
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

def calculate_mean_beta(prev_lines):
    """
    Expects a list of lists, and calculates per-sample mean beta in that window.
    
    Returns list of mean betas.
    """
    # calculate means (of 4th column till the end)
    mean_betas = []
    
    for n in range(3, len(prev_lines[0])):
        n_betas = [float(x[n]) for x in prev_lines if x[n] != 'NA']
        if n_betas:
            mean_betas.append(f'{statistics.mean(n_betas):.4f}')
        else:
            # when n_betas is an empty [], statistics.mean([]) errors
            # (i.e., sample has all NA values for that window)
            mean_betas.append('NA')
    
    return mean_betas

def calculate_frac_positive(mean_betas, meaningful_groups):
    """
    For each group in `meaningful_groups`, slice `mean_betas` appropriately
    then calculate the fraction of samples with beta >= 0.3.
    
    Note that numbering in `meaningful_groups` starts at 3 and column indices of
    `mean_betas` starts at 0, hence the offset of 3.
    """
    frac_positive = []
    for mg in meaningful_groups:
        # apply leftward offset of 3 for all values
        bool_betas = [float(mean_betas[x - 3]) >= 0.3 for x in meaningful_groups[mg] \
                      if mean_betas[x - 3] != 'NA']
        
        # guard against empty `bool_betas` due to NAs everywhere
        if len(bool_betas):
            frac_positive.append(f'{sum(bool_betas) / len(bool_betas):.4f}')
        else:
            frac_positive.append('NA')
     
    return frac_positive


# parse clinical details file
clin_df = pd.read_table(CLINICAL_DETAILS_FILE)
clin_df = clin_df[clin_df['emseq_include'] == True]  # select non-excluded samples

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


print ('hectobp_id', *header[3:], *[f'frac{x}' for x in meaningful_groups], sep='\t')
prev_lines = []
prev_hectobp_id = ''
# start reading file to identify hectobasepair windows with >= 5 CpGs
for row in tsv_reader:
    # remember that everything in `row` are strings. cast to int/float when needed
    scaf, start_pos, end_pos = row[:3]
    
    # since hectobasepair windows are geared towards qPCR/dPCR translation
    # and we will never design primers targeting chrM / spiked in controls...
    if scaf in ['chrM', 'lambda', 'pUC19']: continue
    
    # each window starts at scaf:xxx00 and ends at scaf:xxx99 to make
    # programming easy
    curr_hectobp_id = f'{scaf}:{int(int(start_pos) / 100)}'
    if curr_hectobp_id != prev_hectobp_id:
        # check whether previous window was worth calculating mean betas
        if len(prev_lines) >= 5:
            mean_betas = calculate_mean_beta(prev_lines)
            
            # also print out # of CpGs in the window, i.e., chr1:147_7 means
            # chr1:14700..14799 has 7 CpGs with meth data
            print (f'{prev_hectobp_id}_{len(prev_lines)}', *mean_betas, 
                   *calculate_frac_positive(mean_betas, meaningful_groups), sep='\t')
        
        # flush old data, reset counters
        prev_lines = []
        prev_hectobp_id = curr_hectobp_id
    
    prev_lines.append(row)

# handle last hectobp_id
if len(prev_lines) >= 5:
    mean_betas = calculate_mean_beta(prev_lines)
    print (f'{prev_hectobp_id}_{len(prev_lines)}', *mean_betas, 
           *calculate_frac_positive(mean_betas, meaningful_groups), sep='\t')
