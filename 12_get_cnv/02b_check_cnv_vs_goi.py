#!/usr/bin/env python3

docstring = """
> 02b_check_cnv_vs_goi.py <

Uses the cnvkit library to create a composite plot for the genes of interest
involved in NSq-EAC disease progression.

Overlay the plot with risk labeling (do high-risk samples tend to be associated
with CNVs?) and with symbols indicating higher- or lower-than-expected
methylation beta values (do CNVs alter methylation levels?).
"""

import csv
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

import cnvlib


# read the table produced by 02a for samples with adequate coverage
goi_df = pd.read_csv('./goi_mean_betas.tsv', sep='\t', index_col=0)

# N samples do not have CNV, and meanX rows are not plottable too
# for the plotting df, take the difference in the per-sample mean beta in every
# loci vs. the mean-of-means for that state
#   e.g., goi_df B02's CDKN2A is 0.380, meanB is 0.367.
#         plot_df stores B02 as 0.380 - 0.367 = 0.013
plot_df = pd.concat(goi_df.filter(regex=f'^{x}', axis='rows') - goi_df.loc[f'mean{x}'] for x in 'BLHC')
short_to_long_goi = {x.split('\n')[0]:x for x in list(plot_df.columns)}

# add the filenames into an array
included_samples = []
for short_id in plot_df.index:
    temp = Path('01_run_cnvkit') / f'{short_id}.cns'
    assert temp.is_file(), f'{short_id}.cns file does not exist!'
    included_samples.append(temp)

# read risk labels from the clinical details table. use short_id column as index
clin_df = pd.read_csv('../00_common/emseq-rnaseq_clin_details.240716.tsv',
                      sep='\t', index_col='short_id')
# only interested in risk labeling for this plot
clin_df = clin_df[clin_df.index.isin(plot_df.index)]['risk']

# merge risk info into plot_df
plot_df = pd.concat([plot_df, clin_df], axis='columns')

# read gene-of-interest table
goi_coords = {}
tsv_reader = csv.reader(open('./genes_of_interest.tsv'), delimiter='\t')
for line in tsv_reader:
    goi = line[0]
    long_goi = short_to_long_goi[goi]
    coords = line[2]
    goi_coords[long_goi] = coords

# use cnvlib to read the .cns files and plot stuff out
segments = [cnvlib.read(f) for f in included_samples]
plt.rcParams.update({'font.sans-serif': ['Arial'], 'font.size': 6})
fig, axs = plt.subplots(figsize=(8, 6), ncols=len(goi_coords), nrows=1)
for n, long_goi in enumerate(goi_coords):
    cnvlib.do_heatmap(segments,
                      show_range=goi_coords[long_goi],
                      do_desaturate=True,
                      title=long_goi,
                      ax=axs[n])
    
    # remove colormap--very hacky!
    axs[n]._colorbars[0].remove()
    axs[n]._colorbars = []
    
    # overlay mean beta as text onto plot
    xmin, xmax = axs[n].get_xlim()
    # u2191 is upwards arrow; u2193 is downwards arrow, u2194 is left-right
    for m, delta in enumerate(plot_df[long_goi]):
        # change delta to a unicode symbol
        if delta > 0.02:
            delta_symbol = ' \u2191'
            delta_color = '#4daf4a'
        elif delta < -0.02:
            delta_symbol = ' \u2193'
            delta_color = '#e41a1c'
        else:
            delta_symbol = ' \u2194'
            delta_color = '#333333'
        
        axs[n].text(xmax, m + 0.5, delta_symbol, fontsize=7, fontweight='bold',
                    ha='left', va='center', color=delta_color)
    
    # hide x ticks and x labels
    axs[n].axes.get_xaxis().set_visible(False)
    
    #  depending on which subplot...
    if n == 0:
         # ... for the first one, hide ticks, and add risk labeling
         axs[n].tick_params(tick1On=False)
         x_ticks = axs[0].get_xticks()
         y_ticks = axs[0].get_yticks()
         for m, risk in enumerate(plot_df['risk']):
             risk_color = '#e9a3c9' if risk == 'H' else '#a1d76a'
             rect = Rectangle((x_ticks[0] - 0.012, y_ticks[m] - 0.5), 0.004, 1,
                              color=risk_color, clip_on=False)
             axs[0].add_patch(rect)
                       
         
    elif n > 0:
        # ... for all others, hide y ticks and y labels
        axs[n].axes.get_yaxis().set_visible(False)
    
fig.tight_layout()
fig.savefig('raw_fig3a.pdf')
