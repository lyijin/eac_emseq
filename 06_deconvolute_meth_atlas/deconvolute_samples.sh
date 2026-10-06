#!/bin/bash

# > deconvolute_samples.sh <
#
# Runs https://github.com/nloyfer/meth_atlas on the post-CpG-filtered (but
# pre-sample-filtered) dataset, to eyeball and exclude samples with weird
# methylation patterns that aren't expected from the staging of the disease.

### PRE-PROCESS FILES FOR METH_ATLAS ###
# script was written to filter for HM450 CpGs from the whole-genome EM-seq data;
# and to produce an input file that can be fed into the `meth_atlas`'s
# "deconvolve.py" script
python3 filter_array-specific_CpGs.py HM450.hg38.manifest.tsv.gz ../04_filter_cpgs/all.beta.filt.tsv.gz -v
zcat all.beta.filt.HM450.csv.gz | wc -l  # 459943, makes sense as it's roughly 450k
# produces "all.beta.filt.HM450.csv.gz", _provided_ in this repo

### DOWNLOAD EXTERNAL FILES; NOT PROVIDED IN THIS REPO ###
# get wanding zhou's processed HM450 manifest files (hg38)
wget https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz

# get the `nloyfer/meth_atlas` repo from github
git clone https://github.com/nloyfer/meth_atlas
cd meth_atlas
# test files (provided by `meth_atlas` should work)
python3 deconvolve.py examples.csv --atlas_path reference_atlas.csv
# overwrite the original "deconvolve.py" with the one in this repo. i made
# minor modifications to how samples were ordered on the x-axis, from pure
# alphabetical to one incorporating the increasing N/B/L/H/C disease severity
mv ../deconvolve.py .
# then point it to the file produced by the pre-processing step
python3 deconvolve.py ../all.beta.filt.HM450.csv.gz --atlas_path reference_atlas.csv

# assuming things went well, there will be two output files in `meth_atlas/`;
# move it to the folder one above
mv all.beta.filt_deconv_output.csv ../all.beta.filt.HM450_deconv_output.csv
mv all.beta.filt_deconv_plot.png ../all.beta.filt.HM450_deconv_plot.png
