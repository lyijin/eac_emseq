#!/bin/bash

# > metilene_analysis.cell_types.sh <
#
# commands to quickly regenerate relevant comparisons by cell-type groupings.
#
# there are four groups--columnar, squamous, lymphoid and erythroid. as
# "squamous" is NSq samples + a few NDBEs + a single HGD, this is probably
# the best baseline for the other three groups to be compared against. the
# three separate pairwise comparisons can then be "united" in a single heatmap
# that tracks whether certain promoters are differentially methylated in one
# other cluster, or multiple other clusters.
#
# implementation-wise, there's no point in carrying out Squamous vs X AND
# X vs Squamous, as both resulting tables are identical with the exception of
# group means swapping columns (one column is for group mean 1, next column is
# for group mean 2).
#
# long story short, an example demostrating nomenclature for comparisons are
#   X_vs_Squamous.dmr.tsv, where
#     X is other three cell types
#     Y is samples in the "Squamous" group ("baseline" state)
# delta beta values reported in the tables are [X - Y],
#   i.e., positive values are hypermethylated biomarkers in [one of the three]
#         than [Squamous]

# create output folder (and ignore errors if it exists)
mkdir -p metilene_output
zcat ../07_filter_samples/all.beta.filt.samp_excl.tsv.gz > all.beta.filt.samp_excl.tsv
rm -f metilene_analysis.log

### custom comparisons by cell-type groupings ###
#
# heatmap groupings are stored in "heatmap_annots.tsv". challenge is to parse
# the tsv for sample names, add an X in front of the non-Squamous group, and
# add a Y for the Squamous group
#
# ... the non-Squamous samples vs. Squamous
#     [varies]                     [B04 B12 B21 B26 B28 B30 H15 N01 N02 N03 N04 N05 N06 N07 N08 N09 N10]
for a in Columnar Lymphoid Erythroid; do
    head -1 all.beta.filt.samp_excl.tsv > header.tmp
    for b in `grep "${a}" heatmap_annots.tsv | cut -f 1`; do
        sed -i "s/${b}/X${b}/" header.tmp
    done

    for c in `grep 'Squamous' heatmap_annots.tsv | cut -f 1`; do
        sed -i "s/${c}/Y${c}/" header.tmp
    done

    cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > non-squamous_vs_squamous.tmp
    ~/tools/metilene/metilene_linux64 -a X -b Y -t 10 non-squamous_vs_squamous.tmp > tmp
    sort -k1,1 -k2,2n tmp | gzip > metilene_output/${a}_vs_Squamous.dmr.tsv.gz
    rm -f tmp
    echo "[$(date)] ${a} vs. Squamous. Samples were" >> metilene_analysis.log
    echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log
done

# repeat for Columnar
# ... Lymphoid/Erythroid vs. Columnar
#     [varies]               [B02 B03 B06 B08 B09 B13 B15 B17 B18 B20 B23 B24 B29 C03 C11 H10 H14 L07 L11 L14 L15 L17]
for a in Lymphoid Erythroid; do
    head -1 all.beta.filt.samp_excl.tsv > header.tmp
    for b in `grep "${a}" heatmap_annots.tsv | cut -f 1`; do
        sed -i "s/${b}/X${b}/" header.tmp
    done

    for c in `grep 'Columnar' heatmap_annots.tsv | cut -f 1`; do
        sed -i "s/${c}/Y${c}/" header.tmp
    done

    cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > non-columnar_vs_columnar.tmp
    ~/tools/metilene/metilene_linux64 -a X -b Y -t 10 non-columnar_vs_columnar.tmp > tmp
    sort -k1,1 -k2,2n tmp | gzip > metilene_output/${a}_vs_Columnar.dmr.tsv.gz
    rm -f tmp
    echo "[$(date)] ${a} vs. Columnar. Samples were" >> metilene_analysis.log
    echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log
done

# and for Lymphoid
# ... Erythroid                                                                vs. Lymphoid
#    [C02 C06 C07 C09 C10 C13 H01 H02 H03 H04 H05 H06 L02 L10 L12 L18 L20 L22]     [B11 B14 B16 B19 B22 B25 B27 B31 B32 B33 C04 C08 C12 C14 H08 H09 H11 H12 H13 H16 H17 L03 L05 L06 L08 L09 L13 L19 L21]
for a in Erythroid; do
    head -1 all.beta.filt.samp_excl.tsv > header.tmp
    for b in `grep "${a}" heatmap_annots.tsv | cut -f 1`; do
        sed -i "s/${b}/X${b}/" header.tmp
    done

    for c in `grep 'Lymphoid' heatmap_annots.tsv | cut -f 1`; do
        sed -i "s/${c}/Y${c}/" header.tmp
    done

    cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > non-lymphoid_vs_lymphoid.tmp
    ~/tools/metilene/metilene_linux64 -a X -b Y -t 10 non-lymphoid_vs_lymphoid.tmp > tmp
    sort -k1,1 -k2,2n tmp | gzip > metilene_output/${a}_vs_Lymphoid.dmr.tsv.gz
    rm -f tmp
    echo "[$(date)] ${a} vs. Lymphoid. Samples were" >> metilene_analysis.log
    echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log
done

# final cleanup
rm -f header.tmp
rm -f non-squamous_vs_squamous.tmp
rm -f non-columnar_vs_columnar.tmp
rm -f non-lymphoid_vs_lymphoid.tmp
rm -f all.beta.filt.samp_excl.tsv
