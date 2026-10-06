#!/bin/bash

# > metilene_pairwise_comparisons.sh <
#
# commands to quickly regenerate biologically relevant comparisons.
#
# there's no point in carrying out N vs B AND B vs N, as both resulting tables
# are identical with the exception of group means swapping columns (one column
# is for group mean 1, next column is for group mean 2).
#
# nomenclature for comparisons are
#   allX_vs_allY.dmr.tsv, where
#     X is more serious state along disease progression scale
#     Y is less serious ("baseline" state)
# delta beta values reported in the tables are [X - Y],
#   i.e., positive values are hypermethylated biomarkers in [more serious state]
#         than [less serious state]

# create output folder (and ignore errors if it exists)
mkdir -p metilene_output
zcat ../07_filter_samples/all.beta.filt.samp_excl.tsv.gz > all.beta.filt.samp_excl.tsv
rm -f metilene_analysis.log

### pairwise comparisons across disease states ###
#
# N vs B/L/H/C
for a in B L H C; do ~/tools/metilene/metilene_linux64 -a ${a} -b N -t 10 all.beta.filt.samp_excl.tsv > tmp && sort -k1,1 -k2,2n tmp | gzip > metilene_output/all${a}_vs_allN.dmr.tsv.gz && rm -f tmp && echo "[$(date)] Completed ${a} vs N." >> metilene_analysis.log; done

# B vs L/H/C
for a in L H C; do ~/tools/metilene/metilene_linux64 -a ${a} -b B -t 10 all.beta.filt.samp_excl.tsv > tmp && sort -k1,1 -k2,2n tmp | gzip > metilene_output/all${a}_vs_allB.dmr.tsv.gz && rm -f tmp && echo "[$(date)] Completed ${a} vs B." >> metilene_analysis.log; done

# L vs H/C
for a in H C; do ~/tools/metilene/metilene_linux64 -a ${a} -b L -t 10 all.beta.filt.samp_excl.tsv > tmp && sort -k1,1 -k2,2n tmp | gzip > metilene_output/all${a}_vs_allL.dmr.tsv.gz && rm -f tmp && echo "[$(date)] Completed ${a} vs L." >> metilene_analysis.log; done

# H vs C
for a in C; do ~/tools/metilene/metilene_linux64 -a ${a} -b H -t 10 all.beta.filt.samp_excl.tsv > tmp && sort -k1,1 -k2,2n tmp | gzip > metilene_output/all${a}_vs_allH.dmr.tsv.gz && rm -f tmp && echo "[$(date)] Completed ${a} vs H." >> metilene_analysis.log; done

# HC grouped together, then vs N/B/L
for a in N B L; do ~/tools/metilene/metilene_linux64 -a C -b ${a} -t 10 <(sed '1 s/H/CC/g' all.beta.filt.samp_excl.tsv) > tmp && sort -k1,1 -k2,2n tmp | gzip > metilene_output/allHC_vs_all${a}.dmr.tsv.gz && rm -f tmp && echo "[$(date)] Completed HC vs ${a}." >> metilene_analysis.log; done

### pairwise comparisons across risk categories (for same disease states) ###
#
# i.e., "high-risk & adjacent" vs. "low-risk"
#
# this is a bit hacky. metilene compares two groups of samples; group membership
# is defined by the first character of the sample name. using "H" or "L" to
# denote high- and low-risk is not possible, as they already mean "HGD" and
# "LGD" in this dataset.
#
# rename "high-risk" samples to start with "X" instead (append to start of
# string); while "low-risk" starts with Y. risk classifications can be found in
# ../00_common/emseq-rnaseq_clin_details.240716.tsv.
# doing things on the command line can be quite hacky, but here goes!
#
# note: N   H vs. L is considered "group 5" in 4/7/24 email
#       B   H vs. L is considered "group 6" in 3/9/24 email
#       BLH H vs. L is considered "group 1" in 4/7/24 email

for a in N B L H BL BLH; do
    # copy header line out into separate file
    head -1 all.beta.filt.samp_excl.tsv > header.tmp
    
    # high-risk samples of that disease state now start with X
    while read b; do
        c=`echo ${b} | sed 's/^/X/'`
        sed -i "s/${b}/${c}/" header.tmp
    done < <(cat ../00_common/emseq-rnaseq_clin_details.240716.tsv | awk -F'\t' '$4 == "H" && $13 == "adj"' | cut -f 8 | grep -P "^[${a}]")
    
    # low-risk samples of that disease state start with Y
    # (reason why this is a separate step is because some high-risk samples are
    # not considered "adj" and thus not renamed to X. these samples will appear
    # "low-risk" even when they aren't, hence actual low-risks have to be
    # renamed to something distinguishable)
    while read d; do
        e=`echo ${d} | sed 's/^/Y/'`
        sed -i "s/${d}/${e}/" header.tmp
    done < <(cat ../00_common/emseq-rnaseq_clin_details.240716.tsv | awk -F'\t' '$4 == "L"' | cut -f 8 | grep -P "^[${a}]")
    
    # replace old header with new, and also replace the "NA"s in the original
    # table into '-'
    cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv | sed 's/NA/-/g') > high_vs_low_risk.tmp
    
    # run metilene on this frankenfile
    ~/tools/metilene/metilene_linux64 -a X -b Y -t 10 high_vs_low_risk.tmp > tmp
    sort -k1,1 -k2,2n tmp | gzip > metilene_output/high${a}_vs_low${a}.dmr.tsv.gz
    rm -f tmp
    echo "[$(date)] Completed high-risk vs. low-risk ${a}. Samples were" >> metilene_analysis.log
    echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log
done

# cleanup, remove all temp files (for now)
rm -f high_vs_low_risk.tmp

### custom comparisons for high-risk samples ###
#
# if sampling error occurred, it's more "OK" to fail to catch errors that
# progressed to L (ablation), than errors than progressed to C (radical
# oesophagectomy).
#
# as these samples have inconsistent annotations in the clin_details sheet,
# ye olde hardcoding stuff to work is probably the simplest way forwards.
#
# enough blah blah. same program logic as above.
#
# are there DMRs that distinguish...
#
# ... high-risk B, worst diagnosis HC vs. high-risk B, worst diagnosis L
#    [B16 B17 B26 B31 B32]               [B02 B11 B20 B23 B25]
head -1 all.beta.filt.samp_excl.tsv > header.tmp
for b in B16 B17 B26 B31 B32; do
    c=`echo ${b} | sed 's/^/X/'`
    sed -i "s/${b}/${c}/" header.tmp
done

for d in B02 B11 B20 B23 B25; do
    e=`echo ${d} | sed 's/^/Y/'`
    sed -i "s/${d}/${e}/" header.tmp
done

cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > higher_vs_lower_risk.tmp
~/tools/metilene/metilene_linux64 -a X -b Y -t 10 higher_vs_lower_risk.tmp > tmp
sort -k1,1 -k2,2n tmp | gzip > metilene_output/highBworstHC_vs_highBworstL.dmr.tsv.gz
rm -f tmp
echo "[$(date)] high-risk B, worst diagnosis HC vs. high-risk B, worst diagnosis L. Samples were" >> metilene_analysis.log
echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log

# ... high-risk L, worst diagnosis C vs. high-risk L, worst diagnosis H
#    [L02 L05 L07 L08 L10 L11 L17 L21]   [L09 L12 L13 L20 L22]
head -1 all.beta.filt.samp_excl.tsv > header.tmp
for b in L02 L05 L07 L08 L10 L11 L17 L21; do
    c=`echo ${b} | sed 's/^/X/'`
    sed -i "s/${b}/${c}/" header.tmp
done

for d in L09 L12 L13 L20 L22; do
    e=`echo ${d} | sed 's/^/Y/'`
    sed -i "s/${d}/${e}/" header.tmp
done

cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > higher_vs_lower_risk.tmp
~/tools/metilene/metilene_linux64 -a X -b Y -t 10 higher_vs_lower_risk.tmp > tmp
sort -k1,1 -k2,2n tmp | gzip > metilene_output/highLworstC_vs_highLworstH.dmr.tsv.gz
rm -f tmp
echo "[$(date)] Completed high-risk L, worst diagnosis C vs. high-risk L, worst diagnosis H. Samples were" >> metilene_analysis.log
echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log

# ... high-risk BL, worst diagnosis C               vs. high-risk BL, worst diagnosis LH
#    [B16 B17 B26 B31 L02 L05 L07 L08 L10 L11 L17 L21]  [B02 B11 B20 B23 B25 B32 L09 L12 L13 L20 L22]
head -1 all.beta.filt.samp_excl.tsv > header.tmp
for b in B16 B17 B26 B31 L02 L05 L07 L08 L10 L11 L17 L21; do
    c=`echo ${b} | sed 's/^/X/'`
    sed -i "s/${b}/${c}/" header.tmp
done

for d in B02 B11 B20 B23 B25 B32 L09 L12 L13 L20 L22; do
    e=`echo ${d} | sed 's/^/Y/'`
    sed -i "s/${d}/${e}/" header.tmp
done

cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > higher_vs_lower_risk.tmp
~/tools/metilene/metilene_linux64 -a X -b Y -t 10 higher_vs_lower_risk.tmp > tmp
sort -k1,1 -k2,2n tmp | gzip > metilene_output/highBLworstC_vs_highBLworstLH.dmr.tsv.gz
rm -f tmp
echo "[$(date)] Completed high-risk BL, worst diagnosis C vs. high-risk BL, worst diagnosis LH. Samples were" >> metilene_analysis.log
echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log

# ... high-risk BL, worst diagnosis HC                                      vs. high-risk BL, worst diagnosis L + low-risk BL. also known as "group 3" in 4/7/24 email
#    [B16 B17 B26 B31 B32 L02 L05 L07 L08 L09 L10 L11 L12 L13 L17 L20 L21 L22]  [B02 B03 B04 B06 B08 B09 B11 B12 B13 B14 B15 B18 B20 B21 B22 B23 B24 B25 B27 B28 B30 B33 L03 L06 L14 L15 L18 L19]
head -1 all.beta.filt.samp_excl.tsv > header.tmp
for b in B16 B17 B26 B31 B32 L02 L05 L07 L08 L09 L10 L11 L12 L13 L17 L20 L21 L22; do
    c=`echo ${b} | sed 's/^/X/'`
    sed -i "s/${b}/${c}/" header.tmp
done

for d in B02 B03 B04 B06 B08 B09 B11 B12 B13 B14 B15 B18 B20 B21 B22 B23 B24 B25 B27 B28 B30 B33 L03 L06 L14 L15 L18 L19; do
    e=`echo ${d} | sed 's/^/Y/'`
    sed -i "s/${d}/${e}/" header.tmp
done

cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > higher_vs_lower_risk.tmp
~/tools/metilene/metilene_linux64 -a X -b Y -t 10 higher_vs_lower_risk.tmp > tmp
sort -k1,1 -k2,2n tmp | gzip > metilene_output/highBLworstHC_vs_highBLworstLlowBL.group3.dmr.tsv.gz
rm -f tmp
echo "[$(date)] Completed high-risk BL, worst diagnosis HC vs. high-risk BL, worst diagnosis L + low-risk BL (group 3). Samples were" >> metilene_analysis.log
echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log

# ... high-risk BLH, worst diagnosis HC                                                                             vs. high-risk B, worst diagnosis L + low-risk BLH. also known as "group 2" in 4/7/24 email
#    [B16 B17 B26 B31 B32 L02 L05 L07 L08 L09 L10 L11 L12 L13 L17 L20 L21 L22 H01 H02 H04 H09 H10 H11 H13 H14 H15 H17]  [B02 B03 B04 B06 B08 B09 B11 B12 B13 B14 B15 B18 B20 B21 B22 B23 B24 B27 B28 B30 B31 B33 L03 L06 L14 L15 L18 L19 H03 H05 H06 H08 H12 H16]
head -1 all.beta.filt.samp_excl.tsv > header.tmp
for b in B16 B17 B26 B31 B32 L02 L05 L07 L08 L09 L10 L11 L12 L13 L17 L20 L21 L22 H01 H02 H04 H09 H10 H11 H13 H14 H15 H17; do
    c=`echo ${b} | sed 's/^/X/'`
    sed -i "s/${b}/${c}/" header.tmp
done

for d in B02 B03 B04 B06 B08 B09 B11 B12 B13 B14 B15 B18 B20 B21 B22 B23 B24 B27 B28 B30 B31 B33 L03 L06 L14 L15 L18 L19 H03 H05 H06 H08 H12 H16; do
    e=`echo ${d} | sed 's/^/Y/'`
    sed -i "s/${d}/${e}/" header.tmp
done

cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > higher_vs_lower_risk.tmp
~/tools/metilene/metilene_linux64 -a X -b Y -t 10 higher_vs_lower_risk.tmp > tmp
sort -k1,1 -k2,2n tmp | gzip > metilene_output/highBLHworstHC_vs_highBworstLlowBLH.group2.dmr.tsv.gz
rm -f tmp
echo "[$(date)] Completed high-risk BLH, worst diagnosis HC vs. high-risk B, worst diagnosis L + low-risk BLH (group 2). Samples were" >> metilene_analysis.log
echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log

# ... high-risk BLH, worst diagnosis C                                                      vs. high-risk BL, worst diagnosis LH + low-risk BLH. also known as "group 4" in 4/7/24 email
#    [B16 B17 B26 B31 L02 L05 L07 L08 L10 L11 L17 L21 H01 H02 H04 H09 H10 H11 H13 H14 H15 H17]  [B02 B03 B04 B06 B08 B09 B11 B12 B13 B14 B15 B18 B20 B21 B22 B23 B24 B25 B27 B28 B30 B32 B33 L03 L06 L09 L12 L13 L14 L15 L18 L19 L20 L22 H03 H05 H06 H08 H12 H16]
head -1 all.beta.filt.samp_excl.tsv > header.tmp
for b in B16 B17 B26 B31 L02 L05 L07 L08 L10 L11 L17 L21 H01 H02 H04 H09 H10 H11 H13 H14 H15 H17; do
    c=`echo ${b} | sed 's/^/X/'`
    sed -i "s/${b}/${c}/" header.tmp
done

for d in B02 B03 B04 B06 B08 B09 B11 B12 B13 B14 B15 B18 B20 B21 B22 B23 B24 B25 B27 B28 B30 B32 B33 L03 L06 L09 L12 L13 L14 L15 L18 L19 L20 L22 H03 H05 H06 H08 H12 H16; do
    e=`echo ${d} | sed 's/^/Y/'`
    sed -i "s/${d}/${e}/" header.tmp
done

cat header.tmp <(sed 1d all.beta.filt.samp_excl.tsv) > higher_vs_lower_risk.tmp
~/tools/metilene/metilene_linux64 -a X -b Y -t 10 higher_vs_lower_risk.tmp > tmp
sort -k1,1 -k2,2n tmp | gzip > metilene_output/highBLHworstC_vs_highBLworstLHlowBLH.group4.dmr.tsv.gz
rm -f tmp
echo "[$(date)] Completed high-risk BLH, worst diagnosis C vs. high-risk BL, worst diagnosis LH + low-risk BLH (group 4). Samples were" >> metilene_analysis.log
echo "[$(date)]" $(cat header.tmp | grep -Po "X.*?\t") vs $(cat header.tmp | grep -Po "Y.*?\t") >> metilene_analysis.log

# final cleanup
rm -f header.tmp
rm -f higher_vs_lower_risk.tmp
rm -f all.beta.filt.samp_excl.tsv
