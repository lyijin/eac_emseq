#!/bin/bash

# > rename_sorted_bam.sh <
#
# CNVkit is much easier run when filenames inherit the short IDs (Nxx, Bxx, etc)
# hence this shell script aims to carry out the renaming with symbolic links
# pointing to the sorted bam files with the old numeric IDs.

for a in /scratch1/lie128/03b_sorted_dedup_grch38p13_lambda_puc/*.sorted.bam
#for a in ../../02_bismark_per-pos_cov/03b_sorted_dedup_grch38p13_lambda_puc/*.sorted.bam
do
  num_id=`basename ${a} | cut -c1-2`
  short_id=`cut -f 5,8 ../../00_common/emseq-rnaseq_clin_details.240716.tsv | grep ^${num_id} | cut -f 2`
  
  ln -s ${a} ${short_id}.bam
  ln -s ${a}.bai ${short_id}.bam.bai
done
