# `celfie_hg38_input/` folder #

This folder documents the steps needed to run a hg38-equivalent of the CelFiE-ISH pipeline. There are three key files used by the `deconvolution` script. I'll do my best to track and explain issues that I encountered (but did not affect the script running to completion). But let's start by looking at the deconvoluted form of the `deconvolution` command (I have preferentially used the expanded flag e.g., `--genomic_intervals` instead of `-i` to ease understanding/explanation).

`deconvolution --model celfie-ish --mixture <epiread_file> --cpg_coordinates <cpg_coordinates_file> --genomic_intervals <interval_file> -b --atlas_file <atlas_file> --epiformat old_epiread_A --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 --outfile <out_file>`

Note: no SNP-related files were used in generating input files.

## The `--mixture` file ##

`--mixture` is an epiread/pat file, converted from BAM, produced by mapping your reads against your genome of interest, that you want to deconvolute into the buckets outlined in the atlas file.

Catch is, the BAM has to be sorted and indexed. Assuming you have BAM files containing mapped reads ending with `*.deduplicated.bam` (in our case, produced by `deduplicate_bismark`):

```bash
# in these bash commands, ${a} comes from a loop that i ran over all the
# "*.deduplicated.bam" files in a folder. Feel free to set ${a} to a file of
# yours while testing these commands out, e.g.,
# $ set a=05.deduplicated.bam
$ samtools sort ${a} -o ${a/.deduplicated.bam/}.sorted.bam -@ 8
$ samtools index ${a/.deduplicated.bam/}.sorted.bam
```

> [!IMPORTANT]
> The epiread format demanded by CelFiE-ISH **IS NOT** identical to any standard epiread format produced by `biscuit epiread`.

Why this isn't in flashing bright neon lights on their page (and instructions provided to generate those files), I don't know.

Others have picked on it too https://github.com/methylgrammarlab/deconvolution_models/issues/12. Caused some grief on our end, had to stare at the files in the `demo/` folder in https://github.com/methylgrammarlab/deconvolution_models to figure things out.

The commands listed here is to generate something acceptable for `--epiformat old_epiread_A`; there are two other flags for alternate inputs, but I lack the will to figure out the proper commands for `--epiformat old_epiread` and `--epiformat pat` (and `pat` requires other files to be numbered differently too, can of worms that I won't open).

Assuming you have sorted BAM files ending with `*.sorted.bam`, and installed the `biscuit` tool somewhere locally:

```bash
# change directory to "epiread_files/", less messy that way
$ cd ../epiread_files

# the -O and -A flags ARE ABSOLUTELY NECESSARY. YES, "-O" IS OLD AND DEPRECATED
# BUT CELFIE-ISH DEMANDS THIS
$ ~/tools/biscuit epiread -@ 8 -O -A <your_fave_hg38_genome_fasta> ${a} > ${a/.sorted.bam/}.epiread

# the shell script (provided in the "epiread_files/" folder) adds two additional
# columns into the epiread file. absolutely necessary & absolutely UNDOCUMENTED
$ ./convert_biscuit_epiread_to_celfie.sh ${a/.sorted.bam/}.epiread > ${a}.tmp

# sort / bgzip / index / remove temp files
$ sort -k1,1 -k2,2n ${a}.tmp > ${a/.sorted.bam/}.sorted.epiread
$ bgzip -f ${a/.sorted.bam}.sorted.epiread
$ tabix -f -p bed ${a/.sorted.bam}.sorted.epiread
$ rm -f ${a}.tmp ${a/.sorted.bam/}.epiread; done
```

## The `--cpg_coordinates` file ##

This is the universe of CpG positions from the genome you're using. The Python3 script is provided in this folder.

```bash
# run the python script to produce plaintext BED
$ python3 create_cpg_coords_bed.py <your_fave_hg38_genome_fasta> > sorted_CpGs.hg38.bed
# bgzip & tabix
$ bgzip sorted_CpGs.hg38.bed
$ tabix -p bed sorted_CpGs.hg38.bed.gz
# this is what it should look like, roughly
$ ls -l sorted_CpGs.hg38.bed*
-rw-r--r-- 1 lie128 lie128 202M Feb  4 11:33 sorted_CpGs.hg38.bed.gz
-rw-r--r-- 1 lie128 lie128 1.7M Feb  4 14:36 sorted_CpGs.hg38.bed.gz.tbi
```

Do use the same hg38 genome that the epiread files were generated on. Consistency is king.

Please read the comments in `create_cpg_coords_bed.py`--it documents a quirk with chromosomal coordinates for CelFIE-ISH to work. In general, input files generated here should NOT be used in other tools/scripts as they're NOT standards compliant unfortunately.

`sorted_CpGs.hg38.bed.gz` IS NOT provided in the repo, it's too large for GitHub; ~200 MB gzipped. But here's a quick glance at that file, as a sanity check for yours:

```bash
$ zcat sorted_CpGs.hg38.bed.gz | head -5
chr1    10468   10469   CpG1
chr1    10470   10471   CpG2
chr1    10483   10484   CpG3
chr1    10488   10489   CpG4
chr1    10492   10493   CpG5
```

## The `--genomic_intervals` and `--atlas_file` files ##

Both files here are linked to each other, easier to lump them together and explain both at the same time.

The `--atlas_file` outlines which genomic loci is differentially methylated in a particular tissue / set of tissues. In other words, methylation information OUTSIDE of the atlas loci is pointless, as it doesn't help deconvolute the mix of epireads. As an analogy, let's say you bought a stack of magazines in a garage sale, and you wanted to deconvolute the pile of magazines by publisher and by year of publication. The key pages (in the atlas file) would most likely be "cover" and "first few pages", as the cover would likely have the issue date and the name of the publisher should be couple of pages past the cover.

In order to ignore all other pages but the first few... that's where `--genomic_intervals` comes in. Technically speaking, the loci in `--genomic_intervals` *should* be a **superset** of those in the `--atlas_file`. To convince yourself that that's the case, check out https://github.com/methylgrammarlab/deconvolution_models/tree/main/demo, and look for the "beta_atlas.txt" file and "U250.tsv".

Gah you know what I'll note it here anyways.

```bash
$ head beta_atlas.txt | cut -f 1-7  # note: wide file with per-CpG meth/cov
CHROM   START   END     Adipocytes_METH  Adipocytes_COV  Endothelium_METH  Endothelium_COV
chr1    1045636 1045637 39      42      272     281
chr1    1045689 1045690 48      56      340     383
chr1    1045697 1045698 56      59      351     389
chr1    1045700 1045701 55      56      363     386
chr1    1045787 1045788 60      72      458     510
chr1    1095821 1095822 14      16      107     113
chr1    1095931 1095932 47      55      279     313
chr1    1095940 1095941 60      62      313     322
chr1    1095947 1095948 49      63      290     336
$ sort -k1,1 -k2,2n U250.tsv | head  # note: file is not sorted initially
chr1    1045636 1045789 Colon-Ep
chr1    1095821 1096180 Colon-Ep
chr1    1197515 1197867 Eryth-prog
chr1    1336479 1336761 Liver-Hep
chr1    1435887 1436130 Lung-Ep-Alveo
chr1    1491549 1491951 Thyroid-Ep
chr1    1492054 1492375 Thyroid-Ep
chr1    1547175 1547251 Lung-Ep-Alveo
chr1    1560429 1560532 Heart-Cardio
chr1    1560959 1561055 Heart-Cardio
```

Notice that the first two rows of the U250 file contains all ten positions of the atlas. If you're inclined to dig deeper, check out those positions on a hg19 genome browser, and you might start noticing some... positional weirdness as well (which is why the hg38 versions inherit this weirdness).

The unmodified files hosted on https://github.com/nloyfer/UXM_deconv/tree/main/supplemental that will get molded into usable forms initially looks like:

```bash
$ head -5 Atlas.U250.l4.hg38.full.tsv | cut -f 1-12  # unfortunately the atlas file isn't per-CpG, hence the hoopla later
chr     start   end     startCpG        endCpG  target  name    direction       Adipocytes      Bladder-Ep      Blood-B Blood-Granul
chr1    185606  185761  2174    2180    Gallbladder     chr1:185606-185761      U       0.0     0.0     0.0     0.0
chr1    981709  981807  12952   12960   Skeletal-Musc   chr1:981709-981807      U       0.04    0.064   0.0     0.009
chr1    1003660 1003750 13886   13894   Epid-Kerat      chr1:1003660-1003750    U       0.0     0.0     0.0     0.0
chr1    1110257 1110409 18265   18270   Colon-Ep        chr1:1110257-1110409    U       0.0     0.102   0.0     0.011
$ sort -k1,1 -k2,2n Markers.U250.hg38.tsv | cut -f 1-7 | head -5  # this can be used, with minimal modification
#chr    start   end     startCpG        endCpG  target  region
chr1    185606  185761  2174    2180    Gallbladder     chr1:185606-185761
chr1    981709  981807  12952   12960   Skeletal-Musc   chr1:981709-981807
chr1    1003660 1003750 13886   13894   Epid-Kerat      chr1:1003660-1003750
chr1    1110257 1110409 18265   18270   Colon-Ep        chr1:1110257-1110409
```

### Downloading the raw Loyfer 2023 beta files ###

To Loyfer's credit, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458 mentioned in their paper really does have everything in hg19 and hg38.

The files I wanted were `*.hg38.beta` files in "GSE186458_RAW.tar". The beta files here aren't plaintext, it's a binary filetype that `nloyfer/wgbs_tools` is able to interpret (like `*.bam` with `samtools`). It's about 56 MB each, and about 10% of the entire mega-tar.

To get those download links in particular, I parsed it out of the "series matrix file" https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/matrix/ with

```bash
$ mkdir raw_beta_files
$ cd raw_beta_files

# replace tab with newlines (for grep later)    grep-with-regex                 remove dups
$ sed 's/\t/\n/g' GSE186458_series_matrix.txt | grep -Po "ftp.*?\.hg38\.beta" | sort | uniq | wc -l
253
# this matches the # samples in GSE186458

$ sed 's/\t/\n/g' GSE186458_series_matrix.txt | grep -Po "ftp.*?\.hg38\.beta" | sort | uniq | head -3
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652176/suppl/GSM5652176_Adipocytes-Z000000T7.hg38.beta
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652177/suppl/GSM5652177_Adipocytes-Z000000T9.hg38.beta
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652178/suppl/GSM5652178_Adipocytes-Z000000T5.hg38.beta

# the <() bash-ism creates a temp file with the sed command; feeds line-by-line into loop
$ while read a; do wget ${a}; done < <(sed 's/\t/\n/g' GSE186458_series_matrix.txt | grep -Po "ftp.*?\.hg38\.beta" | sort | uniq)
# LET'S GOOOOOOOO (and wait 2 hours)
```

### Filtering atlas-specific CpGs and compiling a hg38 atlas ###

Time to cut out the atlas-specific CpGs from these raw files!

Note: this section assumes that `@nloyfer/wgbs_tools` have been installed locally. My binary was at `~/tools/wgbs_tools/wgbstools`.

```bash
$ cd ..
$ mkdir filtered_beta_files
$ cd filtered_beta_files

# `wgbstools view` does not like header rows, need to remove it. sorting for
# the A E S T H E T H I C
$ sed 1d ../Markers.U250.hg38.tsv | sort -k1,1 -k2,2n > Markers.U250.hg38.bed

# ok so how many CpGs am i expecting
$ cut -f 8 Markers.U250.hg38.bed | sed 's/CpGs//' | awk '{s+=$1} END {print s}'
75158

# initial test on a single file to see whether output is expected
$ for a in ../raw_beta_files/GSM6810048_CNVS-NORM-WBC-110036152-WGBS-Rep1.hg38.beta; do b=`basename ${a}`; ~/tools/wgbs_tools/wgbstools view -L Markers.U250.hg38.bed ${a} > ${b/.hg38.beta/}.atlas_cpg.tsv; done
$ wc -l GSM6810048_CNVS-NORM-WBC-110036152-WGBS-Rep1.atlas_cpg.tsv
75158 GSM6810048_CNVS-NORM-WBC-110036152-WGBS-Rep1.atlas_cpg.tsv
# yeaaaaaaaaaaa

# repeat for all 253 files
$ for a in ../raw_beta_files/*.hg38.beta; do b=`basename ${a}`; ~/tools/wgbs_tools/wgbstools view -L Markers.U250.hg38.bed ${a} > ${b/.hg38.beta/}.atlas_cpg.tsv; done
# wait another 2 hours

# note: not all files have all 75,158 positions; but they aren't too far off
# thankfully. assume missingness as 0 meth 0 cov later
$ wc -l * | grep -v 75158
    74708 GSM5652190_Liver-Endothelium-Z000000RB.atlas_cpg.tsv
    74770 GSM5652226_Cortex-Neuron-Z0000042H.atlas_cpg.tsv
    74956 GSM5652344_Bladder-Epithelial-Z0000043F.atlas_cpg.tsv
    74816 GSM5652364_Gastric-body-Epithelial-Z000000ST.atlas_cpg.tsv
    74845 GSM5652367_Gastric-antrum-Epithelial-Z000000SR.atlas_cpg.tsv
    74855 GSM5652373_Colon-Left-Epithelial-Z000000VA.atlas_cpg.tsv
    75003 GSM5652381_Small-int-Epithelial-Z000000UY.atlas_cpg.tsv
    74722 GSM6810033_CNVS-NORM-WBC-110032789-WGBS-Rep1.atlas_cpg.tsv
    74802 GSM6810039_CNVS-NORM-WBC-110033745-WGBS-Rep1.atlas_cpg.tsv

# do note that while files have 75158 lines, there's only 75152 unique positions
# (likely from overlapping loci in the markers bed file)
$ cat GSM5652240_Pancreas-Duct-Z0000043U.atlas_cpg.tsv | cut -f 1-3 | sort -k1,1 -k2,2n | wc -l
75158
$ cat GSM5652240_Pancreas-Duct-Z0000043U.atlas_cpg.tsv | cut -f 1-3 | sort -k1,1 -k2,2n | uniq | wc -l
75152

# compress files for future upload of folder contents; but NOT "raw_beta_files/"
$ gzip *.atlas_cpg.tsv
```

Once that's done, run Alice's script `compile_atlas.R` (thanks!) on the compressed `*.atlas_cpg.tsv.gz` files, to produce the `Atlas.U250.per-CpG.hg38.tsv` file.

Due to space constraints, the ginormous `*.tsv.gz` files in `filtered_beta_files/` and `*.beta` files in `raw_beta_files/` are not uploaded here. (But none of these files are needed if you're not interested in replicating our pain.)

```bash
$ du -d 1
142M    ./filtered_beta_files
14G     ./raw_beta_files
```

### Accounting for CelFiE-ISH off-by-one weirdness ###

Now that we've generated correct BED files for the atlas, it's now time to slightly modify the atlas and marker files to feed CelFiE-ISH.

```bash
# marker file from nloyfer/UXM
$ sed 1d Markers.U250.hg38.tsv | sort -k1,1 -k2,2n | awk -F"\t" 'BEGIN {OFS=FS}{$2=$2-1; print}' > Markers.U250.hg38.celfie-ish.tsv
$ cut -f 1-8 Markers.U250.hg38.celfie-ish.tsv | head -5
chr1    185605  185761  2174    2180    Gallbladder     chr1:185606-185761      6CpGs
chr1    981708  981807  12952   12960   Skeletal-Musc   chr1:981709-981807      8CpGs
chr1    1003659 1003750 13886   13894   Epid-Kerat      chr1:1003660-1003750    8CpGs
chr1    1110256 1110409 18265   18270   Colon-Ep        chr1:1110257-1110409    5CpGs
chr1    1160441 1160800 20290   20304   Colon-Ep        chr1:1160442-1160800    14CpGs

# atlas file produced by Alice's R script
$ cat Atlas.U250.per-CpG.hg38.tsv | awk -F"\t" 'BEGIN {OFS=FS} FNR==1{print} FNR>1{$3=$3-1; print}' > Atlas.U250.per-CpG.hg38.celfie-ish.tsv
$ cut -f 1-10 Atlas.U250.per-CpG.hg38.celfie-ish.tsv | head -8
CHROM   START   END     Adipocytes_METH Adipocytes_COV  Bladder-Ep_METH Bladder-Ep_COV  Blood-B_METH    Blood-B_COV     Blood-Granul_METH
chr1    185605  185606  20      31      72      87      67      80      72
chr1    185609  185610  24      32      82      89      83      83      91
chr1    185691  185692  105     120     172     190     161     173     192
chr1    185727  185728  153     155     238     242     238     241     253
chr1    185730  185731  143     155     226     244     241     243     257
chr1    185759  185760  118     132     178     197     194     205     220
chr1    981708  981709  92      121     126     164     114     129     101
# note that the first loci has 6 CpGs, starting "185605" and ending "185761" in
# the marker file, but in the atlas file, the 6th CpG ends "185760". THIS IS
# **INTENTIONAL**, and follows the same pattern of the hg19 coords in CelFiE-ISH.
```

## Running CelFiE-ISH with hg38 input files ##

With all the prep work done, this command should finally work.

```bash
$ cd ..
$ deconvolution --model celfie-ish \
                --mixture epiread_files/05.sorted.epiread.gz \   # example input mixture file
                --cpg_coordinates celfie_hg38_input/sorted_CpGs.hg38.bed.gz \
                --genomic_intervals celfie_hg38_input/Markers.U250.hg38.celfie-ish.tsv -b \
                --atlas_file celfie_hg38_input/Atlas.U250.per-CpG.hg38.celfie-ish.tsv \
                --epiformat old_epiread_A --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 \
                --outfile celfie_hg38_output/05.deconvolute.txt  # example output mixture file
```

## Misc notes ##

Most of the large `*.tsv` files in this folder are gzipped-compressed, please decompress them whenever the commands ask for the plain files.
