# Methylation changes in Barrett's metaplasia, dysplasia, and cancer #

## Background ##

The incidence rates of esophageal adenocarcinoma (EAC), one of cancers with the highest case-fatality ratios, has increased by sevenfold over the last half century. Clinical management of EAC, dependent on the accurate diagnosis of metaplasia or dysplasia in Barrett’s esophagus (BE), is hampered by the unreliable histopathological diagnosis of dysplasia and potential sampling errors.

This repository stores the code written to obtain whole-genome methylation profiles from samples across the disease progression spectrum, notably supplementing the scant existing methylation evidence for low-grade (LGD) and high-grade dysplasia (HGD).

## Technical description ##

This repo stores the Python/R scripts written for this project; but is not complete, as the larger data files cannot be uploaded to GitHub.

Combine this repo with data hosted on CSIRO DAP (https://data.csiro.au/collection/78585), which stores the larger data files needed to run some of the scripts in this repo.

For completeness, the raw BAM files have been uploaded to SRA, TODO INSERT SRA LINK HERE. Note that none of the scripts here depend on the BAM files--the BAM files were used to generate the per-position cov files with `bismark_methylation_extractor`, and everything downstream uses those cov files.

American English is (was meant to be?) used throughout--not exactly my preference--as the project team agreed to produce figures with E for "esophagus" (vs. O-equivalents for "oesophagus") for consistency across our manuscripts. Apologies if you notice inconsistency in spelling, it's hard to switch my ~~favored~~ *favoured* spelling preferences overnight.

## Publication ##

Currently on bioRxiv, TODO INSERT BIORXIV LINK.

## Folder structure ##

Combining contents across CSIRO DAP + GitHub, the repo should look like this. Folders sans comments are present in this repo.

Rule of thumb is that scripts in later folders may require those in previous folders to be run to generate intermediate files; but never the reverse.

```bash
$ tree -L 1
.
├── 00_common
├── 01_txdb
├── 02_bismark_per-pos_cov     ## on CSIRO DAP
├── 03_bismark_per-cpg_cov     ## on CSIRO DAP
├── 04_filter_cpgs             ## on CSIRO DAP
├── 05_qc_post-CpG_filter
├── 06_deconvolute_meth_atlas  ## "*.csv.gz" on CSIRO DAP
├── 07_filter_samples          ## "*.tsv.gz" on CSIRO DAP
├── 08_qc_post-sample_filter   ## "*.rds" on CSIRO DAP
├── 10_call_dmr_metilene       ## metilene_output/ on CSIRO DAP
├── 11_analyse_dmr_metilene    ## data/ and filtered_dmrs/ on CSIRO DAP
├── 12_get_cnv                 ## 01_run_cnvkit/ on CSIRO DAP
├── 13_meth_vs_expr
├── 14_deconv_celfieish
├── 15_dmr_by_cell_types       ## filtered_dmrs/ and metilene_output/ on CSIRO DAP
├── 16_expr_by_cell_types
├── data                       ## on CSIRO DAP
└── README.md
```

## Short descriptions of folders ##

### Data used in multiple folders ###

`data`: stores publicly available resources (genomes, SNP tables, ...).

`00_common`: stores scripts/tables used in other subfolders.

`01_txdb`: family of scripts that provides genomic annotations for genomic loci, courtesy Jason Ross.

### Results from the `bismark` pipeline + QC ###

`02_bismark_per-pos_cov`: `bismark_methylation_extractor --bedGraph --gzip` on the BAM files uploaded to SRA, keeping the *.cov.gz files and discarding everything else.

`03_bismark_per-cpg_cov`: `coverage2cytosine --merge_CpG --gzip --genome_folder ../data/grch38p13_lambda_puc/ --dir per-cpg_cov_files/` on the cov files in `02`. In plain English, this script (part of `bismark`) checks the genome and merges meth/unmeth reads from both cytosines in the same CpG together. CpG methylation in humans is mostly symmetrical, and helps reduce the dataset by half (to a more manageable ~29M CpGs).

`04_filter_cpgs`: the shell script `filter_cpgs.sh` contains the commands used & relevant outputs used to pick out CpGs that were decently covered in a majority of samples. Drops the valid CpGs to 25M from 29M.

`05_qc_post-CpG_filter`: check whether CpG filtering worked well. Bunch of sequential R scripts.

`06_deconvolute_meth_atlas`: check whether samples were correctly labelled. See `deconvolute_samples.sh` for the commands used.

`07_filter_samples`: remove samples that weren't. See `filter_samples.sh` for the commands used.

`08_qc_post-sample-filter`: basically a repeat of `05`, with post-filtered samples. Bunch of sequential R scripts as well.

### FUN STARTS HERE ###

`09` being excluded is INTENTIONAL (really), for me to mentally bin `0*/` folders as the preprocessing stuff, and `1*/` folders as the proper analysis.

`10_call_dmr_metilene`: uses `metilene` to call DMRs. 

`11_analyse_dmr_metilene`: bunch of R scripts to process the `metilene` calls and separate them into hypermethylated and hypomethylated markers. For the `0` series scripts, two comparisons have separate scripts (HC vs. N and HC vs. B) as there were existing literature that dealt with these two comparisons, and I wanted to check whether our comparisons overlapped others. For the `1` series of scripts, I tried different ways to get reasonable functional enrichments linked to the DMR. What tends to work well for methylation arrays (top ~1% most variable --> functional annotation) **DOES NOT** work well for sequencing-based data.

`12_get_cnv`: deals with anything copy number related in the manuscript. has a separate `README.md` to describe the running of `cnvkit` on our samples. Note that the tool runs on coord-sorted BAM files, which is NOT uploaded to SRA (it's ~3 TB). The sorted files can be generated with `samtools sort && samtools index`, but unsorting sorted files are... a bit more convoluted.

`13_meth_vs_expr`: courtesy RNAseq data being made available by Udumanne et al. (2026), first published in bioRxiv (https://www.biorxiv.org/content/10.64898/2026.08.25.746910v1). Correlated methylation levels with expression levels. Far less complex, no README, just a bunch of scripts that can be run sequentially.

`14_deconv_celfieish`: deconvoluting bulk methylation patterns with CelFiE-ISH. Running this tool was far from straightforward (don't get me wrong, I really appreciate the authors' efforts in creating this tool, I'm nowhere as smart), please refer to the `README.md` in this folder to find out how I crowbarred things to work on hg38. **There are hg38-equivalent atlas files you can download from this folder, to shortcut this process!!**

`15_dmr_by_cell_type`: do existing markers make sense on a cell-type-based view of EAC methylation patterns? Run the scripts sequentially, no README needed.

`16_expr_by_cell_type`: similar to `13` really, but with a cell-type lens. Just one script in this folder.

The figure/table numbering in each folder, produced by the respective scripts, should correspond correctly to the figures/tables in the manuscript. Hopefully the changes aren't too drastic during peer review, would save me heaps of time/effort renumbering stuff.
