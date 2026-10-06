# `12_deconv_celfieish/` folder #

This folder documents the ~~pain~~ effort involved in getting CelFiE-ISH to run on our data, which was mapped to a hg38 genome.

## Background ##

CelFiE-ISH (Unterman et al., Genome Biol, 2024; https://link.springer.com/article/10.1186/s13059-024-03275-x) works on hg19 inputs. I don't understand why the preference for a genome old enough to vote next year (2027), especially given its recent publication date.

The atlas powering CelFiE-ISH (Loyfer et al., Nature, 2023) thankfully has a hg38 version--it is a bit hidden though, and credit to Alice McAtamney for finding them. They can be found at https://github.com/nloyfer/UXM_deconv/tree/main/supplemental, specifically `Markers.U250.hg38.tsv`. In the same folder, `Atlas.U250.l4.hg38.full.tsv` looks tantalisingly useful--but upon inspection, the atlas stores the average beta in the marker loci. This is **NOT** the per-CpG beta/methylation/coverage info that CelFiE-ISH requires, and we were not able to find this anywhere (are they? please let me know if so).

The hacky method to generate the per-CpG atlas is to simply `liftOver` the hg19 version provided by CelFiE-ISH --> hg38 and call it a day, which was what I did before biting the perfection bullet. Weirdness included:

1. one CpG got lifted into an ALT chromosome, and crashes the `deconvolution` command when it's present; needed a `grep -v` to remove it manually.

2. some CpGs got lifted into the wrong area, sequence context was just wrong. I have receipts. One example was

   ```
   chr1:145688480 (hg19) --> chr1:145746604 (hg38)
        CCGTCTGAA                 CGGGGGGAG
   ```

   others have noted similar happenings, see https://groups.google.com/a/soe.ucsc.edu/g/genome/c/vQzWSWOQN50.

3. ~0.5% CpGs were unmapped (to be specific, hg19 atlas has 61,706 CpGs; liftOver + filtering for valid hg38 CpGs = 61,466 CpGs). I know, I know. This doesn't sound like it's a lot, but it's not **PERFECT** (and perfect ends up being 75,158 CpGs--I don't know why, but I guess perfection paid off?).

> [!IMPORTANT]
> Details on how the hg38-equivalent files were generated (with command line instructions & Python/R scripts) are available in the `celfie_hg38_input/` folder.


## Folder structure ##

Make sure you mimic the structure in this folder (including subfolders).

Quicker hack might be to `git clone` this entire repo (or this particular subfolder with `git sparse-checkout` but I've never tried that properly sorry), and get everything set up properly.


## Software required ##

### CelFiE-ISH ###

Repo is at https://github.com/methylgrammarlab/deconvolution_models. We used

```bash
$ pip install git+https://github.com/methylgrammarlab/deconvolution_models
```

### biscuit ###

Repo is at https://github.com/huishenlab/biscuit. We downloaded the precompiled binary for v1.7.1 (most recent version ~early 2026). Newer versions should work too.

```bash
$ wget "https://github.com/huishenlab/biscuit/releases/download/v1.7.1.20250908/biscuit_1_7_1_linux_amd64" -O "biscuit"
$ chmod 755 biscuit
```

> [!IMPORTANT]
> Note that the old epiread output from biscuit STILL REQUIRES MODIFICATION before it can be used as input. Shell script & steps needed to do so is in the `celfie_hg38_input/` folder.

### wgbs_tools ###

Repo is at github.com/nloyfer/wgbs_tools. Follow the instructions provided to install the tool.

```bash
$ git clone https://github.com/nloyfer/wgbs_tools.git
$ cd wgbs_tools
$ python setup.py
```


## Running CelFiE-ISH with hg38 input files ##

Assuming things went swimmingly, this command should finally work.

```bash
# run this command from this folder
$ deconvolution --model celfie-ish \
                --mixture epiread_files/05.sorted.epiread.gz \   # example input mixture file
                --cpg_coordinates celfie_hg38_input/sorted_CpGs.hg38.bed.gz \
                --genomic_intervals celfie_hg38_input/Markers.U250.hg38.celfie-ish.tsv -b \
                --atlas_file celfie_hg38_input/Atlas.U250.per-CpG.hg38.celfie-ish.tsv \
                --epiformat old_epiread_A --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 \
                --outfile celfie_hg38_output/05.deconvolute.txt  # example output mixture file
```


## Merging output files ##

As CelFiE-ISH runs on individual files, merge them on the command line so that I can start eyeballing patterns in the data with the best data exploration tool: Microsoft Excel (really!). Plus, it's also much easier to ask R to read a single tsv than run a loop and merge 86 files together in R.

```bash
#      [file containing a single line for the header. add "emseq_sample_no <TAB>" to the start of the line                                                                               ]  [file containing content rows: add sample number via the FILENAME variable, then sed away unwanted bits                         ]
$ cat <(head -1 celfie_hg38_input/Atlas.U250.per-CpG.hg38.celfie-ish.tsv | sed 's/\t/\n/g' | grep '_METH' | sed 's/_METH//' | tr "\n" "\t" | sed 's/$/\n/' | sed 's/^/emseq_sample_no\t/') <(awk '{print FILENAME"\t"$0}' celfie_hg38_output/*.deconvolute.txt | sed 's|^celfie_hg38_output/||' | sed 's/.deconvolute.txt//') > all_86_samples.celfie-ish_deconv.tsv
```
