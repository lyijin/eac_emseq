# `13_get_cnv/` folder #

This folder documents the steps to replicate the running of CNVkit on our data (surprisingly painless!). GATK was far more painful--and the biggest con that I eventually realised, was that there wasn't an in-built way to visualise copy number changes across multiple samples; whereas CNVkit has a full page devoted to explaining what it could do https://cnvkit.readthedocs.io/en/stable/plots.html.

Was very appealing to the lazy coder in me; CNVkit it is.


## Background ##

Cite & read this: CNVkit (Talevich et al., PLOS Comp Biol, 2016; https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873).

No other gotchas!


## Folder structure ##

Make sure you mimic the structure in this folder (including subfolders).

Quicker hack might be to `git clone` this entire repo (or this particular subfolder with `git sparse-checkout` but I've never tried that properly sorry), and get everything set up properly.


## Software required ##

### CNVkit ###

Repo is at https://github.com/etal/cnvkit/. Opted to use the `miniconda` install method.

```bash
$ conda create -n cnvkit -c bioconda
$ conda activate cnvkit
$ conda install cnvkit
```

> [!IMPORTANT]
> The `cnvkit` conda package contains its own preferred version of Python, R and associated libraries. `import cnvkit` (Python) would work in this environment, but this environment does not have e.g., R's tidyverse. Some scripts in this folder needs to be run with the `cnvkit` environment; some scripts have to be run OUTSIDE the environment. Ye be warned.

### Additional flat files ###

Well technically not software, but download these files somewhere. Only relevant if you're working on human genomes. It's such a privilege to work with model organisms backed with massive genomic resources \o/

```bash
$ wget https://github.com/etal/cnvkit/raw/refs/heads/master/data/refFlat_hg38.txt
$ wget https://github.com/etal/cnvkit/raw/refs/heads/master/data/access-10kb.hg38.bed
```

## Prepping BAM files for `cnvkit` ##

CNVkit runs on coord-sorted BAM files. It is HIGHLY HIGHLY recommended to rename your files at this stage to what you'd like to see on the y-axis of the resulting plots. I lost hair trying to fight `cnvkit` in wanting the y-axis to only have the abbreviated sample names, the library's... not particularly customisable in its visualisations (I do appreciate them making the viz easy to generate though; it's just hard to tweak the defaults).

To do so, I used symbolic links that looked like XXX.bam to point to the actual, expanded, full-fat filename.

As an example:

```bash
$ ls -l 00_sorted_bam/ | head -5
lrwxrwxrwx 1 lie128 lie128 110 Apr  9 12:19 B02.bam -> /scratch1/lie128/03b_sorted_dedup_grch38p13_lambda_puc/11_GTACACCT-CATGAGGA_R1_val_1_bismark_bt2_pe.sorted.bam
lrwxrwxrwx 1 lie128 lie128 114 Apr  9 12:19 B02.bam.bai -> /scratch1/lie128/03b_sorted_dedup_grch38p13_lambda_puc/11_GTACACCT-CATGAGGA_R1_val_1_bismark_bt2_pe.sorted.bam.bai
lrwxrwxrwx 1 lie128 lie128 110 Apr  9 12:19 B03.bam -> /scratch1/lie128/03b_sorted_dedup_grch38p13_lambda_puc/12_CGGCATTA-TGACTGAC_R1_val_1_bismark_bt2_pe.sorted.bam
lrwxrwxrwx 1 lie128 lie128 114 Apr  9 12:19 B03.bam.bai -> /scratch1/lie128/03b_sorted_dedup_grch38p13_lambda_puc/12_CGGCATTA-TGACTGAC_R1_val_1_bismark_bt2_pe.sorted.bam.bai
```


## Running `cnvkit` ##

RTFM: https://cnvkit.readthedocs.io/en/stable/quickstart.html.

I preferred having the output files contained in a folder.

```bash
$ mkdir 01_run_cnvkit
$ cd 01_run_cnvkit
$ cnvkit.py batch ../00_sorted_bam/B*.bam ../00_sorted_bam/L*.bam ../00_sorted_bam/H*.bam ../00_sorted_bam/C*.bam --normal ../00_sorted_bam/N*.bam --method wgs --drop-low-coverage --processes 12 --fasta ../../data/grch38p13_lambda_puc/grch38p13_lambda_puc.fa --annotate ../../data/refFlat_hg38.txt --access ../../data/access-10kb.hg38.bed --short-names --target-avg-size 20000
```

Additional notes:

My normal samples are all prefixed by N__. `cnvkit` generates a method-specific coverage file per normal sample--but taking a step back, why do we need to do this at all? It's because any blah-seq method would have biases in coverage tied to genome context. The normal samples thus becomes the baseline "expected" coverage, but it also means that normals, by definition, **contain NO CNVs**.

`--method wgs`, because, uh, no baits were used?

`--drop-low-coverage` was recommended to reduce the noisy fluctuations in copy number in low-coverage bins, which plague samples with low coverage (I do have some unfortunately).

`--processes 12` choo choo multiprocessing.

`--fasta [...]/grch38p13_lambda_puc.fa` self explanatory. But why, I hear a plaintive cry, does the genome have lambda and pUC? It's EM-seq-specific, as they are controls to check for conversion rates. The non-human chromosomes disappear anyway, because of the `--access` thing (read on).

`--annotate [...]/refFlat_hg38.txt` annotation file provided by cnvkit. Hope you've downloaded this file somewhere, I've provided the file location earlier in this README.

`--access [...]/access-10kb.hg38.bed` not every bit of the human genome is accessible. This file tells cnvkit where (chrom:start-end) it should be computing coverages for. And this is where non-human weird sequences disappear, as it's... not defined as being accessible.

`--short-names` summarises gene names. Brevity is the soul of wit. Unfortunately not something I adhere to when writing READMEs.

`--target-avg-sze 20000` this controls how large the genomic buckets are to calculate coverages for. This is something that needs a bit of trial-and-error--docs recommend to start at 5000 for WGS. Number has to go higher if samples have poorer coverage. Google and decide for yourself yeah. I tried 10 kb but some of the samples looked very stripey (very noisy) in the visualisations, much better at 20 kb. Not much difference with 40 kb.

Wait overnight for code to run.

This command produces all the `*.cnn *.cnr *.cns` files in this folder. I've gzipped them for upload reasons (please decompress if you want to replicate what I did, and for the downstream commands to work).


## Creating the viz used in the manuscript ##

There's more viz that can be generated with standard tools, but I cared more about the heatmap as I wanted to compare/contrast across samples (and always arranging them in the disease progression continuum).

I very quickly realised that mean coverage plays a HUGE role in producing meaningful results for any CNV analysis. For our samples, I decided to draw the line at coverage > 12. I did store the mean coverages in a flat file (`heatmap_inclusion.tsv`), used bash to read the file and plot a heatmap on samples with coverage above that threshold.

```bash
$ cd ..
$ mkdir 02_cnvkit_viz
$ cd 02_cnvkit_viz
$ b=`for a in $(sed 1d ../heatmap_inclusion.tsv | awk '{if ($4 > 12) print}' | grep -v '^N' | cut -f1); do printf -- "../01_run_cnvkit/${a}.cns "; done` && cnvkit.py heatmap ${b} --no-shift-xy -o heatmap.BLHC.pdf
```

This command produces the sole pdf in the folder (Supplementary Fig. 2 of the manuscript). Remember, N samples are all defined as not having any CNVs (so no point plotting them out).

I used `--no-shift-xy` because CNVkit was a bit atrocious in auto-detecting the sexes of our samples (the shading of X/Y depends on sex). Not sure why it kept detecting some of the Ms as Fs, and there's no way to override the auto-detection (you can, but you can only declare ALL samples as being ALL of a certain sex. I think).


## Producing the other figures in this manuscript ##

I used the `cnvkit.py call` command to call predicted copy numbers of the amplified/deleted segments.

```bash
$ cd ..
$ mkdir 03_cnvkit_call_cnv
$ cd 03_cnvkit_call_cnv
$ for a in $(sed 1d ../heatmap_inclusion.tsv | awk '{if ($4 > 12) print}' | grep -v '^N' | cut -f1); do cnvkit.py call ../01_run_cnvkit/${a}.cns -o ${a}.call.cns; done
```

Once this is prepped, the scripts in the root folder can be run sequentially.

```bash
$ cd ..
# i hate writing Rmd files, my hack is to render R --> HTML (and why i comment the heck out of my R scripts)
## RUN THESE R SCRIPTS OUTSIDE THE CNVKIT CONDA ENVT ##
$ for a in 01_*.R 03_*.R; do Rscript -e "rmarkdown::render('"${a}"', output_format=rmarkdown::html_document(df_print='paged', highlight='zenburn'))"; done

## RUN THE FIRST PYTHON SCRIPT OUTSIDE THE CNVKIT CONDA ENVT. GONNA ASSUME YOU HAVE NUMPY AND PANDAS ##
$ python3 02a_calc_goi_median_beta.py

## RUN THE SECOND PYTHON SCRIPT INSIDE CNVKIT CONDA ENVT. INSIDE!!! ##
$ python3 -i 02b_check_cnv_vs_goi.py
```

These commands generate the other figs in the manuscript (Fig. 3A, 3B, Supplementary Fig. 3).

`matplotlib` is powerful but shows its age, it's clunky and unintuitive sometimes. I used Affinity Design to get some REALLY hard-to-code legends & lines demarcating disease stages into the figs. Everything is visual sugar, no changes were done to the actual meat of the figure.
