# `01_txdb`: convert gencode GTF into R-friendly sqlitedb/RData format #

Git cloned from https://bitbucket.csiro.au/scm/~ros259/txdb.git, and modified to work on GENCODE v43.

## GENCODE details ##

GENCODE release 43 link\
https://www.gencodegenes.org/human/release_43.html

File downloaded from GENCODE\
http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

## What do ##

If you have the required libraries (check out `01_make_gencode_db.R`), run the script in RStudio or on the commandline ("Rscript \<script_name\>").

Run `02_test_annotate_ranges.R` to confirm things are working.

## Output ##

Produces the `.sqlite` and `.RData` files in the `../data/` folder.

HTML outputs show library versions used in the scripts.
