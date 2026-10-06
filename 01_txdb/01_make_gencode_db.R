#!/usr/bin/env Rscript

"> 01_make_gencode_db.R <

When given a converted genome and GENCODE annotations, create an annotated
genome in the form sqlitedb and RData files.

Original repo (that processed hg19 stuff)
https://bitbucket.csiro.au/users/ros259/repos/txdb/browse
" -> doc

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicFeatures)
  library(rtracklayer)
  library(this.path)
  library(txdbmaker)
})


#####  Paths  #####
setwd(this.path::here())
path_base <- this.path::here()
path_data <- file.path(path_base, "../data")


#####  Functions  #####
extractAttributes <- function(df, threads=10) {
  
  attributesForRow <- function(x) {
    # Should always have a key-value pair split by space
    y <- strsplit(x, " ")
    z <- sapply(y, '[[', 2)
    names(z) <- sapply(y, '[[', 1)
    z
  }
  
  # Sometimes there is a space after the semi-colon
  x <- strsplit(df$attribute, "; ?")
  
  attributes <- parallel::mclapply(x, attributesForRow, mc.preschedule = TRUE,
                                   mc.cores = threads)
  attribute_set <- unique(unlist(lapply(attributes, names)))
  
  m <- matrix(NA, nrow=length(attributes), ncol=length(attribute_set),
              dimnames=list(1:length(attributes), attr=attribute_set))
  
  for(i in 1:length(attributes)) {
    a <- attributes[[i]]
    m[i, names(a)] <- as.vector(a)
  }
  data.frame(m, stringsAsFactors = FALSE)
}


rangesTxdb <- function(txdb=gencode, gtf, feature="gene") {
  
  stopifnot(feature %in% c("gene", "transcript"))
  
  if(feature == "gene") {
    gr <- genes(txdb)
  }
  if(feature == "transcript") {
    gr <- transcripts(txdb)
  }
  
  df <- gtf[gtf$feature == feature, ]
  df_attr <- extractAttributes(df)
  
  stopifnot(identical(nrow(df), length(gr)))
  
  merge_by <- NA
  last_mcol <- mcols(gr)[, ncol(mcols(gr))]
  for(i in 1:ncol(df_attr)) {
    is_same <- last_mcol %in% df_attr[, i]
    if(sum(is_same) == nrow(mcols(gr))) {
      merge_by <- i
      break
    }
  }
  
  attr_index <- match(last_mcol, df_attr[, i])
  mcols(gr) <- cbind(mcols(gr), df_attr[attr_index, ])
  gr
}


UCSCToGRanges <- function(UCSC_table, GRCh, seqnames_field="chrom",
                          start_field="chromStart", end_field="chromEnd") {
  
  gr <- makeGRangesFromDataFrame(UCSC_table, keep.extra.columns=TRUE,
                                 seqnames.field = seqnames_field,
                                 start.field=start_field,
                                 end.field=end_field,
                                 starts.in.df.are.0based=TRUE)
  # Prune so that only chr1... chrM is kept
  gr <- gr[seqnames(gr) %in% seqlevels(GRCh), ]
  seqlevels(gr) <- seqlevels(GRCh)
  seqinfo(gr) <-  GRCh
  gr
}


#####  Load files and make TxDb  #####

# Only use release 19 for hg19. http://www.gencodegenes.org/releases/19.html
#release <- 19
#release <- 36
release <- 43
genome_assembly <- "GRCh38.p13"

source_url <- paste("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                    release, "/gencode.v", release, ".annotation.gtf.gz", sep="")

gtf_suffix <- paste("gencode.v", release, ".annotation", sep="")
gtf_file <- file.path(path_data, paste(gtf_suffix, 'gtf.gz', sep="."))
rda_file <- file.path(path_data, paste(gtf_suffix, "_", Sys.Date(), '.RData', sep=""))
sqlite_file <- file.path(path_data, paste(gtf_suffix, 'sqlite', sep="."))

if(!file.exists(gtf_file)) {
  download.file(source_url, destfile = gtf_file, method="wget")
}

# Make the GENCODE GRanges object
GRCh = Seqinfo(genome=genome_assembly)
if(seqnames(GRCh)[1] == "1") {
  GRCh = GRCh[names(GRCh)[1:25]]
  seqnames(GRCh) = paste("chr", sub("MT", "M", seqnames(GRCh)), sep="")
}

gencode <- makeTxDbFromGFF(gtf_file,
                           format = "gtf", dataSource = source_url,
                           organism = "Homo sapiens",
                           chrominfo = GRCh)

saveDb(gencode, file=sqlite_file)


#####  Annotate  #####

gtf <- read.delim(gtf_file,
                  stringsAsFactors=FALSE, comment.char="#", header = FALSE,
                  colClasses = c("character", "character", "character",
                                 "integer", "integer", "character",
                                 "character", "character", "character"),
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attribute"))

attr(gtf, "source_url") <- source_url
attr(gtf, "gtf_file") <- gtf_file


gencode_gene <- rangesTxdb(txdb=gencode, gtf, feature="gene")
gencode_tx <- rangesTxdb(txdb=gencode, gtf, feature="transcript")

#gencode_cds <- cds(gencode)
#gencode_exon <- exons(gencode)
#gencode_intergenic <- gaps(gencode_gene)
#gencode_introns <- unlist(intronsByTranscript(gencode))


#####  CpGislands  #####

# Now create a UCSC session to hoover in CpGislands
session <- browserSession("UCSC")
if(release <= 19) {
  this_genome <- "hg19"  
} else {
  this_genome <- "hg38"
}
genome(session) <- this_genome

#ucscTables(this_genome, "CpG Islands")
#cpgIslands <- getTable(ucscTableQuery(session, track="CpG Islands", table="cpgIslandExt"))
cpgIslands <- getTable(ucscTableQuery(session, table="cpgIslandExt"))

## The start positions need to be converted into 1-based positions,
## to adhere to the convention used in Bioconductor:
cpgIslands <- UCSCToGRanges(cpgIslands, GRCh = GRCh)

#CpG island shores & 5kb
shore_size <- 2000

cpgShores <- setdiff(resize(cpgIslands, width(cpgIslands) + shore_size * 2,
                            fix="center"), cpgIslands)
cpgShores <- trim(cpgShores)

cpg5kb <- resize(cpgIslands, 5000)
cpg5kb <- trim(cpg5kb)


save(gencode, gencode_gene, gencode_tx, cpgIslands, cpgShores, cpg5kb,
     file=rda_file)  # gencode_cds, gencode_exon, gencode_intergenic, gencode_introns

#####  ENCODE cis-Regulatory Elements (cCREs)  #####

# Hoover in ENCODE Registry of candidate cis-Regulatory Elements (cCREs) 

#ucscTables(this_genome, "ENCODE cCREs")
#encodeCcre <- getTable(ucscTableQuery(session, track="ENCODE cCREs", table="encodeCcreCombined"))
encodeCcre <- getTable(ucscTableQuery(session, table="encodeCcreCombined"))

## The start positions need to be converted into 1-based positions,
## to adhere to the convention used in Bioconductor:
encodeCcre <- UCSCToGRanges(encodeCcre, GRCh = GRCh)


#####  Repeats  #####

# Hoover in Repeat sequences

#ucscTables(this_genome, "Simple Repeats")
#ucscTables(this_genome, "RepeatMasker")

#simple <- getTable(ucscTableQuery(session, track="Simple Repeats", table="simpleRepeat"))
#rmsk <- getTable(ucscTableQuery(session, track="RepeatMasker", table="rmsk"))
simple <- getTable(ucscTableQuery(session, table="simpleRepeat"))
rmsk <- getTable(ucscTableQuery(session, table="rmsk"))


## The start positions need to be converted into 1-based positions,
## to adhere to the convention used in Bioconductor:

simple <- UCSCToGRanges(simple, GRCh = GRCh)
rmsk <- UCSCToGRanges(rmsk, GRCh = GRCh, seqnames_field = "genoName",
                      start_field = "genoStart", end_field = "genoEnd")


#####  Save  #####

save(gencode, gencode_gene, gencode_tx, cpgIslands, cpgShores, cpg5kb, encodeCcre, simple, rmsk,
     file=rda_file)  # gencode_cds, gencode_exon, gencode_intergenic, gencode_introns


# for replicability purposes
sessionInfo()
