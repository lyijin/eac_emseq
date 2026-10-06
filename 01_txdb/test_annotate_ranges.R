library(testthat)
library(GenomicRanges)

path_base = file.path("", "datasets", "work", "hb-multicancer-dpcr", "work",
                      "TxDB")
path_data <- file.path(path_base, "data")
path_code <- file.path(path_base, "scripts")

source(file.path(path_code, "annotate_ranges_functions.R"))

query <- GRanges(seqnames = c("chr2", rep("chr1", 9)),
                 ranges = IRanges(start=c(1, 1, 1, 1, 81, 1, 201, 50, 42, 50),
                                  width=c(100, 100, 100, 100, 20, 20, 100, 100, 1, 20)),
                 strand = c("+", "+", "-", "*", "+", "+", "+", "-", "+", "-"))

subject <- GRanges(seqnames = c(rep("chr1", 4)),
                 ranges = IRanges(start=c(60, 60, 50, 40),
                                  width=c(100, 50, 50, 30)),
                 strand = c("+", "+", "+", "-"))


test_that("coverage", {
    
    # Case for non-matching chromosome, should have 0 overlap
    expect_identical(coverageRatio(query[1], subject, ratio = FALSE), 0)
    
    # Cases for same region but different strand information in the subject
    expect_identical(coverageRatio(query[2:4], subject, ratio = FALSE), rep(61, 3))
    expect_identical(coverageRatio(query[2:4], subject[c(1, 3, 4)], ratio = FALSE), rep(61, 3))
    
    # Test for correct calcution of ratio versus width output
    expect_identical(coverageRatio(query[5:10], subject, ratio = FALSE), c(20, 0, 0, 100, 1, 20))
    expect_identical(coverageRatio(query[5:10], subject, ratio = TRUE), c(1, 0, 0, 1, 1, 1))

})