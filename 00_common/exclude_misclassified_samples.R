#!/usr/bin/env Rscript

"> exclude_misclassified_samples.R <

After considering exploratory RNA-seq and EM-seq results, some samples look
obviously misclassified. Decision was made to exclude these samples. To reduce
duplication in code, this script aims to be 'source'-d into individual scripts,
so that every script would have the same samples excluded.

If the samples to be excluded gets changed, changing this single script is far
easier than modifying all individual scripts.

Currently, decision has been made to exclude:
  B01 B05 B07 B10
  L01 L04 L16
  H07 H18
  C01

B samples
  B01 - Contains 95% NSq and only 5% at maximum of IM. Both EM/RNA-seq look
        similar to NSq. Decision - EXCLUDE
  B05 - Hospital pathology report diagnosis was IM/BE. Pathologists' diagnoses
        IM or CM (cardiac mucosa/simple columnar mucosa). EMseq 'upper GI'.
        Decision - probably CM with pseudo-goblet cells. EXCLUDE
  B07 - Hospital pathology report diagnosis was EAC but cancer is not present in
        this tissue. Decision - EXCLUDE as tissue contains cardiac columnar
        mucosa (CM) only. EMseq consistent with CM.
  B10 - Patient previously had LGD. Uncertain if this sample is LGD or BE only. 

L samples
  L01 - Hospital pathology report diagnosis was LGD. Pathologists' diagnoses CM.
        EMseq 'upper GI'. Decision - EXCLUDE 
  L04 - Hospital pathology report diagnosis was LGD. Pathologists' diagnoses
        IM/indefinite for dysplasia or CM. EMseq 'upper GI'. Decision - EXCLUDE
        probably CM only.
  L16 - Hospital pathology report diagnosis was IM/BE. Pathologists' diagnoses
        LGD, IM and CM (cardiac mucosa/simple columnar mucosa). EMseq 'upper GI'.
        Decision - probably CM with pseudo-goblet cells. EXCLUDE

H samples
  H07 - Pathologist's assessment is suspicious of cancer. Confirmed diagnosis
        is HGD. Uncertain if HGD or EAC.
  H18 - Hospital pathology report diagnosis was HGD. Pathologists' diagnoses IM
        or IM+dysplasia except Duncan.
        Reassessment by Duncan - I think it is all columnar lined mucosa, no
        intestinal metaplasia, no dysplasia or malignancy. EXCLUDE

C samples
  C01 - Hospital pathology report diagnosis was EAC. Pathologists' diagnoses EAC,
        except Duncan.
        Reassessment by Duncan - Just oesophageal submucosa and muscularis
        mucosa. No mucosa (no epithelium), no IM, dysplasia or malignancy.
        EXCLUDE
" -> doc

excluded_samples <- c('B01', 'B05', 'B07', 'B10',
                      'L01', 'L04', 'L16',
                      'H07', 'H18',
                      'C01')

# df can be most objects, really. so long sample names are in rownames()
exclude_samples_by_row <- function(df) {
  excluded_rows <- rownames(df) %in% excluded_samples
  
  df[!excluded_rows, ]
}

# df can be most objects, really. so long sample names are in colnames()
exclude_samples_by_col <- function(df) {
  excluded_cols <- colnames(df) %in% excluded_samples
  
  if ('data.table' %in% class(df)) {
    df <- df[, !excluded_cols, with=FALSE]
  } else {
    df <- df[, !excluded_cols]
  }
  
  df
}
