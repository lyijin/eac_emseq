#!/usr/bin/env python3

docstring = """
> exclude_misclassified_samples.py <

After considering exploratory RNA-seq and EM-seq results, some samples look
obviously misclassified. Decision was made to exclude these samples. To
facilitate removal of misclassified samples from tsv files used downstream
in a uniform and documented manner, this script takes in plaintext tabular
data and returns the same data with the unwanted columns removed.

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

Autodetects gzip compression.
""".strip()

MISCLASSIFIED_SAMPLES = ['B01', 'B05', 'B07', 'B10', 
                         'L01', 'L04', 'L16',
                         'H07', 'H18',
                         'C01']

import argparse
import csv
from pathlib import Path
import glob
import gzip
import sys
import time

def diag_print(*args):
    """
    Helper function to print diagnostic-level verbose output to stderr.
    """
    print (f'[{time.asctime()}]', *args, file=sys.stderr)


parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('beta_tsv', metavar='tsv_files', type=Path,
                    help='Whole genome methylation data, in tsv format.')
parser.add_argument('--verbose', '-v', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# gzip autodetection
if args.beta_tsv.suffix == '.gz':
    # assume gzip-compressed file
    tsv_reader = csv.reader(gzip.open(args.beta_tsv, 'rt'), delimiter='\t')
else:
    # assume it's plaintext
    tsv_reader = csv.reader(open(args.beta_tsv, 'r'), delimiter='\t')

if args.verbose:
    diag_print(f'Reading {args.beta_tsv.name}...')

# check header for misclassified samples, and which column they're at. assumes
# columns have unique names (fairly safe assumption...)
header = next(tsv_reader)
unwanted_cols = {}
for ms in MISCLASSIFIED_SAMPLES:
    if ms in header:
        unwanted_cols[ms] = header.index(ms)

assert unwanted_cols, \
    'Script aborted, no misclassified samples detected in tabular file.'

if args.verbose:
    diag_print(f'Removing columns {unwanted_cols}...')

# remove unwanted cols. pop() works best from largest integer to smallest
unwanted_col_vals = sorted(unwanted_cols.values(), reverse=True)
for ucv in unwanted_col_vals:
    header.pop(ucv)

print (*header, sep='\t')

for row in tsv_reader:
    if len(row) < unwanted_col_vals[0]: continue
    
    for ucv in unwanted_col_vals:
        row.pop(ucv)
    
    print (*row, sep='\t')
