#!/usr/bin/env python3

docstring = """
> create_cpg_coords_bed.py <

CelFiE-ISH requires a input bed file for the --cpg_coordinates flag, and that
bed file is generated from grepping the input FASTA genome (plain-text) for
CG dinucleotides, and sequentially numbering them with CpG1, CpG2, ...

RANT: the coords that CelFiE-ISH desires isn't very logical unfortunately.
BED files are typically 0-based, i.e., the first CpG in hg38 should be
chr1 \t 10468 \t 10470 (and end minus start should be 2). If this was meant
to be 1-based, then the first CpG should be chr1 \t 10469 \t 10470... but what
CeLFiE-ISH wants is chr1 \t 10468 \t 10469.

MAKE IT MAKE SENSE RAAAAAAHHHHHHHH.

I mean, I think get it. CelFiE-ISH is a Python script so arrays are 0-based,
and when the genome is saved into an array, the first CpG can be accessed with
array['chr1'][10468] and array['chr1'][10469] respectively. But ugh I wish
people cared more about logical consistencies in input files.

ANYWAY TL;DR DON'T USE THIS SCRIPT AS A GENERIC WAY TO PARSE CpGs FROM A GENOME
(OR MODIFY IT TO ACCOUNT FOR THE INTENTIONAL OFF-BY-ONE ERROR). YE BE WARNED.
"""
import argparse
import re

import parse_fasta

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('fasta_file', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='plain-text genome FASTA.')

args = parser.parse_args()

# read all sequences into memory (yep, memory-inefficient), but this removes
# all newlines (if they're there) separating the genome sequences
sequences = parse_fasta.get_all_sequences(args.fasta_file, 'fasta')

# need a CpG counter for fourth column
cpg_counter = 0

# iterate through sequences, find CGs, and print em out (chrom-by-chrom)
for chrom in sequences:
    cgs = re.finditer(r'CG', sequences[chrom])
    start_positions = [x.start() for x in cgs]
    
    for sp in start_positions:
        cpg_counter += 1
        print (chrom, sp, sp+1, f'CpG{cpg_counter}', sep='\t')
        # to produce a correct 0-based BED referring to a CpG dinucleotide,
        # replace "sp+1" with "sp+2"
