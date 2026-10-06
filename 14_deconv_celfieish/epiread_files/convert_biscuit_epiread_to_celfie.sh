#!/bin/bash

# > convert_biscuit_epiread_to_celfie.sh <
#
# CelFiE-ISH has an input format called "old_epiread_A", and that is what's
# used in their vignettes and in their demo/mixture.epiread.gz file--but
# what is NOT obvious is that it's a modified form of biscuit's output
# (and to be precise, it's the older biscuit output of `biscuit epiread -O -A`).
# This difference is not documented, and hugely annoying.
#
# Example of their input file:
#
# chr1  1045481  1045700  A00160:227:HL2H3DMXX:2:1187:27172:34585  2  -  1045481,1045505,1045508,1045519,1045537,1045636,1045689,1045697,1045700  C-CCT-CTT  .  .  Pancreas-Beta
# chr1  1045605  1045700  A00172:278:HHVCKDMXX:1:1156:22471:23625  1  +  1045636,1045689,1045697,1045700 CCTC  .  .  Adipocytes
# chr1  1045636  1045636  A00167:234:HL3C3DMXX:2:2437:3821:10175   1  -  1045636  C  .  .  Neuron
#
# Example of unmodified `biscuit epiread -O -A` output (with no SNP-containing
# bed file provided):
#
# chr1  A00121:904:HHN5KDSX7:4:2565:27697:27696_1:N:0:TCCTCATG+AGGTGTAC  1  +  54353   T
# chr1  A00121:904:HHN5KDSX7:3:1569:26241:1814_1:N:0:TCCTCATG+AGGTGTAC   2  +  99165,99217  C-
# chr1  A00121:904:HHN5KDSX7:1:2211:20907:24925_1:N:0:TCCTCATG+AGGTGTAC  2  +  108774  C
#
# This shell script adds in columns 2/3 (start/end), to facilitate sorting
# and tabix-ing downstream (which is again not mentioned in their docs...).
# Underlying logic is: check column 5. Leftmost value is start, rightmost end.
#
# Script expects arguments as follows: 
#   convert_biscuit_epiread_to_celfie.sh <non-compressed epiread file>
# (only the first argument is recognised)
#
# Prints to standard output.
#
# Adapted code courtesy Alice McAtamney.

awk 'BEGIN {{ FS=OFS="\t" }}
     {{
       n=split($5,a,","); if(n<1) next;
       start=a[1]; end=a[n];
       print $1,start,end,$2,$3,$4,$5,$6
     }}' ${1}
