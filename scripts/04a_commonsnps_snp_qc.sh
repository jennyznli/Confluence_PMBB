#!/bin/bash

## ALL ANCESTRIES: COMMON SNPS QC ##
module load plink/2.0

# Directory setup
HDIR="/project/knathans_tecac/jenny/breast"
WDIR="$HDIR/analysis/ALL/step1"
cd "$WDIR"
INPUT_FILE="common_snps_snp"

### SNP LEVEL FILTERS ###
plink2 --bfile "$INPUT_FILE" \
      --hwe 0.000001 \
      --maf 0.05 \
      --make-bed \
      --write-snplist --write-samples --no-id-header \
      --out "common_snps_qc"

### SAMPLE LEVEL FILTER ###
# plink2 --bfile "$INPUT_FILE" \
#      --mind 0.05 \
#      --make-bed \
#      --out "common_snps_mind"
# none removed
