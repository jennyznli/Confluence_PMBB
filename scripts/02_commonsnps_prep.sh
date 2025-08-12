#!/bin/bash

### ASSIGN SEX TO COMMON SNPS ### 

module load plink/1.9-20210416

# directory setup
HDIR="/project/knathans_tecac/jenny/breast/analysis"
WDIR="$HDIR/ALL/step1"
ID_DIR="$HDIR/ids"
mkdir -p "$WDIR"
cd "$WDIR"

PHE_DIR="/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0"
COM_DIR="/static/PMBB/PMBB-Release-2024-3.0/Imputed/common_snps_LD_pruned"

# input files
COV="$PHE_DIR/PMBB-Release-2024-3.0_covariates.txt"
SNP_FILE="$COM_DIR/PMBB-Release-2024-3.0_genetic_imputed.commonsnps"

cp "$SNP_FILE.bim" common_snps.bim
cp "$SNP_FILE.bed" common_snps.bed

awk 'NR > 1 {print $1, $6}' "$COV" > "$ID_DIR/all_ids_sex.txt"

awk '
BEGIN { OFS="\t" }
NR==FNR {
    sex_code = ($2 == "Male" ? 1 : ($2 == "Female" ? 2 : 0));
    sex[$1] = sex_code;
    next
}
{
    if ($2 in sex) {
        $5 = sex[$2]
    } else {
        $5 = 0  # Unknown sex
    }
    print $1, $2, $3, $4, $5, $6
}' "$ID_DIR/all_ids_sex.txt" "$SNP_FILE.fam" > common_snps.fam

echo "Sex distribution in updated file:"
awk '{print $5}' common_snps.fam | sort | uniq -c

