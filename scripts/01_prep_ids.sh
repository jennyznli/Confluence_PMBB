#!/bin/bash

### PREPARE ID FILES ###

HDIR="/project/knathans_tecac/jenny"
WDIR="$HDIR/breast/analysis/grafanc"
INPUT_FILE="$HDIR/grafanc/results/ancestry_final.txt"
ID_DIR="$HDIR/breast/analysis/ids"

CASES_FILE="$HDIR/phenotype/breast_cancer_filtered_patients_ids.txt"
CONTROLS_FILE="$HDIR/phenotype/non_cancer_female_nocrep_ids.txt"
UNREL_FILE="/static/PMBB/PMBB-Release-2024-3.0/Exome/IBD/PMBB-Release-2024-3.0_genetic_exome.3rd_degree_unrelated.txt"

mkdir -p "$ID_DIR"
cd "$WDIR"

### INTERSECTIONS ###
tail -n +2 "$INPUT_FILE" | cut -f1 | sort -u > "grafanc_ids.txt"
comm -12 "grafanc_ids.txt" <(sort $UNREL_FILE) > imp_3rd_ids.txt 

comm -12 imp_3rd_ids.txt <(sort "$CASES_FILE")  > case_3rd_ids.txt
comm -12 imp_3rd_ids.txt <(sort "$CONTROLS_FILE")  > control_3rd_ids.txt
cat case_3rd_ids.txt control_3rd_ids.txt | sort -u > case_control_3rd_ids.txt

### ANCESTRY SUMMARIES ###
awk '{print $1 "\t" "Case"}' case_3rd_ids.txt > sample_status.tmp
awk '{print $1 "\t" "Control"}' control_3rd_ids.txt >> sample_status.tmp

# breast case/controls
awk -F'\t' 'BEGIN {OFS="\t"} NR==FNR {status[$1] = $2; next} 
    FNR==1 {print $0, "Status"} FNR>1 && $1 in status {print $0, status[$1]}' \
    sample_status.tmp "$INPUT_FILE" > ancestry_breast.txt

# all 3rd unrelated samples 
awk -F'\t' 'NR==FNR {unrelated[$1] = 1; next} 
    FNR==1 {print} FNR>1 && $1 in unrelated {print}' \
    imp_3rd_ids.txt "$INPUT_FILE" > ancestry_3rd.txt

rm sample_status.tmp

### ANCESTRY-SPECIFIC IDS ###
unique_ancestries=$(cut -f3 ancestry_breast.txt | tail -n +2 | sort -u)

for ancestry in $unique_ancestries; do
    # breast case/controls
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && $NF=="Case" {print $1}' ancestry_breast.txt > "${ancestry}_case_ids.txt"
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && $NF=="Control" {print $1}' ancestry_breast.txt > "${ancestry}_control_ids.txt"
    
    # breast combined case/control
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && ($NF=="Case" || $NF=="Control") {print $1}' ancestry_breast.txt > "${ancestry}_case_control_ids.txt"
    
    # 3rd degree unrelated
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc {print $1}' ancestry_3rd.txt > "${ancestry}_3rd_ids.txt"
    
    # no filtering
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc {print $1}' "$INPUT_FILE" > "${ancestry}_ids.txt"
done

### IDS FOR --KEEP ### 
for file in *_ids.txt; do
    if [[ -f "$file" ]]; then
        pid_file="${file/_ids.txt/_pids.txt}"
        awk '{print "0\t" $1}' "$file" > "$pid_file"
    fi
done

cp *_ids.txt *_pids.txt $ID_DIR
