#!/bin/bash

# Usage: ./shadow_summary_species.sh consolidated_output.tsv

OUTPUT="$1"

# List of species names
species_list=("bonobo" "gorilla" "chimpanzee" "orangutan" "gibbon" "macaque" "baboon")

# Base folder path where the files are located
BASE_PATH="/path/to/folder"   # adjust this to your folder

# Temporary file to store results
TMP=$(mktemp)

for sp in "${species_list[@]}"; do
    FILE="$BASE_PATH/${sp}.tsv"  # assuming file names like bonobo.tsv, gorilla.tsv, etc.

    if [[ ! -f "$FILE" ]]; then
        echo "File not found: $FILE"
        continue
    fi

    TOTAL=$(awk -F'\t' 'NR>1 {print $15; print $16}' "$FILE" | wc -l)

    awk -F'\t' -v total="$TOTAL" -v species="$sp" 'NR>1 {
        if($15>0){pos++; up++} else if($15<0) neg++; else zero++
        if($16>0){pos++; down++} else if($16<0) neg++; else zero++
    }
    END {
        if(pos>0){up_perc = up/pos*100; down_perc = down/pos*100} else {up_perc=0; down_perc=0}
        printf "%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%.2f\t%.2f\n", species, total, pos, pos/total*100, neg, neg/total*100, zero, zero/total*100, up_perc, down_perc
    }' "$FILE" >> "$TMP"

done

# Add header and save to output file
sed '1iSpecies\tTotal\tPositive\tPositive%\tNegative\tNegative%\tZero\tZero%\tPos_upstream%\tPos_downstream%' "$TMP" > "$OUTPUT"

rm "$TMP"

echo "Consolidation done. Output saved in $OUTPUT"
