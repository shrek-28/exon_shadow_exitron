#!/usr/bin/env bash
set -euo pipefail

cd "/c/Users/admin/Desktop/Coding/exon_shadow_exitron/shadow_sequence_content"

PYTHON_PATH="/c/Users/admin/AppData/Local/Programs/Python/Python312/python.exe"
PERL_SCRIPT_PATH="/c/Users/admin/Desktop/Coding/exon_shadow_exitron/shadow_sequence_content/GC_Stretch_finder.pl"

for species in bonobo borangutan gibbon gorilla sorangutan human chimpanzee
do
  for region in upstream downstream
  do
    input_path="shadow_fasta/${region}_${species}_filtered_combined_non_zero.fa"
    output_path="output_files/${species}_${region}.csv"

    "$PYTHON_PATH" seq_processing_with_stretch.py \
      --fasta "$input_path" \
      --out "$output_path" \
      --perl "$PERL_SCRIPT_PATH"

  done
done

echo "All jobs finished."
