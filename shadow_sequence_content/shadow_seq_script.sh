cd "/c/Users/admin/Desktop/Coding/exon_shadow_exitron/shadow_sequence_content"
PYTHON_PATH="/c/Users/admin/AppData/Local/Programs/Python/Python312/python.exe"
for species in bonobo borangutan gibbon gorilla sorangutan human chimpanzee
do
  for region in upstream downstream
  do
    input_path="shadow_fasta/${region}_${species}_filtered_combined_non_zero.fa"
    output_path="output_files/${species}_${region}.csv"
    "$PYTHON_PATH" sequence_processing.py \
      --fasta "$input_path" \
      --out "$output_path"
  done
done
echo "All jobs finished."