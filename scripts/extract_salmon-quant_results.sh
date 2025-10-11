#!/usr/bin/bash

# Extract transcript names from first sample (header)
FIRST_DIR=$(ls -d *_salmon_quant_output | head -n 1)
cut -f 1 "${FIRST_DIR}/quant.sf" > names.txt

# Process each sample directory
for DIR in *_salmon_quant_output
do
    PREFIX=${DIR/_salmon_quant_output/}
    
    # Extract TPM values (column 4), skip header, add sample name as header
    tail -n +2 "${DIR}/quant.sf" | cut -f 4 > temp_tpm.txt
    echo "$PREFIX" | cat - temp_tpm.txt > "${PREFIX}_tpms.tsv"
    
    # Extract NumReads values (column 5), skip header, add sample name as header
    tail -n +2 "${DIR}/quant.sf" | cut -f 5 > temp_counts.txt
    echo "$PREFIX" | cat - temp_counts.txt > "${PREFIX}_counts.tsv"
    
    rm temp_tpm.txt temp_counts.txt
done

# Combine all samples
paste names.txt *_tpms.tsv > allSamples_tpms.tsv
paste names.txt *_counts.tsv > allSamples_counts.tsv

# Clean up individual sample files (optional)
# rm *_tpms.tsv *_counts.tsv names.txt

echo "Processing complete!"
echo "Output files:"
echo "  - allSamples_tpms.tsv"
echo "  - allSamples_counts.tsv"
