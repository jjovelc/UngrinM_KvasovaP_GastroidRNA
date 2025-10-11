#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --mem=16G

####### Set environment variables ###############

start_time=$(date +%s)
echo "Job started on: $(date)"

####### Run your script #########################

DIR="/work/vetmed_shared_dbs/transcriptomes"
TRANSCRIPTOME_IDX="Mus_musculus.GRCm39.cdna.all_idx"  # Transcriptome index directory
READ1=$1                          # Read 1 file (for paired-end)
READ2=$2                          # Read 2 file (for paired-end)
PREFIX="${READ1/_R1.fastq.gz*/}"
OUTPUT_DIR="${PREFIX}_salmon_quant_output"               # Output directory for this sample

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

eval "$(conda shell.bash hook)"
conda activate salmon

# For paired-end bulk RNA-seq
salmon quant \
    -i $DIR/$TRANSCRIPTOME_IDX \
    -l A \
    -1 $READ1 \
    -2 $READ2 \
    -o $OUTPUT_DIR \
    --validateMappings \
    --gcBias \
    --seqBias \
    -p 24

end_time=$(date +%s)
echo "Job ended on: $(date)"

# Calculate the difference in times to get the duration
duration=$((end_time - start_time))

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Job duration: $hours hours, $minutes minutes, $seconds seconds"
