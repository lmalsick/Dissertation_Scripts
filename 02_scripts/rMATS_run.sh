#!/usr/bin/env bash

# run_rmats.sh: wrapper to run rMATS-turbo for differential splicing
# Usage: bash rMATS_run.sh <group1_bam_list.txt> <group2_bam_list.txt> <annotation.gtf> <output_dir> [threads] [read_length]
# Example: bash run_rmats.sh wt_bams.txt mut816_bams.txt /path/to/annotation.gtf rmats_output_816_vs_744 20 144

set -euo pipefail

# Required inputs
group1_bams=$1   # text file: comma-separated or space-separated BAM paths for group 1 (e.g. WT)
group2_bams=$2   # text file: BAM list for group 2 (e.g. mutant)
gtf=$3          # genome annotation GTF file
output_dir=$4   # directory to store rMATS results

# Optional inputs	hreads=${5:-20}
read_length=${6:-150}

# Create output directory
mkdir -p "$output_dir"

echo "Running rMATS-turbo with:" 
 echo "  Group1 BAMs: $group1_bams"
 echo "  Group2 BAMs: $group2_bams"
 echo "  GTF:         $gtf"
 echo "  Output:      $output_dir"
 echo "  Threads:     $threads"
 echo "  Read length: $read_length"

# Concatenate BAM lists into comma-separated strings
b1=$(paste -sd "," "$group1_bams")
b2=$(paste -sd "," "$group2_bams")

tmp_dir="$output_dir/tmp"
mkdir -p "$tmp_dir"

# Run rMATS-turbo
rmats.py \
  --b1 "$b1" \
  --b2 "$b2" \
  --gtf "$gtf" \
  --od "$output_dir" \
  --tmp "$tmp_dir" \
  --readLength $read_length \
  --cstat 0.0001 \
  --task both \
  --libType fr-unstranded \
  --nthread $threads

# After completion
echo "rMATS-turbo analysis completed. Results in $output_dir"
