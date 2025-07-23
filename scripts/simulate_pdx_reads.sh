#!/bin/bash

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
echo "Script directory: $SCRIPT_DIR"
echo "Project directory: $PROJECT_DIR"

REF_DIR="${PROJECT_DIR}/tests/data/references"
OUT_DIR="${PROJECT_DIR}/tests/data/fastp_input/synthetic_pdx"

echo "Using human reference: $HUMAN_REF"
echo "Using mouse reference: $MOUSE_REF"

HUMAN_REF="${REF_DIR}/hg38_chr22.fa"
MOUSE_REF="${REF_DIR}/mm39_chr19.fa"
TOTAL_READS=50000           # Total read pairs per sample
READ_LENGTH=150

mkdir -p "$OUT_DIR"

simulate_sample() {
  SAMPLE_NAME=$1
  HUMAN_PERCENT=$2

  HUMAN_READS=$(echo "$TOTAL_READS * $HUMAN_PERCENT / 100" | bc | awk '{print int($1)}')
  MOUSE_READS=$(echo "$TOTAL_READS - $HUMAN_READS" | bc)

  R1_OUT="${OUT_DIR}/${SAMPLE_NAME}_R1_001.fastq.gz"
  R2_OUT="${OUT_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"

  TMP_HUMAN_R1="tmp_human_${SAMPLE_NAME}_R1.fq"
  TMP_HUMAN_R2="tmp_human_${SAMPLE_NAME}_R2.fq"
  TMP_MOUSE_R1="tmp_mouse_${SAMPLE_NAME}_R1.fq"
  TMP_MOUSE_R2="tmp_mouse_${SAMPLE_NAME}_R2.fq"

  echo "🧬 Simulating sample: $SAMPLE_NAME"
  echo "   - Human reads: $HUMAN_READS ($HUMAN_PERCENT%)"
  echo "   - Mouse reads: $MOUSE_READS ($((100-HUMAN_PERCENT))%)"

  wgsim -N $HUMAN_READS -1 $READ_LENGTH -2 $READ_LENGTH -e 0 $HUMAN_REF $TMP_HUMAN_R1 $TMP_HUMAN_R2
  wgsim -N $MOUSE_READS -1 $READ_LENGTH -2 $READ_LENGTH -e 0 $MOUSE_REF $TMP_MOUSE_R1 $TMP_MOUSE_R2

  cat $TMP_HUMAN_R1 $TMP_MOUSE_R1 | gzip > $R1_OUT
  cat $TMP_HUMAN_R2 $TMP_MOUSE_R2 | gzip > $R2_OUT

  rm $TMP_HUMAN_R1 $TMP_HUMAN_R2 $TMP_MOUSE_R1 $TMP_MOUSE_R2

  echo "✅ Created: $R1_OUT and $R2_OUT"
  echo
}

# Create pdx90 (90% human reads) and pdx70 (70% human reads)
#simulate_sample "pdx90_S1" 90
#simulate_sample "pdx70_S2" 70
simulate_sample "human_S3" 100
