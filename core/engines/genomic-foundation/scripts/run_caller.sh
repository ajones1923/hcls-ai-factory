#!/usr/bin/env bash
# E1: GPU variant-caller option for the Genomics Engine. Parabricks supports both;
# DeepVariant is the default (accuracy), HaplotypeCaller is the GATK-compatible option.
set -euo pipefail
CALLER="${CALLER:-deepvariant}"          # deepvariant | haplotypecaller
REF="${REF:-/ref/GRCh38.fa}"; BAM="${BAM:-/out/HG002.bam}"
OUT="${OUT:-/out/HG002.vcf.gz}"; NUM_GPUS="${NUM_GPUS:-1}"
case "$CALLER" in
  deepvariant)
    pbrun deepvariant     --ref "$REF" --in-bam "$BAM" --out-variants "$OUT" --num-gpus "$NUM_GPUS" ;;
  haplotypecaller)
    pbrun haplotypecaller --ref "$REF" --in-bam "$BAM" --out-vcf      "$OUT" --num-gpus "$NUM_GPUS" ;;
  *) echo "unknown CALLER='$CALLER' (use deepvariant|haplotypecaller)" >&2; exit 2 ;;
esac
