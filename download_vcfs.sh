#!/bin/bash
#
# download 1000 genomes phase 3 VCFs (GRCh38, biallelic SNV/INDEL) into a vcfs/ directory.
# requires aria2c: apt install aria2  or  conda install -c conda-forge aria2
# total download is approximately 15 GB.

set -e

OUTDIR="${1:-vcfs}"
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"

mkdir -p "$OUTDIR"

URLS_FILE=$(mktemp)
trap "rm -f $URLS_FILE" EXIT

for chr in $(seq 1 22) X; do
    echo "${BASE_URL}/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    echo "${BASE_URL}/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi"
done > "$URLS_FILE"

aria2c \
    --input-file="$URLS_FILE" \
    --dir="$OUTDIR" \
    --max-concurrent-downloads=4 \
    --max-connection-per-server=4 \
    --continue=true \
    --auto-file-renaming=false

echo "done. files saved to ${OUTDIR}/"
