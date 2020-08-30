#!/bin/bash -l
use PLINK2
use Anaconda
python /seq/epiprod/ENCODE/genotyping/picardtesting/fingerprint_maps/build_fingerprint_maps.py \
    --recomb_file=/plink.chrX.GRCh38.map.shapeit.txt \
    --chrom="X" \
    --VCF_file="ALL.chrX_GRCh38.genotypes.20170504.vcf" \
    --LD_script="ldsc/ldsc.py" \
    --int_dir="intermediates_hg38" \
    --out_dir="output_75_20_hg38" \
    --prune_cutoff=0.20 \
    --clump_cutoff=0.75
