#!/bin/bash -l

# sample run script for building a fingerprint map for chromosome 22 (hg38)
# adjust paths below to match your environment
# run from inside the fpmaps conda env (conda activate fpmaps)

# single chromosome
python3 build_fingerprint_maps.py \
    --recomb_file="references/genetic_map_hg38/plink.chr22.GRCh38.map.shapeit.txt" \
    --chrom 22 \
    --VCF_file="vcfs/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" \
    --LD_script="ldsc/ldsc.py" \
    --int_dir="intermediates_hg38" \
    --out_dir="output_75_20_hg38" \
    --prune_cutoff=0.20 \
    --clump_cutoff=0.75

# all autosomes + X in parallel (4 jobs), using {chrom} placeholder in file paths
# python3 build_fingerprint_maps.py \
#     --recomb_file="references/genetic_map_hg38/plink.chr{chrom}.GRCh38.map.shapeit.txt" \
#     --chrom all \
#     --VCF_file="vcfs/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" \
#     --LD_script="ldsc/ldsc.py" \
#     --int_dir="intermediates_hg38" \
#     --out_dir="output_75_20_hg38" \
#     --prune_cutoff=0.20 \
#     --clump_cutoff=0.75 \
#     --jobs 4
