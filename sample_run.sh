#!/bin/bash -l

use Anaconda 
use PLINK2 
use VCFtools
python build_fingerprint_maps.py \
    --recomb_directory='../1000GP_Phase3/' \
    --chromosome='1' \
    --VCF_file='../VCFs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf' \
    --LD_script='../ldsc/ldsc.py'
