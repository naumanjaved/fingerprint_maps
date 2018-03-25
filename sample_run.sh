#!/bin/bash -l

use Anaconda 
#use PLINK2 
use VCFtools
python build_fingerprint_maps.py \
    --recomb_directory='../1000GP_Phase3/' \
    --chromosome='1' \
    --VCF_file='../chr1.vcf' \
    --LD_script='../ldsc/ldsc.py'
