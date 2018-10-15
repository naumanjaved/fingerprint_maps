
# build_fingerprint_maps 
`build_fingerprint_maps` is a tool for building haplotype maps for use with Picardtools fingerprinting software. A haplotype map is a collection of "blocks" of SNPs which are in tight linkage with SNPs of the same block and low linkage with SNPs of different blocks.

## Dependencies

In order to download `build_fingerprint_maps`, you should clone this repository via the command
```  
git clone https://github.com/naumanjaved/fingerprint_maps.git
```
In order to run `build_fingerprint_maps`, you must have working installations of:

1. Python (>=2.7)
1. PLINK2 - https://www.cog-genomics.org/plink2
2. VCFTools - http://vcftools.sourceforge.net/man_latest.html
3. Anaconda(https://anaconda.org/anaconda/python) or the following modules:
     a. `subprocess` 
     b. `os` 
     c. `itertools`
     d. `numpy`
     e. `sys` 
     f. `argparse`
     g. `traceback`
     h. `time`
     i. `datetime` 
4. LDSC - https://github.com/bulik/ldsc

## Required Files
1. 1000 Genomes Phase 3 VCFs(hg19) - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
2. 1000 Genomes Phase 3 Recombination maps(hg19) - http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/

## How to run
See run.sh to see a sample run script.
Run python build_fingerprint_maps.py -h to see a list of command line options.

# Use with Picardtools
Before using with Picardtools, append the file called "header" to the beginning of a map file.

## Support

Email javed@broadinstitute.org for issues.

## Authors

Nauman Javed(Broad Institute)


