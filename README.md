
# build_fingerprint_maps 
`build_fingerprint_maps` is a tool for building haplotype maps for use with Picardtools fingerprinting software. A haplotype map is a collection of "blocks" of SNPs which are in tight linkage with SNPs of the same block and low linkage with SNPs of different blocks.

In order to download `build_fingerprint_maps`, you should clone this repository via the command
```
git clone https://github.com/naumanjaved/fingerprint_maps.git
```

## Precomputed map files
* [hg19 with "chr" prefix before chromosome names](https://github.com/naumanjaved/fingerprint_maps/blob/master/map_files/hg19_chr.map)
* [hg19 without "chr" prefix before chromosome names](https://github.com/naumanjaved/fingerprint_maps/blob/master/map_files/hg19_nochr.map)
* [hg38 with "chr" prefix before chromosome names](https://github.com/naumanjaved/fingerprint_maps/blob/master/map_files/hg38_chr.map)
* [hg38 without "chr" prefix before chromosome names](https://github.com/naumanjaved/fingerprint_maps/blob/master/map_files/hg38_nochr.map)


## Dependencies

In order to run `build_fingerprint_maps`, you must have working installations of:

1. Python (>=2.7)
1. [PLINK2](https://www.cog-genomics.org/plink2)
2. [VCFTools](http://vcftools.sourceforge.net/man_latest.html)
3. [Anaconda](https://anaconda.org/anaconda/python) or the following modules:
     a. `subprocess` 
     b. `os` 
     c. `itertools`
     d. `numpy`
     e. `sys` 
     f. `argparse`
     g. `traceback`
     h. `time`
     i. `datetime` 
4. [LDSC(LDScore regression)](https://github.com/bulik/ldsc)

## Required Files
Fingerprint maps uses VCFs from 1000 Genomes Phase 3 and recombination maps(SHAPEIT format). These can be found here: 
* [hg19 VCFs](ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) 
* [hg38 VCFs(liftover)](ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/)
* [hg19 recombination maps from SHAPEIT](references/genetic_map_b37.tar.gz) 
* [hg38 recombination maps from SHAPEIT liftover](references/genetic_map_hg38.tar.gz)

See run.sh to see a sample run script.
Run python build_fingerprint_maps.py -h to see a list of command line options.

# Use with Picardtools
Before using with Picardtools, append the appropriate header file to the beginning of the map file:
* [hg19 with "chr" prefix](https://github.com/naumanjaved/fingerprint_maps/blob/master/headers/header_hg19_chr)
* [hg19 without "chr" prefix](https://github.com/naumanjaved/fingerprint_maps/blob/master/headers/header_hg19_nochr)
* [hg38 with "chr" prefix](https://github.com/naumanjaved/fingerprint_maps/blob/master/headers/header_hg38_chr)
* [hg38 without "chr" prefix](https://github.com/naumanjaved/fingerprint_maps/blob/master/headers/header_hg38_nochr)

## Support

Email javed@broadinstitute.org for issues.

## Authors

Nauman Javed(Broad Institute)


