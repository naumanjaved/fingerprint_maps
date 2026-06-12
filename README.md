# build_fingerprint_maps

`build_fingerprint_maps` is a tool for building haplotype maps for use with
[Picard](http://broadinstitute.github.io/picard/) fingerprinting software
`CrosscheckFingerprints`. A haplotype map is a collection of "blocks" of SNPs
which are in tight linkage with SNPs of the same block and low linkage with
SNPs of different blocks.

## Precomputed map files

Precomputed maps with headers (headers do not contain entries for scaffolds or contigs):

* [hg19 with "chr" prefix](map_files/hg19_chr.map)
* [hg19 without "chr" prefix](map_files/hg19_nochr.map)
* [hg38 with "chr" prefix](map_files/hg38_chr.map)
* [hg38 without "chr" prefix](map_files/hg38_nochr.map)

The `map_files/` directory also contains pre-computed maps with relaxed intra-
and inter-block correlation thresholds. Map names contain the parameters used
(e.g. `hg19_75_20_chr.map` uses `--clump_cutoff=0.75 --prune_cutoff=0.20`).

These precomputed maps were generated with the prior python 2 implementation,
the python 2 version of LDSC, and an earlier 1000 Genomes release (hg19, with
the hg38 maps built from SNPs lifted over from hg19). They are kept because
they are the exact maps used in the paper. Maps built with the current python 3
pipeline on the native GRCh38 callset will differ from the published hg38 maps
(see the note under Usage).

## Installation

Everything runs in a single dedicated conda environment. LDSC's python3 port
pins old package versions (numpy 1.21, pandas 1.3, scipy 1.7, bitarray <3) that
require python 3.9, so the env is built on python 3.9:

```bash
git clone https://github.com/naumanjaved/fingerprint_maps.git
cd fingerprint_maps

# create the dedicated env and install python deps
conda create -n fpmaps -c conda-forge python=3.9
conda activate fpmaps
pip install -r requirements.txt
```

PLINK and bcftools are command-line tools, not python packages. Either add them
to the env or make sure they are on your `PATH`:

```bash
conda install -n fpmaps -c bioconda plink bcftools
```

LDSC itself is a script, not a package. Clone the python3 port and point
`--LD_script` at its `ldsc.py`:

```bash
git clone -b ldsc39 https://github.com/CBIIT/ldsc.git
```

IMPORTANT: always `conda activate fpmaps` before running the pipeline.
`build_fingerprint_maps.py` runs `ldsc.py` with the same interpreter it is
launched with, so the active env (and its pinned bitarray <3) is what LDSC uses.

## Dependencies

1. Python 3.9 (pinned by LDSC's dependencies)
2. [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
3. [bcftools](https://samtools.github.io/bcftools/)
4. [LDSC (LD Score regression, python3 port)](https://github.com/CBIIT/ldsc)
5. python packages in `requirements.txt` (PyVCF3 plus LDSC's numpy/pandas/scipy/bitarray)

## Required data files

**Recombination maps** are bundled in this repo and need to be extracted:

```bash
cd references
tar xzf genetic_map_b37.tar.gz    # hg19
tar xzf genetic_map_hg38.tar.gz   # hg38
```

**1000 Genomes VCFs** (one per chromosome) must be downloaded separately.

Install aria2c first (available on conda-forge, no compiler or perl needed):

```bash
conda install -c conda-forge aria2
```

For **GRCh38 (hg38)**, a download script is included that fetches all
chromosomes (~15 GB total) from the
[1000 Genomes 2019 GRCh38 release](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)
(biallelic SNV/INDEL, phased) using [aria2c](https://aria2.github.io/):

```bash
./download_vcfs.sh          # downloads to vcfs/ by default
./download_vcfs.sh mydir    # or specify a custom directory
```

For **GRCh37 (hg19)**, download the phase3 (20130502) b37 callset instead, e.g.
for a single chromosome:

```bash
aria2c -x 8 -d vcfs \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

Both compressed (`.vcf.gz`) and uncompressed VCF inputs are supported. No
separate reference genome FASTA is needed.

## Usage

See `sample_run.sh` for a sample run script, and run
`python3 build_fingerprint_maps.py -h` for the full list of options. Always
`conda activate fpmaps` first.

**hg38 (GRCh38) single chromosome:**

```bash
python3 build_fingerprint_maps.py \
    --chrom 22 \
    --VCF_file="vcfs/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" \
    --recomb_file="references/genetic_map_hg38/plink.chr22.GRCh38.map.shapeit.txt" \
    --LD_script="ldsc/ldsc.py" \
    --int_dir="intermediates" \
    --out_dir="output"
```

**hg19 (GRCh37) single chromosome:**

```bash
python3 build_fingerprint_maps.py \
    --chrom 22 \
    --VCF_file="vcfs/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
    --recomb_file="references/genetic_map_b37/genetic_map_chr22_combined_b37.txt" \
    --LD_script="ldsc/ldsc.py" \
    --int_dir="intermediates" \
    --out_dir="output"
```

To process multiple chromosomes in parallel, use `{chrom}` placeholders in the
VCF and recombination file paths and pass `--jobs`:

```bash
python3 build_fingerprint_maps.py \
    --chrom all \
    --VCF_file="vcfs/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" \
    --recomb_file="references/genetic_map_hg38/plink.chr{chrom}.GRCh38.map.shapeit.txt" \
    --LD_script="ldsc/ldsc.py" \
    --int_dir="intermediates" \
    --out_dir="output" \
    --jobs 4
```

`--jobs N` processes N chromosomes at once. PLINK and LDSC each use multiple
threads and a fair bit of memory per chromosome, so a large `N` can
oversubscribe CPU/RAM (this only affects speed, not correctness). Pick `N`
based on your machine's available cores and memory.

### Note on hg38 results vs the published maps

The published `hg38_*.map` files in `map_files/` were generated from the
original phase3 SNPs lifted over from hg19 to hg38. The `download_vcfs.sh`
script above instead fetches the native 2019 GRCh38 callset, which was
re-genotyped and re-phased and therefore contains a different SNP set. Running
on these native GRCh38 VCFs will produce a valid map that only partially
overlaps the published hg38 maps (roughly half the SNP positions in practice).
This is expected and is due to the input data, not the pipeline: running on the
original phase3 b37 VCFs reproduces the published hg19 map almost exactly.

## Use with Picard

The above maps are to be used as inputs to `CrosscheckFingerprints` and
`ExtractFingerprint` in the Picard suite.

For most cases where each file you want to compare with
`CrosscheckFingerprints` contains data for only a single fingerprint, you
should run with the `CROSSCHECK_BY=FILE` flag. Picard with default settings can
be strict about properly formatted headers and read names, so if a validation
error arises, try running with `VALIDATION_STRINGENCY=LENIENT` (after ensuring
that the formatting error does not indicate a legitimate problem with the input
bam file).

When comparing many files, it is recommended to precompute VCFs containing
extracted fingerprints using `ExtractFingerprint`. This avoids
`CrosscheckFingerprints` having to redundantly compute fingerprints for the
same file each time it is used for a comparison.

## Custom map files

If you create a custom map file, make sure to prepend the appropriate header
file. Below are headers for hg19 and hg38 with entries for reference
chromosomes:

* [hg19 with "chr" prefix](headers/header_hg19_chr)
* [hg19 without "chr" prefix](headers/header_hg19_nochr)
* [hg38 with "chr" prefix](headers/header_hg38_chr)
* [hg38 without "chr" prefix](headers/header_hg38_nochr)

## Support

Open an issue or email nauman.javed48@gmail.com for support.

## Authors

* Nauman Javed (Broad Institute) -- wrote the map generation scripts
* Yossi Farjoun (Broad Institute) -- wrote `CrosscheckFingerprints` and `ExtractFingerprints`

## Citation

If you use the above tool/maps with `CrosscheckFingerprints` in your
publication, please cite the Picard-tools repo as well as:

> Javed, N., Farjoun, Y., Fennell, T.J. et al. Detecting sample swaps in
> diverse NGS data types using linkage disequilibrium. *Nat Commun* 11, 3697
> (2020). [doi:10.1038/s41467-020-17453-5](https://doi.org/10.1038/s41467-020-17453-5)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file
for details.
