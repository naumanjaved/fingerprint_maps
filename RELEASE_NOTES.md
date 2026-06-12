## fingerprint_maps v2.0.0 — Python 3 modernization

First major update since the paper. The pipeline now runs on Python 3 with
modern, actively maintained tooling, while the originally published maps remain
unchanged in `map_files/`.

### Highlights
- Full Python 2 → Python 3 port of the map-building pipeline
- VCF filtering now uses **bcftools** (supports VCF v4.3 / GRCh38 releases)
- LD score step uses the **Python 3 LDSC port** (CBIIT/ldsc), run with the
  active interpreter — no separate environment needed
- Optional **parallel processing** across chromosomes via `--jobs`
- Faster I/O: gzip module for LD score files, set-based lookups, single-pass
  reads, batched writes
- Clearer failures (e.g. explicit error when clumping yields no output)

### Breaking changes
- Requires **Python 3.9** and a single conda env (see README)
- **bcftools** replaces **vcftools**
- LDSC must be the Python 3 port (CBIIT/ldsc), not the original Python 2 repo

### Project
- Added MIT **LICENSE**
- Added `requirements.txt` (PyVCF3 + pinned LDSC deps: numpy/pandas/scipy,
  bitarray<3)
- Added `download_vcfs.sh` for the 2019 GRCh38 1000 Genomes release (aria2c)
- Rewrote the README: single-env install, hg19/hg38 download + usage examples,
  `--jobs` guidance

### Validation
- Rebuilding chr22 from the original phase3 b37 inputs reproduces the published
  hg19 map almost exactly (795/798 shared positions), confirming the port is
  faithful.

### Notes
- The precomputed `hg38_*.map` files were built from phase3 SNPs lifted over
  from hg19. Running the current pipeline on the native 2019 GRCh38 callset
  produces a valid map that only partially overlaps them — this is expected and
  due to the input data, not the code.
