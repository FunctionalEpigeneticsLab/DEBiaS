# A command-line tool for filtering variant loci in single cell or (pseudo)bulk bisulfite sequencing data

## Function
- Filter variant loci based on BED file to obtain a comprehensive "blacklist" including the adjacent positions in the corresponding CpG dinucleotide.
- Remove blacklisted loci in methylation calling ALLC/COV files or BAM files, which can be from single cells or (pseudo)bulk samples.

## Prerequisites
- Python >=3.9
- Include [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) in your path: `export PATH=/path/to/dir_for_bedtools:$PATH`
- Install pybedtools and pysam in the same environment: `pip install pybedtools`; `pip install pysam`

## Installation
```bash
git clone <this_repository>
cd <path_to_this_repository>
pip install -e .
```

## Usage

### Subcommands
This tool supports the following subcommands:

- `filter_mc`: Filter methylation calls.
- `filter_bam`: Filter BAM alignments.

### Parameters

**`--input`**: a CSV file with 3 columns:  
- `sample` (1st column) is the sample name, used as prefix for the output. 
- `vcf_file` (2nd column) should be a `.vcf` or `.vcf.gz` file prefiltered to only include variants that are different from the reference genome bases (`0/1` or `1/1`). 
- `allc_file` (3rd column when calling `filter_mc`) should be a `.allc.tsv` or `allc.tsv.gz` (output from Methylpy or AllCools) or `.cov` or `.cov.gz` file (output from Bismark). These files should be tab delimited, and the first two columns should be chromosome and cytosine position.

- `bam_file`  (3rd column when calling `filter_bam`) should be a `.bam` file.

**`--has_header`**: the allc_file contains header line or not (default: False)

**`--cpg_bedfile`**: a 0-based BED file containing all the CpG coordinates matching the reference genome build used in the input. The BED files for GRCh38 and T2T-CHM13v2.0 are provided here in `reference`.

**`--blacklist_provided`**: when this option is enabled, the 2nd column provided in 'input' should be a blacklist (the first 2 column should be chromosome and position) that does not need further matching to a 'cpg_bedfile'. Enabling this option will skip generating a blacklist. Thus, 'cpg_bedfile' is not needed in this senario.


### Examples
Filter methylation calls:
```bash
debias filter_mc --input samplesheet.csv --cpg_bedfile Human_hg38.CpGsites.0based.sorted.bed --output_dir debias_outs --ncpus $SLURM_CPUS_PER_TASK
```

Filter BAM alignments:
```bash
debias filter_bam --input samplesheet.csv --cpg_bedfile Human_hg38.CpGsites.0based.sorted.bed --output_dir debias_outs --ncpus $SLURM_CPUS_PER_TASK
```
