# A command-line tool for filtering variant loci in single cell or (pseudo)bulk bisulfite sequencing data

> *Recommended workflow for handeling variants:*  
> single-cell BAM files -> pseudobulk BAM files -> split by chromosomes -> SNP calling -> SNP filtering -> joining chromosomes -> `DEBiaS`: extract adjacent positions within the same CpG -> SNP masking for pseudobulk or single-cell BAM files / SNP masking for ALLC files -> downstream analysis

## Functions
- Filter variant loci based on BED file to obtain a comprehensive "blacklist" including the adjacent positions in the corresponding CpG dinucleotide.
- Remove blacklisted loci in methylation calling ALLC/COV files or BAM files, which can come from either single cells or (pseudo)bulk samples.

## Prerequisites
- Python >=3.9
- Install pybedtools and pysam in the same environment: `pip install pybedtools && pip install pysam`
- Include [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) in your path: `export PATH=/path/to/dir_for_bedtools:$PATH`

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
- `filter_bam`: Filter BAM files.
- `generate_blacklist`: Generate blacklist BED files from VCF files without performing filtering.

### Parameters

**`--input`**: a CSV file with columns depending on the subcommand:
- For `filter_mc` and `filter_bam`: 3 columns required:
  1. **sample**: Sample name, used as the output file prefix.
  2. **vcf_file**: Path to a `.vcf` or `.vcf.gz` file, pre-filtered to include only variants with non-reference genotypes (`0/1` or `1/1`).
  3. **allc_file**/**bam_file**:  
     - For `filter_mc`: Path to a `.allc.tsv`, `.allc.tsv.gz` (Methylpy/AllCools), or `.cov`, `.cov.gz` (Bismark) file. These must be tab-delimited with chromosome and cytosine position as the first two columns.
     - For `filter_bam`: Path to a `.bam` file.
- For `generate_blacklist`: 2 columns required:
  1. **sample**: Sample name, used as the output file prefix.
  2. **vcf_file**: Path to a `.vcf` or `.vcf.gz` file, pre-filtered to include only variants with non-reference genotypes (`0/1` or `1/1`).

**`--has_header`**: the allc_file contains header line or not (default: False)

**`--cpg_bedfile`**: a 0-based BED file containing all the CpG coordinates in the reference genome.  
- To generate such a BED file from any genome FASTA, we provide a python [script](scripts/extract_cpg_bed.py). 
- Prerequisite: BioPython installed (`pip install biopython`)  
- Example usage:  
	```bash
	python extract_cpg_bed.py -i genome.fa -o Human_hg38.CpGsites.0based.sorted.bed --threads 10
	```
- Resources: It takes ~8G, 2min for the human genome.

**`--blacklist_provided`**: When enabled, the second column in the input file should be a blacklist in BED format (with chromosome and position as the first two columns, and the third column set to position + 1). This blacklist will be used directly, without further matching to a `cpg_bedfile`, and no additional blacklist will be generated. In this scenario, `cpg_bedfile` is not required.

### Examples
1. Filter methylation calls:
```bash
debias filter_mc --input samplesheet.allc.csv --cpg_bedfile Human_hg38.CpGsites.0based.sorted.bed --output_dir debias_outs --ncpus $SLURM_CPUS_PER_TASK
```

`samplesheet.allc.csv`
```csv
sample,vcf_file,allc_file
test_sample1,/path/to/test_input_files/donor1.filtered.vcf.gz,/path/to/test_input_files/cell1.allc.tsv.gz
test_sample2,/path/to/test_input_files/donor2.filtered.vcf.gz,/path/to/test_input_files/cell2.allc.tsv.gz
```

2. Filter BAM:
```bash
debias filter_bam --input samplesheet.bam.csv --cpg_bedfile Human_hg38.CpGsites.0based.sorted.bed --output_dir debias_outs --ncpus $SLURM_CPUS_PER_TASK
```

`samplesheet.bam.csv`
```csv
sample,vcf_file,bam_file
test_sample1,/path/to/test_input_files/donor1.filtered.vcf.gz,/path/to/test_input_files/cell1.deduplicated.filtered.sorted.bam
test_sample2,/path/to/test_input_files/donor2.filtered.vcf.gz,/path/to/test_input_files/cell2.deduplicated.filtered.sorted.bam
```

3. Generate blacklist BED files only:
```bash
debias generate_blacklist --input samplesheet.blacklist.csv --cpg_bedfile Human_hg38.CpGsites.0based.sorted.bed --output_dir debias_outs --ncpus $SLURM_CPUS_PER_TASK
```

`samplesheet.blacklist.csv`
```csv
sample,vcf_file
test_sample1,/path/to/test_input_files/donor1.filtered.vcf.gz
test_sample2,/path/to/test_input_files/donor2.filtered.vcf.gz
```
