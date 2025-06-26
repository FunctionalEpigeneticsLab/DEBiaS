A tool to filter variant loci for (pseudo)bulk bisulfite sequencing methylation calling data

## Function
- Filter variant loci based on BED file and get the adjacent positions
- Remove blacklisted loci in ALLC/COV files (tab delimited, the first two columns should be chr and pos) using matching SNP positions


## installation
```bash
cd <repository_directory>
pip install -e .
```

## Usage

Command line:

```bash
debias --input samplesheet.csv --cpg_bedfile Human_hg38.CpGsites.0based.sorted.bed --output_dir test_out --mc_file_type cov --ncpus $SLURM_CPUS_PER_TASK
```
