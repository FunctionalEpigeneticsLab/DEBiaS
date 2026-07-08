# A command-line tool for masking variant loci in single-cell or (pseudo)bulk bisulfite sequencing data
**DEBias (Detection and Elimination of Bias related to SNPs in bisulfite Sequencing)**

*Author: Lab for Functional Epigenetics, KU Leuven*

---

*Recommended workflow for handeling variants:*  
<img width="755" height="348.5" alt="image" src="https://github.com/user-attachments/assets/6fd2e615-066b-41f1-a702-a567b4c5583f" />

---

## Functions
Filter variant/blacklisted loci for CG or non-CG (CH) loci from DNA methylation calls. Masking happens when:
- CG: C itself is a variant, or the adjacent position in the corresponding CpG dinucleotide is a variant.
- CH: 1) the cytosine base itself is polymorphic/indel-affected, or 2) the immediate sequence context (strand '+': position i+1; strand '-': position i-1) shows CH to CG conversion.

## Prerequisites
- Python >=3.9
- Install pybedtools and pysam in the same environment: `pip install pybedtools && pip install pysam`
- Include [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) in your path: `export PATH=/path/to/dir_for_bedtools:$PATH`

## Installation
```bash
git clone https://github.com/FunctionalEpigeneticsLab/DEBiaS.git
cd DEBiaS
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

**`--has_header`** (default: False): the allc_file in the input contains header line or not.

**`--cpg_bedfile`**: a 0-based BED file containing all the CpG coordinates in the reference genome.  
> To generate such a BED file from any genome FASTA, we provide a python [script](scripts/extract_cpg_bed.py). 
> - Prerequisite: BioPython installed (`pip install biopython`)  
> - Example usage:  
	`python extract_cpg_bed.py -i genome.fa -o CpGsites.bed --threads 10`
> - Resources: It takes ~8G, 2 min for the human genome.

**`--CH`** (default: False): enable CH-only masking mode.

**`--ch_bedfile`**: a precomputed 0-based CH BED (4 columns: `chrom start end strand`) generated from reference FASTA using [scripts/extract_ch_bed.py](scripts/extract_ch_bed.py). In CH mode, DEBiaS masks a CH cytosine call if:
- the CH cytosine position is variant/indel,
- or the immediate context base is variant/indel (forward context `+1`, reverse context `-1`).

Generate once and reuse across runs:
`python scripts/extract_ch_bed.py -i genome.fa -o CHsites.bed`

**`--blacklist_provided`** (default: False): When this tag is enabled, the 2nd column in the input file should be a blacklist in BED format (with chromosome and position as the first two columns, and the third column set to position + 1). This blacklist will be used directly as the loci to be masked, without further matching to a `cpg_bedfile`, and no additional blacklist will be generated. Therefore, in this scenario, `cpg_bedfile` is not required.

### Examples
**1. Filter methylation calls:**
```bash
debias filter_mc \
	--input samplesheet.allc.csv \
	--cpg_bedfile CpGsites.bed \
	--output_dir debias_outs \
	--ncpus $SLURM_CPUS_PER_TASK
```

`samplesheet.allc.csv`
```csv
sample,vcf_file,allc_file
test_sample1,/path/to/test_input_files/donor1.filtered.vcf.gz,/path/to/test_input_files/cell1.allc.tsv.gz
test_sample2,/path/to/test_input_files/donor2.filtered.vcf.gz,/path/to/test_input_files/cell2.allc.tsv.gz
```

**2. Filter BAM:**
```bash
debias filter_bam \
	--input samplesheet.bam.csv \
	--cpg_bedfile CpGsites.bed \
	--output_dir debias_outs \
	--ncpus $SLURM_CPUS_PER_TASK
```

`samplesheet.bam.csv`
```csv
sample,vcf_file,bam_file
test_sample1,/path/to/test_input_files/donor1.filtered.vcf.gz,/path/to/test_input_files/cell1.deduplicated.filtered.sorted.bam
test_sample2,/path/to/test_input_files/donor2.filtered.vcf.gz,/path/to/test_input_files/cell2.deduplicated.filtered.sorted.bam
```

3. Generate blacklist BED files only:
```bash
debias generate_blacklist \
	--input samplesheet.blacklist.csv \
	--cpg_bedfile CpGsites.bed \
	--output_dir debias_outs \
	--ncpus $SLURM_CPUS_PER_TASK
```

`samplesheet.blacklist.csv`
```csv
sample,vcf_file
test_sample1,/path/to/test_input_files/donor1.filtered.vcf.gz
test_sample2,/path/to/test_input_files/donor2.filtered.vcf.gz
```

4. Filter methylation calls when blacklist BED files are provided:
```bash
debias filter_mc --blacklist_provided \
		--input samplesheet.allc_bl.csv \
		--output_dir allc_filtered \
		--ncpus $SLURM_CPUS_PER_TASK
``` 

`samplesheet.allc_bl.csv`
```csv
sample,vcf_file,allc_file
test_sample1,/path/to/test_input_files/donor1.CpG_variants_blacklist.bed,/path/to/test_input_files/cell1.allc.tsv.gz
test_sample2,/path/to/test_input_files/donor2.CpG_variants_blacklist.bed,/path/to/test_input_files/cell2.allc.tsv.gz
```

5. Filter methylation calls in CH mode:
```bash
debias filter_mc \
	--CH \
	--input samplesheet.allc.csv \
	--ch_bedfile CHsites.bed \
	--output_dir debias_outs \
	--ncpus $SLURM_CPUS_PER_TASK
```

6. Generate CH blacklist BED files only:
```bash
debias generate_blacklist \
	--CH \
	--input samplesheet.blacklist.csv \
	--ch_bedfile CHsites.bed \
	--output_dir debias_outs \
	--ncpus $SLURM_CPUS_PER_TASK
```

## Scalability for CG ALLC masking
<img width="523.5" height="371" alt="image" src="https://github.com/user-attachments/assets/2ca593e8-2298-4fb8-a7c6-9002af29f642" />

