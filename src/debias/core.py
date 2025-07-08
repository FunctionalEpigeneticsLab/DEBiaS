import os
import sys
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from debias.utils import *

def write_blacklist_bed(vcf_file, cpg_chrom_pos_list, output_dir, sample_name):
    """Process VCF file: subset to CpG sites variants and include adjacent cytosine."""
    
    variant_loci = parse_vcf(vcf_file)
    print(f'Finished parsing {vcf_file}', flush=True)
    loci_blacklist = find_overlapping_positions(variant_loci, cpg_chrom_pos_list)
    del(variant_loci)
    blacklist_bedfile = os.path.join(output_dir, f"{sample_name}.CpG_variants_blacklist.bed")
    with open(blacklist_bedfile,'w') as f:
        for chrom,pos in loci_blacklist:
            f.write(f"{chrom}\t{pos}\t{pos+1}\n")
    return blacklist_bedfile

def process_sample(cpg_chrom_pos_list, output_dir, has_header, blacklist_provided, row):

    sample_name,vcf_file,allc_file = row
    output_allc_file = f"{sample_name}.var_masked.allc.tsv"

    if not Path(vcf_file).exists() or not Path(allc_file).exists():
        raise FileNotFoundError(f"Missing: {' '.join(str(f) for f in [vcf, allc] if not f.exists())}")

    if not blacklist_provided:
        blacklist_bedfile = write_blacklist_bed(vcf_file, cpg_chrom_pos_list, output_dir, sample_name)

        # remove those cytosines from the allc file
        print(f'Start processing {allc_file}', flush=True)
        mask_allc_with_bedtools(allc_file, has_header, blacklist_bedfile, output_dir, output_allc_file)
        print(f'Result written to {output_allc_file}', flush=True)

    else:
        print(f'Start processing {allc_file}', flush=True)
        mask_allc_with_bedtools(allc_file, has_header, vcf_file, output_dir, output_allc_file) # here the vcf file should be bed format
        print(f'Result written to {output_allc_file}', flush=True)

def run_parallel(samplesheet, cpg_bedfile, output_dir, ncpus, has_header, blacklist_provided):

    # parse the files that only needs to be read once
    os.makedirs(output_dir, exist_ok=True)
    rows = parse_samplesheet(samplesheet)

    if cpg_bedfile is not None:
        print('Parsing input BED file...')
        cpg_chrom_pos_list = parse_bed(cpg_bedfile)
    else:
        cpg_chrom_pos_list = None

    # run in parallel
    print('Start processing samples...')
    process_func = partial(process_sample, cpg_chrom_pos_list, output_dir, has_header, blacklist_provided)

    with Pool(ncpus) as pool:
        results = pool.map(process_func, rows)

def process_bam_sample(cpg_chrom_pos_list, output_dir, blacklist_provided, row):

    sample_name,vcf_file,bam_file = row
    output_bam_file = f"{sample_name}.var_masked.bam"
    if not Path(vcf_file).exists() or not Path(allc_file).exists():
        raise FileNotFoundError(f"Missing: {' '.join(str(f) for f in [vcf, allc] if not f.exists())}")

    if not blacklist_provided:
        blacklist_bedfile = write_blacklist_bed(vcf_file, cpg_chrom_pos_list, output_dir, sample_name)

        # remove those cytosines from the BAM file
        print(f'Start processing {bam_file}', flush=True)
        mask_bam(bam_file, blacklist_bedfile, output_dir, output_bam_file)
        print(f'Result written to {output_bam_file}', flush=True)

    else:
        print(f'Start processing {bam_file}', flush=True)
        mask_bam(bam_file, vcf_file, output_dir, output_bam_file)
        print(f'Result written to {output_bam_file}', flush=True)


def run_bam_parallel(samplesheet, cpg_bedfile, output_dir, ncpus, blacklist_provided):

    os.makedirs(output_dir, exist_ok=True)
    rows = parse_samplesheet(samplesheet)

    if cpg_bedfile is not None:
        print('Parsing input BED file...')
        cpg_chrom_pos_list = parse_bed(cpg_bedfile)
    else:
        cpg_chrom_pos_list = None

    print('Start processing samples...')
    process_func = partial(process_bam_sample, cpg_chrom_pos_list, output_dir, blacklist_provided)

    with Pool(ncpus) as pool:
        results = pool.map(process_func, rows)

# TODO: add process to only produce the blacklists
