import os
import sys
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from debias.utils import *

def process_sample(cpg_bed, output_dir, has_header, blacklist_provided, row):

    sample_name,vcf_file,allc_file = row
    output_allc_file = f"{sample_name}.var_masked.allc.tsv.gz"

    if not Path(vcf_file).exists() or not Path(allc_file).exists():
        raise FileNotFoundError(f"Missing: {' '.join(str(f) for f in [vcf, allc] if not f.exists())}")

    if not blacklist_provided:
        blacklist_bedfile = find_overlapping_positions(vcf_file, cpg_bed, output_dir, sample_name)

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

    if blacklist_provided:
        cpg_bedfile=None
    else:
        if cpg_bedfile is None:
            raise FileNotFoundError(f"Missing cpg_bedfile")

    # run in parallel
    print('Start processing samples...')
    process_func = partial(process_sample, cpg_bedfile, output_dir, has_header, blacklist_provided)

    with Pool(ncpus) as pool:
        results = pool.map(process_func, rows)

def process_bam_sample(cpg_bed, output_dir, blacklist_provided, row):

    sample_name,vcf_file,bam_file = row
    output_bam_file = f"{sample_name}.var_masked.bam"
    if not Path(vcf_file).exists() or not Path(bam_file).exists():
        raise FileNotFoundError(f"Missing: {' '.join(str(f) for f in [vcf, allc] if not f.exists())}")

    if not blacklist_provided:
        blacklist_bedfile = find_overlapping_positions(vcf_file, cpg_bed, output_dir, sample_name)

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

    if blacklist_provided:
        cpg_bedfile=None
    else:
        if cpg_bedfile is None:
            raise FileNotFoundError(f"Missing cpg_bedfile")

    print('Start processing samples...')
    process_func = partial(process_bam_sample, cpg_bedfile, output_dir, blacklist_provided)

    with Pool(ncpus) as pool:
        results = pool.map(process_func, rows)

# TODO: add process to only produce the blacklists
