import os
import sys
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from collections import defaultdict
import pybedtools
from debias.utils import *

def process_sample(output_dir, has_header, blacklist_bedfile, row):
    """Process a single sample using a pre-generated or provided blacklist"""
    sample_name, vcf_file, allc_file = row
    output_allc_file = f"{sample_name}.var_masked.allc.tsv.gz"

    if not Path(allc_file).exists():
        raise FileNotFoundError(f"Missing ALLC file: {allc_file}")

    # remove those cytosines from the allc file
    print(f'Start processing {allc_file}', flush=True)
    mask_allc_with_bedtools(allc_file, has_header, blacklist_bedfile, output_dir, output_allc_file)
    print(f'Result written to {output_allc_file}', flush=True)

def generate_vcf_blacklists(rows, cpg_bedfile, output_dir):
    """
    Generate blacklists for each unique VCF file (only once per VCF).
    Returns a dict mapping vcf_file paths to their corresponding blacklist BED files.
    """
    unique_vcfs = {}
    
    # Collect unique VCF files
    for row in rows:
        sample_name, vcf_file, _ = row
        if vcf_file not in unique_vcfs:
            if not Path(vcf_file).exists():
                raise FileNotFoundError(f"Missing VCF file: {vcf_file}")
            unique_vcfs[vcf_file] = None
    
    # Generate blacklist for each unique VCF using VCF prefix as the output name
    print(f'Generating {len(unique_vcfs)} unique blacklist(s)...')
    for vcf_file in unique_vcfs:
        vcf_prefix = get_vcf_prefix(vcf_file)
        blacklist_bedfile = os.path.join(output_dir, f"{vcf_prefix}.CpG_variants_blacklist.bed")
        
        if not Path(blacklist_bedfile).exists():
            print(f'  Generating blacklist for VCF: {vcf_file}', flush=True)
            vcf_bed = os.path.join(output_dir, f"{vcf_prefix}.vcf.tmp.bed")
            vcf_to_bed(vcf_file, vcf_bed)
            
            vcf_bt = pybedtools.BedTool(vcf_bed)
            bed_bt = pybedtools.BedTool(cpg_bedfile)
            bed_bt.intersect(vcf_bt, u=True).saveas(blacklist_bedfile)
            
            os.remove(vcf_bed)
            print(f'  Blacklist written to {blacklist_bedfile}', flush=True)
        else:
            print(f'  Using existing blacklist for VCF: {vcf_file}', flush=True)
        
        unique_vcfs[vcf_file] = blacklist_bedfile
    
    return unique_vcfs

def process_sample_with_blacklist(output_dir, has_header, row):
    """Wrapper to pass the pre-generated blacklist to process_sample"""
    sample_name, vcf_file, allc_file, blacklist_bedfile = row
    process_sample(output_dir, has_header, blacklist_bedfile, (sample_name, vcf_file, allc_file))

def run_parallel(samplesheet, cpg_bedfile, output_dir, ncpus, has_header, blacklist_provided):

    # parse the files that only needs to be read once
    os.makedirs(output_dir, exist_ok=True)
    rows = parse_samplesheet(samplesheet)

    if blacklist_provided:
        # When blacklist_provided=True, the 2nd column is already a blacklist BED file
        print('Processing samples with provided blacklists...')
        
        # Create rows with the blacklist file path (2nd column is the blacklist)
        rows_with_blacklist = [(row[0], row[1], row[2], row[1]) for row in rows]
        process_func = partial(process_sample_with_blacklist, output_dir, has_header)
        
        with Pool(ncpus) as pool:
            results = pool.map(process_func, rows_with_blacklist)
    else:
        if cpg_bedfile is None:
            raise FileNotFoundError("Missing cpg_bedfile - required for blacklist generation")
        
        # Generate blacklists once per unique VCF file
        vcf_to_blacklist = generate_vcf_blacklists(rows, cpg_bedfile, output_dir)
        
        # Now process all samples in parallel, using pre-generated blacklists
        print('Start processing samples...')
        
        # Create modified rows with blacklist paths
        rows_with_blacklist = []
        for row in rows:
            sample_name, vcf_file, allc_file = row
            blacklist_bedfile = vcf_to_blacklist[vcf_file]
            rows_with_blacklist.append((sample_name, vcf_file, allc_file, blacklist_bedfile))
        
        process_func = partial(process_sample_with_blacklist, output_dir, has_header)
        
        with Pool(ncpus) as pool:
            results = pool.map(process_func, rows_with_blacklist)

def process_bam_sample(output_dir, blacklist_bedfile, row):
    """Process a single BAM sample using a pre-generated or provided blacklist"""
    sample_name, vcf_file, bam_file = row
    output_bam_file = f"{sample_name}.var_masked.bam"
    
    if not Path(bam_file).exists():
        raise FileNotFoundError(f"Missing BAM file: {bam_file}")

    # remove those cytosines from the BAM file
    print(f'Start processing {bam_file}', flush=True)
    mask_bam(bam_file, blacklist_bedfile, output_dir, output_bam_file)
    print(f'Result written to {output_bam_file}', flush=True)

def process_bam_sample_with_blacklist(output_dir, row):
    """Wrapper to pass the pre-generated blacklist to process_bam_sample"""
    sample_name, vcf_file, bam_file, blacklist_bedfile = row
    process_bam_sample(output_dir, blacklist_bedfile, (sample_name, vcf_file, bam_file))

def run_bam_parallel(samplesheet, cpg_bedfile, output_dir, ncpus, blacklist_provided):

    os.makedirs(output_dir, exist_ok=True)
    rows = parse_samplesheet(samplesheet)

    if blacklist_provided:
        # When blacklist_provided=True, the 2nd column is already a blacklist BED file
        print('Processing BAM files with provided blacklists...')
        rows_with_blacklist = [(row[0], row[1], row[2], row[1]) for row in rows]
        process_func = partial(process_bam_sample_with_blacklist, output_dir)
        
        with Pool(ncpus) as pool:
            results = pool.map(process_func, rows_with_blacklist)
    else:
        if cpg_bedfile is None:
            raise FileNotFoundError("Missing cpg_bedfile - required for blacklist generation")
        
        # Generate blacklists once per unique VCF file
        vcf_to_blacklist = generate_vcf_blacklists(rows, cpg_bedfile, output_dir)
        
        # Now process all BAM files in parallel, using pre-generated blacklists
        print('Start processing BAM files...')
        
        # Create modified rows with blacklist paths
        rows_with_blacklist = []
        for row in rows:
            sample_name, vcf_file, bam_file = row
            blacklist_bedfile = vcf_to_blacklist[vcf_file]
            rows_with_blacklist.append((sample_name, vcf_file, bam_file, blacklist_bedfile))
        
        process_func = partial(process_bam_sample_with_blacklist, output_dir)
        
        with Pool(ncpus) as pool:
            results = pool.map(process_func, rows_with_blacklist)

def process_blacklist_sample(cpg_bed, output_dir, row):
    sample_name, vcf_file = row
    
    if not Path(vcf_file).exists():
        raise FileNotFoundError(f"Missing VCF file: {vcf_file}")
    
    if not Path(cpg_bed).exists():
        raise FileNotFoundError(f"Missing CpG BED file: {cpg_bed}")
    
    print(f'Generating blacklist for {sample_name}', flush=True)
    blacklist_bedfile = find_overlapping_positions(vcf_file, cpg_bed, output_dir, sample_name)
    print(f'Blacklist BED file written to {blacklist_bedfile}', flush=True)
    return blacklist_bedfile

def run_blacklist_parallel(samplesheet, cpg_bedfile, output_dir, ncpus):
    """Generate blacklist BED files for each sample"""
    
    os.makedirs(output_dir, exist_ok=True)
    rows = parse_samplesheet(samplesheet)
    
    if cpg_bedfile is None:
        raise FileNotFoundError("Missing cpg_bedfile - required for blacklist generation")
    
    print('Start generating blacklists...')
    process_func = partial(process_blacklist_sample, cpg_bedfile, output_dir)
    
    with Pool(ncpus) as pool:
        results = pool.map(process_func, rows)
    
    print(f'Generated {len(results)} blacklist BED files in {output_dir}')

