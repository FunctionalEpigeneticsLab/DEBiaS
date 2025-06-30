import os
import sys
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from debias.utils import *

def process_sample(cpg_chrom_pos_list, output_dir, mc_file_type, row):
    sample_name,vcf_file,allc_file = row

    if not Path(vcf_file).exists() or not Path(allc_file).exists():
        raise FileNotFoundError(f"Missing: {' '.join(str(f) for f in [vcf, allc] if not f.exists())}")

    # step 1: subset to CpG sites variants and include adjacent cytosine
    # TODO: add option to provide input directly
    variant_loci = parse_vcf(vcf_file)
    print(f'Finished parsing {vcf_file}', flush=True)
    loci_blacklist = find_overlapping_positions(variant_loci, cpg_chrom_pos_list)
    del(variant_loci)
    blacklist_bedfile = os.path.join(output_dir, f"{sample_name}.CpG_variants_blacklist.bed")
    with open(blacklist_bedfile,'w') as f:
        for chrom,pos in loci_blacklist:
            f.write(f"{chrom}\t{pos}\t{pos+1}\n")

    # step 2: remove those cytosines from the allc file
    print(f'Start processing {allc_file}', flush=True)
    output_allc_file = os.path.join(output_dir, f"{sample_name}.var_masked.allc.tsv")
    mask_allc_with_bedtools(allc_file, mc_file_type, blacklist_bedfile, output_dir, output_allc_file)
    print(f'Result written to {output_allc_file}', flush=True)


def run_parallel(samplesheet, cpg_bedfile, output_dir, ncpus, mc_file_type):

    # parse the files that only needs to be read once
    os.makedirs(output_dir, exist_ok=True)
    rows = parse_samplesheet(samplesheet)
    print('Parsing input BED file...')
    cpg_chrom_pos_list = parse_bed(cpg_bedfile)

    # run in parallel
    print('Start processing samples...')
    process_func = partial(process_sample, cpg_chrom_pos_list, output_dir, mc_file_type)

    with Pool(ncpus) as pool:
        #results = list(tqdm(pool.imap(process_func, rows), total=len(rows), desc="Processing"))
        results = pool.map(process_func, rows)

