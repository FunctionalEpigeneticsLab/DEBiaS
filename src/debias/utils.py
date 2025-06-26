#!/usr/bin/env python3

import gzip
import csv
import sys

def parse_samplesheet(input):
    with open(input, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)
        rows = [row for row in reader]
    return rows

def parse_bed(bed_file):
    position_groups = []
    with open(bed_file, 'r') as f:
        for line in f:
            cols = line.split()
            chrom = str(cols[0])
            start = int(cols[1])+1
            end = int(cols[2])
            position_groups.append({(chrom, pos) for pos in range(start, end + 1)})
    return position_groups

def parse_vcf(vcf_file):
    positions = set()
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            chrom = str(cols[0])
            pos = int(cols[1])
            positions.add((chrom, pos))
    return positions

def find_overlapping_positions(vcf_positions, bed_groups):
    
    output_positions = set()
    for group in bed_groups:
        if any(pos in vcf_positions for pos in group):
            output_positions.update(group)
            #print(f'Blacklisted locus: {group}', file=sys.stdout, flush=True)
    sorted_positions = sorted(output_positions, key=lambda x: (x[0], x[1]))
    return sorted_positions

def mask_allc(allc_file, mc_file_type, blacklisted_loci, output_file):
    print(f'# blacklisted loci: {len(blacklisted_loci)}', file=sys.stdout, flush=True)

    has_header=False

    if mc_file_type=='cov':
        has_header=True

    open_func = gzip.open if allc_file.endswith('.gz') else open

    with open_func(allc_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        first = True
        for line in infile:
            if first and has_header:
                outfile.write(line)
                first = False
                continue
            first = False
            fields = line.strip().split('\t')
            chrom = str(fields[0])
            pos = int(fields[1])

            if (chrom, pos) not in blacklisted_loci:
                outfile.write(line)
            else:
                print(f'Excluding {chrom}_{pos}', file=sys.stdout, flush=True)

