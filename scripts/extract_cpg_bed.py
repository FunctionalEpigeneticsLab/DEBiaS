#!/usr/bin/env python3
"""
Extract 0-based CpG locations from a genome FASTA file and output in BED format.

Usage:
    python extract_cpg_bed.py -i genome.fa -o cpg_sites.bed

Requirements:
    - BioPython (pip install biopython)
"""

import argparse
import sys
import re
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def find_cpg_sites(sequence, chromosome):
    """
    Find all CpG sites in a sequence using efficient regex.
    """
    seq_upper = sequence.upper()
    cpg_sites = []
    for match in re.finditer(r'(?=CG)', seq_upper):
        start_pos = match.start()  # 0-based start position
        end_pos = start_pos + 2    # end position (start + 2 for CG dinucleotide)
        cpg_sites.append((chromosome, start_pos, end_pos))
    
    return cpg_sites

def process_chromosome(record, temp_dir, progress_lock):
    """
    Process a single chromosome and write CpG sites to a temporary file.
    """
    chromosome = record.id
    sequence = str(record.seq)
    
    temp_file = os.path.join(temp_dir, f"{chromosome}.bed")
    cpg_sites = find_cpg_sites(sequence, chromosome)
    with open(temp_file, 'w') as f:
        for chr_name, start, end in cpg_sites:
            f.write(f"{chr_name}\t{start}\t{end}\n")
    with progress_lock:
        print(f"Completed chromosome: {chromosome} (length: {len(sequence):,} bp, CpGs: {len(cpg_sites):,})", file=sys.stdout)
    return temp_file, chromosome, len(cpg_sites)

def process_fasta_file(input_file, output_file, max_workers=None):
    """
    Process FASTA file and extract CpG sites to BED format using parallel processing.
    """
    if max_workers is None:
        max_workers = min(32, os.cpu_count() + 4)
    
    print(f"Processing FASTA file: {input_file}", file=sys.stdout)
    print(f"Output BED file: {output_file}", file=sys.stdout)
    print(f"Using {max_workers} worker threads", file=sys.stdout)
    
    total_cpg_count = 0
    progress_lock = Lock()
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            print(f"Created temporary directory: {temp_dir}", file=sys.stdout)
            
            print("Loading FASTA records...", file=sys.stdout)
            records = list(SeqIO.parse(input_file, "fasta"))
            print(f"Found {len(records)} chromosomes/contigs", file=sys.stdout)
            
            temp_files = []
            chromosome_stats = []
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_record = {
                    executor.submit(process_chromosome, record, temp_dir, progress_lock): record 
                    for record in records
                }
                for future in as_completed(future_to_record):
                    try:
                        temp_file, chromosome, cpg_count = future.result()
                        temp_files.append(temp_file)
                        chromosome_stats.append((chromosome, cpg_count))
                        total_cpg_count += cpg_count
                    except Exception as e:
                        record = future_to_record[future]
                        print(f"Error processing chromosome {record.id}: {e}", file=sys.stderr)
                        raise
            
            temp_files.sort(key=lambda x: os.path.basename(x))
            
            print("Combining results into final BED file...", file=sys.stdout)
            
            with open(output_file, 'w') as bed_out:
                for temp_file in temp_files:
                    if os.path.exists(temp_file):
                        with open(temp_file, 'r') as temp_in:
                            bed_out.write(temp_in.read())
    
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print("\n" + "="*50, file=sys.stdout)
    print("PROCESSING SUMMARY", file=sys.stdout)
    print("="*50, file=sys.stdout)
    for chromosome, cpg_count in sorted(chromosome_stats):
        print(f"{chromosome:>15}: {cpg_count:>12,} CpG sites", file=sys.stdout)
    print("-"*50, file=sys.stdout)
    print(f"{'Total':>15}: {total_cpg_count:>12,} CpG sites", file=sys.stdout)
    print(f"BED file written successfully: {output_file}", file=sys.stdout)

def validate_fasta_file(filename):
    """
    Validate that the input file is a valid FASTA file.
    """
    try:
        with open(filename, 'r') as handle:
            first_record = next(SeqIO.parse(handle, "fasta"))
            return True
    except:
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Extract 0-based CpG locations from genome FASTA file to BED format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python extract_cpg_bed.py -i genome.fa -o cpg_sites.bed --threads 16

Performance:
    The script uses parallel processing to handle multiple chromosomes simultaneously.
    Default thread count is auto-detected based on your CPU.

Output format:
    BED format with 0-based coordinates
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input FASTA file'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output BED file for CpG locations'
    )
    
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=None,
        help='Number of threads to use (default: auto-detect)'
    )
    
    args = parser.parse_args()
    
    if not validate_fasta_file(args.input):
        print(f"Error: '{args.input}' is not a valid FASTA file.", file=sys.stderr)
        sys.exit(1)
    
    process_fasta_file(args.input, args.output, args.threads)


if __name__ == "__main__":
    main()

