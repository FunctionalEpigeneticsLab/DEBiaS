#!/usr/bin/env python3
"""
Extract 0-based CH cytosine locations from a genome FASTA into one BED file.

Output columns:
    chrom  start  end  strand

Where `start:end` is the 1-bp cytosine-call position to mask in CH mode:
    strand '+' : forward CH cytosine (reference base C, context at +1)
    strand '-' : reverse CH cytosine represented by reference base G (context at -1)

This single BED is reused across runs with DEBiaS `--CH --ch_bedfile`.

Usage:
    python extract_ch_bed.py -i genome.fa -o ch_sites.bed
"""

import argparse
import sys
from Bio import SeqIO


def write_ch_sites_for_record(record, out_handle):
    """Write CH sites for one FASTA record and return count."""
    chrom = record.id
    seq = str(record.seq).upper()
    n = len(seq)
    count = 0

    for i, base in enumerate(seq):
        if base == 'C' and i + 1 < n and seq[i + 1] != 'G':
            out_handle.write(f"{chrom}\t{i}\t{i + 1}\t+\n")
            count += 1
        elif base == 'G' and i - 1 >= 0 and seq[i - 1] != 'C':
            out_handle.write(f"{chrom}\t{i}\t{i + 1}\t-\n")
            count += 1

    return count


def process_fasta_file(input_file, output_file):
    """
    Stream FASTA records sequentially to keep memory usage low.
    """
    print(f"Processing FASTA file: {input_file}", file=sys.stdout)
    print(f"Output BED file: {output_file}", file=sys.stdout)

    total = 0
    per_chrom = []

    try:
        with open(output_file, 'w') as out:
            for record in SeqIO.parse(input_file, "fasta"):
                ch_count = write_ch_sites_for_record(record, out)
                total += ch_count
                per_chrom.append((record.id, ch_count, len(record.seq)))
                print(
                    f"Completed chromosome: {record.id} "
                    f"(length: {len(record.seq):,} bp, CH sites: {ch_count:,})",
                    file=sys.stdout,
                )
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n" + "=" * 50, file=sys.stdout)
    print("PROCESSING SUMMARY", file=sys.stdout)
    print("=" * 50, file=sys.stdout)
    for chrom, count, _ in per_chrom:
        print(f"{chrom:>15}: {count:>12,} CH sites", file=sys.stdout)
    print("-" * 50, file=sys.stdout)
    print(f"{'Total':>15}: {total:>12,} CH sites", file=sys.stdout)
    print(f"BED file written successfully: {output_file}", file=sys.stdout)


def validate_fasta_file(filename):
    """Validate that the input file is a valid FASTA file."""
    try:
        with open(filename, 'r') as handle:
            next(SeqIO.parse(handle, "fasta"))
            return True
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Extract 0-based CH cytosine locations from genome FASTA to BED",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python extract_ch_bed.py -i genome.fa -o CHsites.bed

Output format:
    4-column BED: chrom, start, end, strand
        """,
    )

    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file for CH locations')
    args = parser.parse_args()

    if not validate_fasta_file(args.input):
        print(f"Error: '{args.input}' is not a valid FASTA file.", file=sys.stderr)
        sys.exit(1)

    process_fasta_file(args.input, args.output)


if __name__ == "__main__":
    main()
