import gzip
import pybedtools
import tempfile
import os
import csv
import sys
import io
import pysam

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

def vcf_to_bed(vcf_file, out_bed_file):
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    with open_func(vcf_file, 'rt') as f, open(out_bed_file, 'w') as out:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            out.write(f"{chrom}\t{pos}\t{pos+1}\n")


def allc_to_bed(allc_file, out_bed_file, has_header):
    input_is_gz = allc_file.endswith('.gz')
    #has_header = (mc_file_type == 'cov')

    open_in = gzip.open if input_is_gz else open

    with open_in(allc_file, 'rt') as f_in, open(out_bed_file, 'w') as f_out:
        for i, line in enumerate(f_in):
            if has_header and i == 0:
                continue  # skip header for BED
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            rest = "\t".join(fields)
            f_out.write(f"{chrom}\t{pos}\t{pos+1}\t{rest}\n")


def mask_allc_with_bedtools(allc_file, has_header, blacklist_bedfile, output_dir, output_file):
    #with tempfile.TemporaryDirectory() as tmpdir:
    allc_bed = os.path.join(output_dir, f"{allc_file}.tmp.bed")
    allc_to_bed(allc_file, allc_bed, has_header)

    # subtract
    allc_bedtool = pybedtools.BedTool(allc_bed)
    mask_result = allc_bedtool.intersect(blacklist_bedfile, v=True)

    with open(os.path.join(output_dir, output_file), 'w') as out:
        for feature in mask_result:
            # feature[3:] are the original columns from allc
            out.write("\t".join(feature[3:]) + "\n")

    os.remove(allc_bed)



def mask_bam(bam_file, vcf_file, output_dir, output_bam_file):
    positions = parse_vcf(vcf_file)
    output_file = os.path.join(output_dir, output_bam_file)
    #modified_output_file = f"{output_dir}/{output_bam_file}.modified_reads.bam"
    with pysam.AlignmentFile(bam_file, "rb") as in_bam, \
            pysam.AlignmentFile(output_file, "wb", header=in_bam.header) as out_bam:
    #    pysam.AlignmentFile(modified_output_file, "wb", header=in_bam.header) as modified_bam:

        for read in in_bam:
            modified = False
            if not read.is_unmapped:
                seq = list(read.query_sequence)
                ref_positions = read.get_reference_positions(full_length=True)

                for query_index, ref_pos in enumerate(ref_positions):
                    if ref_pos is not None and (read.reference_name, ref_pos + 1) in positions:
                        seq[query_index] = 'N'
                        modified = True

                read.query_sequence = ''.join(seq)

            out_bam.write(read)
            #if modified:
            #    modified_bam.write(read)
    # TODO: write to a summary file: how many reads (percentage) is modified due to variants
