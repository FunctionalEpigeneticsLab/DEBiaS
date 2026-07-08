import gzip
import pybedtools
import os
import csv
import sys
import io
import pysam
from pathlib import Path
from collections import defaultdict

def get_vcf_prefix(vcf_file):
    basename = os.path.basename(vcf_file)
    # Remove .vcf or .vcf.gz extension
    if basename.endswith('.vcf.gz'):
        return basename[:-7]
    elif basename.endswith('.vcf'):
        return basename[:-4]
    return basename

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
            pos1 = int(cols[1])+1 # because BED is 0-based
            pos2 = int(cols[2])
            position_groups.append({(chrom, pos1), (chrom, pos2)})
    positions = set().union(*position_groups)
    return positions

def vcf_to_bed(vcf_file, out_bed_file):
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    with open_func(vcf_file, 'rt') as f, open(out_bed_file, 'w') as out:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            chrom = cols[0]
            pos = int(cols[1]) # 1-based position
            out.write(f"{chrom}\t{pos-1}\t{pos}\n") # writing to 0-based BED

def find_overlapping_positions(vcf_file, cpg_bed, output_dir, sample_name):
    vcf_bed=os.path.join(output_dir, f"{sample_name}.vcf.tmp.bed")
    vcf_to_bed(vcf_file, vcf_bed)
    vcf_bt = pybedtools.BedTool(vcf_bed)
    bed_bt = pybedtools.BedTool(cpg_bed) # cpg_bed must be 0-based also
    blacklist_bedfile = os.path.join(output_dir, f"{sample_name}.CpG_variants_blacklist.bed")
    bed_bt.intersect(vcf_bt, u=True).saveas(blacklist_bedfile)
    os.remove(vcf_bed)
    return blacklist_bedfile


def parse_vcf_for_ch_masking(vcf_file):
    snv_alt_by_pos = defaultdict(set)
    affected_positions = set()
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                continue

            chrom = cols[0]
            pos1 = int(cols[1])
            ref = cols[3].upper()
            alts = cols[4].split(',')
            ref_span = max(1, len(ref))
            for p in range(pos1, pos1 + ref_span):
                affected_positions.add((chrom, p))

            for alt in alts:
                alt = alt.upper()
                if alt == '.':
                    continue
                if len(ref) == 1 and len(alt) == 1 and alt in {'A', 'C', 'G', 'T'}:
                    snv_alt_by_pos[(chrom, pos1)].add(alt)

    return snv_alt_by_pos, affected_positions


def build_ch_blacklist_from_bed(vcf_file, ch_bedfile, out_bed_file, output_dir):
    """
    build a CH blacklist BED from VCF + pre-computed CH BED
    rules implemented:
        1) direct variant/indel at CH cytosine position (BED start/end)
        2) context-changing SNP at context base:
           strand '+' context ALT includes G (CH -> CG)
           strand '-' context ALT includes C (CH -> CG on reverse strand)
        3) indel/complex variant at immediate context base
    """
    if not Path(ch_bedfile).exists():
        raise FileNotFoundError(f"Missing CH BED file: {ch_bedfile}")

    snv_alt_by_pos, affected_positions = parse_vcf_for_ch_masking(vcf_file)
    masked = set()

    with open(ch_bedfile, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue

            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[3] if len(cols) >= 4 else '+'

            cyt_pos1 = start + 1

            direct_variant = (chrom, cyt_pos1) in affected_positions

            if strand == '-':
                # reverse CH: context is one base to the left
                ctx_pos1 = start
                context_to_cg = 'C' in snv_alt_by_pos.get((chrom, ctx_pos1), set())
            else:
                # forward CH: context is one base to the right
                ctx_pos1 = end + 1
                context_to_cg = 'G' in snv_alt_by_pos.get((chrom, ctx_pos1), set())

            context_indel_or_complex = (chrom, ctx_pos1) in affected_positions and not snv_alt_by_pos.get((chrom, ctx_pos1), set())

            if direct_variant or context_to_cg or context_indel_or_complex:
                masked.add((chrom, start, end))

    with open(out_bed_file, 'w') as out:
        for chrom, start, end in sorted(masked, key=lambda x: (x[0], x[1])):
            out.write(f"{chrom}\t{start}\t{end}\n")

    return out_bed_file


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
            f_out.write(f"{chrom}\t{pos-1}\t{pos}\t{rest}\n")


def mask_allc_with_bedtools(allc_file, has_header, blacklist_bedfile, output_dir, output_file):
    os.makedirs(output_dir, exist_ok=True)
    allc_bed = os.path.join(output_dir, f"{allc_file}.tmp.bed")
    allc_to_bed(allc_file, allc_bed, has_header)

    # subtract
    allc_bedtool = pybedtools.BedTool(allc_bed)
    mask_result = allc_bedtool.intersect(blacklist_bedfile, v=True)

    # count before and after
    num_total = sum(1 for _ in allc_bedtool)
    num_retained = sum(1 for _ in mask_result)
    num_removed = num_total-num_retained
    percent_removed = round((num_removed / num_total) * 100, 2)

    # write result
    with gzip.open(os.path.join(output_dir, output_file), 'wt') as out:
        for feature in mask_result:
            # feature[3:] are the original columns from allc
            out.write("\t".join(feature[3:]) + "\n")

    # remove temp file
    os.remove(allc_bed)

    # write report
    report_file=os.path.join(output_dir, output_file.replace('.gz','').replace('.tsv','.report.tsv'))
    with open(report_file, "w") as report:
        report.write(f"Total input lines\t{num_total}\n")
        report.write(f"Lines retained\t{num_retained}\n")
        report.write(f"Lines removed\t{num_removed}\n")
        report.write(f"Percentage removed: {percent_removed}%\n")

def mask_bam(bam_file, blacklist_bedfile, output_dir, output_bam_file):
    positions = parse_bed(blacklist_bedfile)

    index_file = bam_file + ".bai"
    if not os.path.exists(index_file):
        pysam.index(bam_file)

    output_file = os.path.join(output_dir, output_bam_file)
    count_all=0
    count_modified=0

    with pysam.AlignmentFile(bam_file, "rb") as in_bam, pysam.AlignmentFile(output_file, "wb", header=in_bam.header) as out_bam:

        for read in in_bam:
            modified = False
            if not read.is_unmapped:
                count_all+=1
                seq = list(read.query_sequence)
                ref_positions = read.get_reference_positions(full_length=True)

                for query_index, ref_pos in enumerate(ref_positions):
                    if ref_pos is not None and (read.reference_name, ref_pos + 1) in positions:
                        seq[query_index] = 'N'
                        modified = True

                read.query_sequence = ''.join(seq)

            out_bam.write(read)
            if modified:
                count_modified+=1
    
    pysam.index(output_file)
    percent_modified = round((count_modified/count_all) * 100, 2)
    
    # write report
    report_file=os.path.join(output_dir, output_bam_file.replace('.bam','.report.tsv'))
    with open(report_file, "w") as report:
        report.write(f"Total input reads\t{count_all}\n")
        report.write(f"Reads unchanged\t{int(count_all-count_modified)}\n")
        report.write(f"Reads modified\t{count_modified}\n")
        report.write(f"Percentage modified: {percent_modified}%\n")

