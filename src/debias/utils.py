import gzip
import pybedtools
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
            pos = int(cols[1]) # VCF position is 1-based exact position
            out.write(f"{chrom}\t{pos-1}\t{pos}\n") # writing to 0-based BED

def find_overlapping_positions(vcf_file, cpg_bed, output_dir, sample_name):
    vcf_bed=os.path.join(output_dir, f"{sample_name}.vcf.tmp.bed")
    vcf_to_bed(vcf_file, vcf_bed)

    vcf_bt = pybedtools.BedTool(vcf_bed)
    bed_bt = pybedtools.BedTool(cpg_bed) # cpg_bed must be 0-based also

    blacklist_bedfile = os.path.join(output_dir, f"{sample_name}.CpG_variants_blacklist.bed")
    bed_bt.intersect(vcf_bt, u=True).saveas(blacklist_bedfile)

    # remove temp file
    os.remove(vcf_bed)

    return blacklist_bedfile

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
    #modified_output_file = f"{output_dir}/{output_bam_file}.modified_reads.bam"

    count_all=0
    count_modified=0

    with pysam.AlignmentFile(bam_file, "rb") as in_bam, pysam.AlignmentFile(output_file, "wb", header=in_bam.header) as out_bam:
    # pysam.AlignmentFile(modified_output_file, "wb", header=in_bam.header) as modified_bam:

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
            #    modified_bam.write(read)
    
    pysam.index(output_file)
    percent_modified = round((count_modified/count_all) * 100, 2)
    
    # write report
    report_file=os.path.join(output_dir, output_bam_file.replace('.bam','.report.tsv'))
    with open(report_file, "w") as report:
        report.write(f"Total input reads\t{count_all}\n")
        report.write(f"Reads unchanged\t{int(count_all-count_modified)}\n")
        report.write(f"Reads modified\t{count_modified}\n")
        report.write(f"Percentage modified: {percent_modified}%\n")

