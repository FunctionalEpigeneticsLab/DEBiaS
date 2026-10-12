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
    #        pos1 = int(cols[1])+1
    #        pos2 = int(cols[2])
    #        position_groups.append({(chrom, pos1), (chrom, pos2)})
    #positions = set().union(*position_groups)
            start = int(cols[1])
            end = int(cols[2])
            for pos1 in range(start + 1, end + 1):
                positions.add((chrom, pos1))

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
            #out.write(f"{chrom}\t{pos-1}\t{pos}\n")
            ref = cols[3]
            start = pos - 1
            out.write(f"{chrom}\t{start}\t{start + len(ref)}\n")

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

            # mark all reference positions spanned by this variant as affected
            ref_span = max(1, len(ref))
            for p in range(pos1, pos1 + ref_span):
                affected_positions.add((chrom, p))

            for alt in alts:
                alt = alt.upper()
                if alt == '.':
                    continue

                # simply record every ALT base at its corresponding genomic position.
                # SNV:     len(ref)=1, len(alt)=1  → records alt[0] at pos1
                # MNV:     len(ref)>1, len(alt)>1  → records each alt base at pos1+i
                # insertion with anchor: alt[0]==ref[0], len(alt)>len(ref) -> records anchor (unchanged) + inserted bases
                # deletion: len(alt) < len(ref)  → records fewer bases than ref span (affected_positions already marks deleted positions)
                for i, base in enumerate(alt):
                    if base in {'A', 'C', 'G', 'T'}:
                        snv_alt_by_pos[(chrom, pos1 + i)].add(base)

    return snv_alt_by_pos, affected_positions


def build_ch_blacklist_from_bed(vcf_file, ch_bedfile, out_bed_file, output_dir):
    """
    build a CH blacklist BED from VCF + pre-computed CH BED
    rules implemented:
        1) direct variant/indel at CH cytosine position (BED start/end)
        2) context-changing SNP at context base: strand '+' context ALT includes G (CH -> CG), OR strand '-' context ALT includes C (CH -> CG on reverse strand)
        3) indel/complex variant at immediate context base (deprecated - not masked now)
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

            if direct_variant or context_to_cg:
                masked.add((chrom, start, end))

    with open(out_bed_file, 'w') as out:
        for chrom, start, end in sorted(masked, key=lambda x: (x[0], x[1])):
            out.write(f"{chrom}\t{start}\t{end}\n")

    return out_bed_file


def allc_to_bed(allc_file, out_bed_file, has_header):

    input_is_gz = allc_file.endswith('.gz')
    open_in = gzip.open if input_is_gz else open

    skipped = 0
    with open_in(allc_file, 'rt') as f_in, open(out_bed_file, 'w') as f_out:
        for i, line in enumerate(f_in):

            if has_header and i == 0:
                continue

            stripped = line.strip()
            if not stripped:
                skipped += 1
                continue

            fields = stripped.split('\t')

            if len(fields) < 2 or fields[0] == '' or fields[1] == '':
                print(f"WARNING: skipping malformed line {i+1} in {allc_file}: {repr(line.rstrip())}",
                      flush=True)
                skipped += 1
                continue

            chrom = fields[0]

            try:
                pos = int(fields[1])
            except ValueError:
                print(f"WARNING: skipping line {i+1} in {allc_file}: "
                      f"non-integer position {repr(fields[1])}: {repr(line.rstrip())}",
                      flush=True)
                skipped += 1
                continue

            # Skip invalid 1-based coordinates
            if pos < 1:
                print(f"WARNING: skipping line {i+1} in {allc_file}: "
                      f"pos={pos} is not a valid 1-based coordinate",
                      flush=True)
                skipped += 1
                continue

            rest = "\t".join(fields)
            f_out.write(f"{chrom}\t{pos-1}\t{pos}\t{rest}\n")

    if skipped > 0:
        print(f"WARNING: {allc_file} — skipped {skipped} malformed line(s) total", flush=True)

def _report_malformed_bed(bed_file, error_msg):
    print(f"\n{'='*60}", flush=True)
    print(f"ERROR: MalformedBedLineError in: {bed_file}", flush=True)
    print(f"pybedtools says: {error_msg}", flush=True)
    print(f"Scanning file for malformed lines...", flush=True)
    print(f"{'='*60}", flush=True)

    issues_found = 0
    try:
        with open(bed_file, 'r') as f:
            for i, line in enumerate(f, start=1):
                stripped = line.rstrip('\n')
                fields = stripped.split('\t')

                reason = None

                if stripped == '' or stripped.isspace():
                    reason = "blank/whitespace line"
                elif len(fields) < 3:
                    reason = f"too few fields ({len(fields)}, need >=3)"
                else:
                    chrom, start, end = fields[0], fields[1], fields[2]
                    if chrom == '':
                        reason = "empty chrom"
                    elif start == '' or end == '':
                        reason = f"empty start or end field: start={repr(start)}, end={repr(end)}"
                    else:
                        try:
                            s, e = int(start), int(end)
                            if s < 0:
                                reason = f"negative start: {s}"
                            elif e < 0:
                                reason = f"negative end: {e}"
                            elif s > e:
                                reason = f"start > stop: {s} > {e}"
                        except ValueError:
                            reason = f"non-integer start/end: start={repr(start)}, end={repr(end)}"

                if reason:
                    print(f"  Line {i:>6d}: [{reason}]  →  {repr(stripped)}", flush=True)
                    issues_found += 1

                    if issues_found >= 20:
                        print(f"  ... (stopping after 20 issues)", flush=True)
                        break

    except Exception as scan_err:
        print(f"  Could not scan file: {scan_err}", flush=True)

    if issues_found == 0:
        print(f"  No obvious issue found", flush=True)

    print(f"{'='*60}\n", flush=True)

def mask_allc_with_bedtools(allc_file, has_header, blacklist_bedfile, output_dir, output_file):
    os.makedirs(output_dir, exist_ok=True)
    allc_bed = os.path.join(output_dir, f"{allc_file}.tmp.bed")
    allc_to_bed(allc_file, allc_bed, has_header)

    allc_bedtool = pybedtools.BedTool(allc_bed)
    mask_result = allc_bedtool.intersect(blacklist_bedfile, v=True)

    # count total — catch malformed lines with context
    try:
        num_total = sum(1 for _ in allc_bedtool)
    except pybedtools.cbedtools.MalformedBedLineError as e:
        _report_malformed_bed(allc_bed, str(e))
        raise

    # count retained — catch malformed lines with context
    try:
        num_retained = sum(1 for _ in mask_result)
    except pybedtools.cbedtools.MalformedBedLineError as e:
        _report_malformed_bed(allc_bed, str(e))
        raise

    num_removed = num_total - num_retained
    percent_removed = round((num_removed / num_total) * 100, 2)

    # write result
    try:
        with gzip.open(os.path.join(output_dir, output_file), 'wt') as out:
            for feature in mask_result:
                out.write("\t".join(feature[3:]) + "\n")
    except pybedtools.cbedtools.MalformedBedLineError as e:
        _report_malformed_bed(allc_bed, str(e))
        raise

    if os.environ.get('DEBIAS_KEEP_TMP'):
        print(f"DEBUG: keeping tmp file: {allc_bed}", flush=True)
    else:
        os.remove(allc_bed)

    report_file = os.path.join(output_dir, output_file.replace('.gz', '').replace('.tsv', '.report.tsv'))
    with open(report_file, "w") as report:
        report.write(f"Total input lines\t{num_total}\n")
        report.write(f"Lines retained\t{num_retained}\n")
        report.write(f"Lines removed\t{num_removed}\n")
        report.write(f"Percentage removed\t{percent_removed}%\n")

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

