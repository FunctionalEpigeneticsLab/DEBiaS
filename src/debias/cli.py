import argparse
from debias.core import run_parallel, run_bam_parallel, run_blacklist_parallel

# multi-command CLI
def main():
    # param common for different commands
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('--output_dir', required=True, help='directory to write output')
    common_parser.add_argument('--ncpus', type=int, default=8, help='number of threads to use')
    common_parser.add_argument('--cpg_bedfile', default=None, required=False, help="input CpG BED file")
    common_parser.add_argument('--blacklist_provided', action='store_true', help="the 2nd column provided in 'input' is a blacklist that does not need further matching to a 'cpg_bedfile'. Enabling this option will skip generating a blacklist, so 'cpg_bedfile' is not needed in this senario.")

    # main parser
    parser = argparse.ArgumentParser(description="A tool to filter variant loci for (pseudo)bulk bisulfite sequencing data")
    subparsers = parser.add_subparsers(dest='command', required=True)

    # filter_mc
    parser_filter_mc = subparsers.add_parser('filter_mc', parents=[common_parser], help='Filter methylation calling data')
    parser_filter_mc.add_argument('--input', required=True, help='csv file with columns: [sample,vcf_file,allc_file]. "sample" will be used as prefix for the output.')
    parser_filter_mc.add_argument('--has_header', action='store_true', help='first line of methylation calling file is header.')
    # TODO: add an option in the workflow to have a input of blacklist without retriving from cpg bed file

    # filter_bam
    parser_filter_bam = subparsers.add_parser('filter_bam', parents=[common_parser], help='Filter BAM file')
    parser_filter_bam.add_argument('--input', required=True, help='csv file with columns: [sample,vcf_file,bam_file]. "sample" will be used as prefix for the output.')

    # generate_blacklist
    parser_generate_blacklist = subparsers.add_parser('generate_blacklist', parents=[common_parser], help='Generate blacklist BED files from VCF files')
    parser_generate_blacklist.add_argument('--input', required=True, help='csv file with columns: [sample,vcf_file]. "sample" will be used as prefix for the output blacklist BED files.')

    args = parser.parse_args()

    # dispatch logic
    if args.command == 'filter_mc':
        print("Running filter_mc")
        run_parallel(samplesheet=args.input, cpg_bedfile=args.cpg_bedfile, output_dir=args.output_dir, ncpus=args.ncpus, has_header=args.has_header, blacklist_provided=args.blacklist_provided) ## add option

    elif args.command == 'filter_bam':
        print("Running filter_bam")
        print(args.input)
        run_bam_parallel(samplesheet=args.input, cpg_bedfile=args.cpg_bedfile, output_dir=args.output_dir, ncpus=args.ncpus, blacklist_provided=args.blacklist_provided)

    elif args.command == 'generate_blacklist':
        print("Running generate_blacklist")
        run_blacklist_parallel(samplesheet=args.input, cpg_bedfile=args.cpg_bedfile, output_dir=args.output_dir, ncpus=args.ncpus)

if __name__ == "__main__":
    main()
