import argparse
from debias.core import run_parallel

def main():
    parser = argparse.ArgumentParser(description='A tool to filter variant loci for (pseudo)bulk bisulfite sequencing methylation calling data')
    parser.add_argument('--input', required=True, help='csv file with columns: [sample,vcf_file,allc_file]. "sample" will be used as prefix for the output.')
    parser.add_argument('--cpg_bedfile', required=True, help="input CpG BED file")
    parser.add_argument('--output_dir', required=True, help='directory to write output')
    parser.add_argument('--ncpus', type=int, default=8, help='number of threads to use')
    parser.add_argument('--mc_file_type', default='allc', help='methylation calling file type. Default is ALLC file.')

    args = parser.parse_args()

    run_parallel(samplesheet=args.input, cpg_bedfile=args.cpg_bedfile, output_dir=args.output_dir, ncpus=args.ncpus, mc_file_type=args.mc_file_type)

if __name__ == "__main__":
    main()
