"""
Given a collection of VCFs, paste the sample names together

This tool assumes:
 - A single sample is present in each vcf.
 - The very first vcf has all of the header definitions needed for all vcfs to be valid.
 - Each sample name is unique.
 - Each vcf has identical loci in the exact same order.
 - Each vcf has identical FORMAT field counts/order.
 - The second to last FORMAT field is PGT. PGT will be swapped with GT.
 - The very first vcf's entry definitions (REF, INFO, etc) is all that's needed to be preserved.
"""
import sys
import gzip
import argparse

import biograph.tools.log as log

def parse_args(args):
    """
    Arg parser
    """
    parser = argparse.ArgumentParser(prog="vcf_sample_paste", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcfs', metavar='VCF', type=str, nargs='+',
                        help='VCF files to paste together')
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="Output file (stdout)")
    args = parser.parse_args(args)
    log.setup_logging()
    if len(args.vcfs) == 1:
        log.error("Can only paste 2 or more VCFs")
        exit(1)
    return args


def replace_gt(sample):
    """
    Swaps GT with PGT - this is a pretty dangerous assumption of [-2]
    """
    data = sample.split(':')
    data[0], data[-2] = data[-2], data[0]
    return ":".join(data)


def main(args):
    """
    Main runner
    """
    args = parse_args(args)
    if args.output is None:
        output = sys.stdout
    else:
        output = open(args.output, 'w')

    def get_fh(file_name):
        """
        Open regular or gzip file
        """
        if file_name.endswith(".gz"):
            return gzip.GzipFile(file_name)
        return open(file_name, 'r')

    file_handlers = []
    header_out = False
    all_samples = []
    for i in args.vcfs:
        fh = get_fh(i)
        samp = None
        while True:
            line = fh.readline()
            if line.startswith("#CHROM"):
                samp = line.strip().split('\t')[9]
                break
            if not header_out:
                output.write(line)
            if not line.startswith("#"):
                log.error("Problem parsing vcf %s\n%s", i, line)
                exit(1)

        header_out = True
        if samp in all_samples:
            log.error("Sample %s present twice in input vcfs (%s, %s)",
                      samp, fh.name, file_handlers[all_samples.index(samp)].name)
            exit(2)
        all_samples.append(samp)
        file_handlers.append(fh)
    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % ("\t".join(all_samples)))
    while True:
        line = file_handlers[0].readline()
        if line == "":
            break
        data = line.strip().split('\t')
        data[9] = replace_gt(data[9])
        output.write("\t".join(data[:10]))
        for fh in file_handlers[1:]:
            line = fh.readline()
            data = replace_gt(line.strip().split('\t')[9])
            output.write("\t" + data)
        output.write('\n')
    output.close()
    log.info("Finished")

if __name__ == '__main__':
    main(sys.argv[1:])
