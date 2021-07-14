"""
Filter a vcf based on minium observations of an allele or minimum number of individuals
"""
import sys
import gzip
import argparse

from collections import Counter
import biograph.tools.log as log

def parse_args(args):
    """
    argument parser
    """
    parser = argparse.ArgumentParser(prog="freq_filter", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--variants", type=str, required=True,
                        help="VCF to filter (- for stdin)")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output file (stdout)")
    parser.add_argument("-m", "--min-observations", type=int, default=1,
                        help="Minimum number of observations (reads) across individuals (%(default)s)")
    parser.add_argument("-i", "--individuals", type=int, default=1,
                        help="Minimum number of individuals (%(default)s)")
    parser.add_argument("-d", "--depth_field", type=str, default="PAD",
                        help="VCF FORMAT field to use for depth.")
    args = parser.parse_args(args)
    log.setup_logging()

    return args


def filter_entry(line, min_obs, min_indiv):
    """
    Check the samples and filter
    returns a bit flag with
    0 - pass
    1 - too-few observations
    2 - too-few individuals
    """
    data = line.strip().split('\t')
    fmt = data[8].split(':')
    dp_idx = fmt.index("PAD")
    obs = 0
    n_indivs = 0
    for i in data[9:]:
        if i == ".":
            continue

        n_indivs += 1
        #This 1 indexing is a problem...
        obs += int(i.split(':')[dp_idx].split(',')[1])
    ret = 0
    if obs < min_obs:
        ret += 1
    if n_indivs < min_indiv:
        ret += 2
    return ret


def main(args):
    """
    Main runner
    """
    args = parse_args(args)
    fh = None
    # pylint has a false positive here
    # pylint: disable=useless-suppression
    if args.variants == "-":
        fh = sys.stdin
    elif args.variants.endswith(".gz"):
        fh = gzip.GzipFile(args.variants)
    else:
        fh = open(args.variants, 'r')
    out = None
    if args.output is None:
        out = sys.stdout
    else:
        out = open(args.output, 'w')
    stats = Counter()
    for line in fh:
        if line.startswith("#"):
            out.write(line)
            continue
        filt = filter_entry(line, args.min_observations, args.individuals)
        stats[filt] += 1
        if filt == 0:
            out.write(line)
    filtered = stats[1] + stats[2] + stats[3]
    total = stats[0] + filtered
    log.info("filtered %d of %d entries (%.2f)", filtered, total, filtered / float(total))
    log.info("Obs-filtered %d", stats[1])
    log.info("Ind-filtered %d", stats[2])
    log.info("Both-filtered %d", stats[3])

if __name__ == '__main__':
    main(sys.argv[1:])
