"""
Discover variants from reference without aligning or genotyping
"""
import sys
import argparse
from datetime import datetime
import os.path

import biograph
import biograph.variants as bgexvar
from biograph.tools.log import setup_logging

def parse_args(args):
    """
    argument parser
    """
    parser = argparse.ArgumentParser(prog="discover", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output VCF file")
    parser.add_argument("-b", "--biograph", required=True,
                        help="Merged BioGraph file containing individuals")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference genome folder")
    parser.add_argument("-s", "--sample", default=None,
                        help="Sample in merged BioGraph to use (%(default)s)")
    parser.add_argument("--cache", default=False, action="store_true",
                        help="Attempt to cache as much as possible in RAM")
    parser.add_argument("--ref-map", default=None,
                        help="Cache reference map in given file")
    parser.add_argument("--filter", default=False, action="store_true",
                        help="Do some simple filtering before outputting")
    args = parser.parse_args(args)
    setup_logging()

    return args

def write_header(fh, args):
    """Writes out the vcf header"""
    fh.write("""##fileformat=VCFv4.1
##fileDate={date}
##source=Spiral Genetics BioGraph v{version} {cmdline}
##reference=file://{reference}
##INFO=<ID=AID,Number=.,Type=Integer,Description="Assembly IDs used in constructing this variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
""".format(date=datetime.now().strftime('%Y%m%d'),
           version=biograph.version(),
           cmdline=" ".join(sys.argv[:]),
           reference=args.reference))

def show_progress(progress):
    """Updates user with current progress."""
    print(f"...{progress*100:.2f}")
    sys.stdout.flush()

def write_asm(fh, ref, scaffold_name, a):
    """Encodes an assembly in VCF format and writes it out"""
    try:
        if a.left_offset == a.right_offset:
            refseq = ""
        else:
            ref_range = ref.make_range(scaffold_name, a.left_offset, a.right_offset, False)
            refseq = str(ref_range.sequence)
        if not a.seq or not refseq:
            position = a.left_offset - 1
            last_base = str(ref.make_range(scaffold_name, position, position + 1, False).sequence)
            refseq = last_base + refseq
            altseq = last_base + str(a.seq)
        else:
            position = a.left_offset
            altseq = a.seq
    except RuntimeError as e:
        print(f"{str(e)}: {str(a)}")
        return
    if not refseq:
        refseq = "."
    if not altseq:
        altseq = "."
    info = f"AID={a.assembly_id}"
    fmt = "GT"
    sample = f"0/1"
    fh.write(f"{scaffold_name}\t{position+1}\t.\t{refseq}\t{altseq}\t.\tPASS\t{info}\t{fmt}\t{sample}\n")

def main(args):
    ''' Discover variants without aligning or genotyping (experimental) '''
    # NOTE: ^^^ this ^^^ docstring is used for the command description in __main__.py

    args = parse_args(args)
    bg = biograph.BioGraph(args.biograph,
                           biograph.CacheStrategy.RAM if args.cache
                           else biograph.CacheStrategy.MMAPCACHE)
    ref = biograph.Reference(args.reference)
    if args.sample is None:
        if len(bg.metadata.samples) == 1:
            args.sample = next(iter(bg.metadata.samples))
    elif args.sample not in bg.metadata.samples:
        raise KeyError("Sample %s not present in BioGraph" % args.sample)
    rm = bg.open_readmap(args.sample)

    if args.ref_map:
        if not os.path.exists(args.ref_map):
            print(f"Generating refmap to {args.ref_map}")
            bgexvar.RefMap.generate_and_save(bg.seqset, ref, args.ref_map, show_progress)
            print("Done generating refmap")
        print(f"Loading refmap from {args.ref_map}")
        rmap = bgexvar.RefMap.load(bg.seqset, ref, args.ref_map)
        print("Done loading refmap")
    else:
        print("Generating refmap in memory")
        rmap = bgexvar.RefMap.generate(bg.seqset, ref, show_progress)
        print("Done generating refmap")

    with open(args.output, 'w') as fh:
        write_header(fh, args)

        for scaffold_name, scaffold_len in ref.scaffold_lens.items():
            fh.write(f"##contig=<ID={scaffold_name},length={scaffold_len}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        disc = bgexvar.ParallelDiscover(rm, ref, rmap)
        disc.add_scaffold("chr13")
        #        disc.add_entire_reference()
        def asm_to_fh(scaffold_name, asm):
            write_asm(fh, ref, scaffold_name, asm)
        print("Starting assembly")
        disc.assemble(asm_to_fh, show_progress, do_filtering=args.filter)
        print("Assembly complete")

if __name__ == '__main__':
    main(sys.argv[1:])
