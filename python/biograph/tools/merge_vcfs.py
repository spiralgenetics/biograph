"""
Runs the tools/process to merge multiple BioGraphAssembly vcfs together
"""
import sys
import shutil
import argparse
import tempfile
import multiprocessing

from biograph.utils import cmd_exe
import biograph.tools.log as log

def program_check():
    """
    This pipeline uses a bunch of external programs.
    Need to ensure every one is available
    """
    log.info("Checking for programs")
    req_progs = ["zcat", "cat", "sed", "bcftools", "vcf-sort", "bgzip", "tabix", "rtg", "parallel"]
    fail = False
    for prog in req_progs:
        ret = cmd_exe("which %s" % (prog))
        if ret.ret_code != 0:
            log.error("Couldn't verify program %s in environment: \n%s", prog, str(ret))
            fail = True
    if fail:
        exit(0)


def make_normalize_cmd(in_vcf, ref_dir, temp_dir, regions=""):
    """
    return the command to run bcftools norm over the input
    and the file name of the temporary output vcf created
    TODO: should we try to operate on raw vcfs, biograph directories, or either?
    """
    out_vcf = tempfile.mkstemp(suffix=".vcf.gz", dir=temp_dir)
    temp_dir = "-t " + temp_dir
    if regions != "":
        regions = "-R " + regions
    cat = ""
    if in_vcf.endswith(".gz"):
        cat = "zcat " + in_vcf
    else:
        cat = "cat " + in_vcf
    cmd = ("{extract_cmd} | sed 's/SVLEN=-/SVLEN=/' |"
           "bcftools norm -c s -m+any -f {reference}/source.fasta {reg} - | "
           "vcf-sort {tmp} | bgzip > {out_vcf} "
           "&& tabix -f -p vcf {out_vcf}")
    cmd = cmd.format(extract_cmd=cat, reference=ref_dir, out_vcf=out_vcf[1], tmp=temp_dir, reg=regions)
    log.debug("cmd: %s", cmd)
    return cmd, out_vcf[1]


def run_parallel(all_cmds, tmp_dir, threads=1):
    """
    Runs muiltple commands at once using Unix parallel
    """
    cmd_txt = tempfile.mkstemp(suffix=".txt", dir=tmp_dir)
    log.debug("parallel commands in %s", cmd_txt[1])
    with open(cmd_txt[1], 'w') as fout:
        fout.write("\n".join(all_cmds))

    log.info("Running %d norm commands over max %d threads", len(all_cmds), threads)
    norm_cmd = "parallel -j {jobs} < {cmd_txt}".format(cmd_txt=cmd_txt[1], jobs=threads)
    log.debug("cmd: %s", norm_cmd)
    ret = cmd_exe(norm_cmd)
    if ret.ret_code != 0:
        log.error("Failure to run `parallel` over norm commands (%s)", norm_cmd)
        log.error(str(ret))
        exit(ret.ret_code)
    log.debug(str(ret))


def run_rtgmerge(out_vcf, in_vcfs):
    """
    return the rtg command to merge the normalized input vcfs
    and write to out_vcf
    Todo: use --add-header
    """
    log.info("Running merge")
    cmd = "rtg vcfmerge -o {out_vcf} {in_vcfs}"
    cmd = cmd.format(out_vcf=out_vcf, in_vcfs=" ".join(in_vcfs))
    log.debug(cmd)
    ret = cmd_exe(cmd)
    if ret.ret_code != 0:
        log.error("Failure to run merge command %s", cmd)
        log.error(str(ret))
        exit(ret.ret_code)
    log.debug(str(ret))


def parse_args(args):
    """
    Parse all arguments
    """
    parser = argparse.ArgumentParser(prog="merge_vcfs", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcfs', metavar='VCF', type=str, nargs='+',
                        help='VCF files to merge (order is preserved)')
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Merged VCF output")
    parser.add_argument("-r", "--reference", type=str, required=True,
                        help="BioGraph reference directory")
    parser.add_argument("-R", "--regions", type=str, default="",
                        help="Bed file of regions to merge")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Maximum number of normalize jobs to create (%(default)s)")
    parser.add_argument("--tmp", type=str, default=tempfile.gettempdir(),
                        help="Temporary directory to write file (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    log.setup_logging(args.debug)
    if len(args.vcfs) == 1:
        log.error("Expected at least two vcfs to merge")
        exit(1)
    return args


def main(args):
    """
    Main modules
    """
    args = parse_args(args)
    program_check()

    # make temporary directory
    tmp_dir = tempfile.mkdtemp(dir=args.tmp)
    log.info("tmp_dir=%s", tmp_dir)

    all_cmds = []
    all_tmp_vcfs = []
    for vcf_fn in args.vcfs:
        cmd, tmp = make_normalize_cmd(vcf_fn, args.reference, tmp_dir, args.regions)
        all_cmds.append(cmd)
        all_tmp_vcfs.append(tmp)

    run_parallel(all_cmds, tmp_dir, args.threads)
    run_rtgmerge(args.output, all_tmp_vcfs)

    # Clean
    shutil.rmtree(tmp_dir)
    log.info("Finished")


if __name__ == '__main__':
    main(sys.argv[1:])
