"""
Wrapper around bgbinary commands as well as pipelining

You'll notice the binary_cmd methods behave differently.
This is so they can be called with `biograph <cmd>` as well as
functionally as part of the `full_pipeline` method.

Function docstrings are used for their descriptions in __main__.py.
"""
import os
import json
import shlex
import shutil
import argparse
import tempfile
import multiprocessing

from collections import OrderedDict

import biograph
from biograph.utils import cmd_exe, format_timedelta
import biograph.tools.log as log

LOGGINGSETUP = False
# Where full_pipeline will store files
VCFDIR = "{bgdir}/analysis"
DISVCF = "{bgdir}/analysis/discovery.vcf.gz"
COVVCF = "{bgdir}/analysis/coverage.vcf.gz"
RESVCF = "{bgdir}/analysis/results.vcf.gz"
GRMDF = "{bgdir}/analysis/grm.jl"
CLSDF = "{bgdir}/analysis/df.jl"
TIMFN = "{bgdir}/qc/timings.json"

def vcf_compress(output):
    """
    Pipe stdout to vcf-sort/bgzip and then tabix
    """
    return f" | vcf-sort | bgzip > {output} && tabix {output}"

def tool_check(tools):
    """
    Check to ensure the environment is setup with the list of tools
    """
    cmd_successful = True
    for i in tools:
        r = cmd_exe(f"which {i}")
        if r.ret_code != 0:
            log.error(f"Couldn't find {i} in environment")
            cmd_successful = False
    if not cmd_successful:
        log.error("At least one program not found in PATH")
    return cmd_successful

def reference_cmd(clargs):
    ''' Build a BioGraph reference from FASTA '''
    log.setup_logging()

    # Sanity check for required reference args. 'bgbinary reference' args that are ignored
    # (eg. --tmp) are intentionally left out.
    parser = argparse.ArgumentParser(
        prog="reference",
        description="Prepare a fasta for use with BioGraph. The specified reference directory will be created and will contain the new reference and database files."
    )
    parser.add_argument("--in", help="Input reference fasta or fasta.gz", required=True)
    parser.add_argument("--refdir", help="Output reference directory", required=True)
    parser.add_argument("-f", "--force", action='store_true', help="Overwrite existing refdir")

    args = parser.parse_args(clargs)

    if isinstance(args, list):
        args = " ".join(args)

    extensions = ('amb', 'ann', 'bwt', 'pac', 'sa') # 'fai' ?

    # 'in' is a reserved word. Doh!
    fasta = vars(args)['in']
    refdir = args.refdir

    for ext in extensions:
        if not os.path.isfile(f"{fasta}.{ext}"):
            raise SystemExit(f"Can't find {fasta}.{ext}. The input fasta must be indexed with BWA.\nPlease run 'bwa index {fasta}' and try again.")

    # use the clargs as specified
    cmd = f"bgbinary reference {' '.join(clargs)}"
    ret = cmd_exe(cmd, cap_stderr=False)
    log.debug(ret)
    if ret.ret_code != 0:
        log.info(ret.ret_code)
        log.info(ret.stderr)
        log.info(ret.stdout)
        return ret

    for ext in extensions:
        shutil.copy(f"{fasta}.{ext}", f"{refdir}/source.fasta.{ext}")

    return ret

def check_create(args, pipe_args):
    """
    ensure the biograph, reference, and out are not specified in the create command
    this is called by pipelines
    pipe_args are the arguments parsed by full_pipeline that need to be checked
    return if checks passed
    """
    successful = True
    args = shlex.split(args)
    for chk in ["--tmp", "--out", "--ref", "--threads"]:
        if chk in args:
            log.error(f"Remove {chk} from --create parameters")
            successful = False

    if os.path.exists(pipe_args.biograph) and '--force' not in args:
        log.error(f"Remove {pipe_args.biograph} file or specify --force")
        successful = False
    return successful

def create_cmd(args, base_cmd="bgbinary create {args}", dryrun=False):
    ''' Convert reads to the BioGraph format '''
    if not LOGGINGSETUP:
        log.setup_logging()
    if isinstance(args, list):
        args = " ".join(args)
    cmd = base_cmd.format(args=args)
    if dryrun:
        return cmd
    ret = cmd_exe(cmd, cap_stderr=False, pipefail=True)
    log.info(ret)
    if not LOGGINGSETUP:
        log.info("Finished create")
    return ret

def check_discovery(args, pipe_args): # pylint:disable=unused-argument
    """
    ensure the biograph, reference, and out are not specified in the create command
    return if checks passed
    """
    args = shlex.split(args)
    successful = True
    for chk in ["--tmp", "--in", "--ref", "--out", "--threads"]:
        if chk in args:
            log.error(f"Remove {chk} from --discovery parameters")
            successful = False
    return successful

def discovery_cmd(args, base_cmd="bgbinary discovery {args}", dryrun=False):
    ''' Discover variants on a BioGraph vs. a reference '''
    if not LOGGINGSETUP:
        log.setup_logging()
    if isinstance(args, list):
        args = " ".join(args)
    cmd = base_cmd.format(args=args)
    if dryrun:
        return cmd
    ret = cmd_exe(cmd, cap_stderr=False, pipefail=True)
    log.debug(ret)
    if not LOGGINGSETUP:
        log.info("Finished discovery")
    return ret

def check_coverage(args, pipe_args): # pylint:disable=unused-argument
    """
    ensure the biograph, reference, and out are not specified in the create command
    """
    args = shlex.split(args)
    successful = True
    for chk in ["-b", "--biograph", "-v", "--variants", "-r", "--reference", "-o", "--output", "-d", "--dataframe"]:
        if chk in args:
            log.error(f"Remove {chk} from --coverage parameters")
            successful = False
    return successful

def coverage_cmd(args, pipe_args, dryrun=False, single_sample=True):
    '''
    Calculate coverage for VCF entries
    If single_sample, we'll be making the command output to inside the biograph,
    Otherwise, there will be a temporary file
    '''
    if single_sample:
        cmd = (f"biograph coverage -b {pipe_args.biograph} -r {pipe_args.reference} "
               f"-v {DISVCF} -d {CLSDF} --threads {pipe_args.threads} {args}") + vcf_compress(COVVCF)
    else:
        if "--ideal-insert" not in args:
            log.info("--ideal-insert not set in coverage args. Assuming 400.")
            args += " --ideal-insert 400"
        cmd = (f"biograph coverage -b {pipe_args.biograph} -r {pipe_args.reference} --placer-max-ambig 1 "
               f"-v {pipe_args.variants} --threads {pipe_args.threads} -d {pipe_args.tmpdf} {args}") + vcf_compress(pipe_args.tmpvcf)
    if dryrun:
        return cmd
    ret = cmd_exe(cmd, cap_stderr=False, pipefail=True)
    log.debug(ret)
    return ret

def check_grms(args, pipe_args):
    """
    ensure the biograph, reference, and out are not specified in the create command
    """
    args = shlex.split(args)
    successful = True
    for chk in ["-i", "--input", "-r", "--reference", "-o", "--output", "-t", "--threads"]:
        if chk in args:
            log.error(f"Remove {chk} in --grm parameters")
            successful = False

    if not os.path.exists(os.path.join(pipe_args.reference, "source.fasta.bwt")):
        log.error(f"Reference directory doen't have bwa index")
        successful = False
    return successful

def grm_cmd(args, pipe_args, dryrun=False):
    """
    Run truvari grm
    """
    cmd = f"truvari anno grm -r {pipe_args.reference}/source.fasta -i {COVVCF} -o {GRMDF} {args}"
    if dryrun:
        return cmd
    ret = cmd_exe(cmd, cap_stderr=False)
    log.debug(ret)
    return ret

def check_qual(args, pipe_args): # pylint:disable=unused-argument
    """
    ensure the biograph, reference, and out are not specified in the create command
    """
    args = shlex.split(args)
    successful = True
    for chk in ["-v", "--vcf", "-m", "--model", "-o", "--out", "--tmp", "-t", "--threads", "--df"]:
        if chk in args:
            log.error(f"Remove {chk} from --qual_classifier parameters")
            successful = False
    return successful

def qual_cmd(args, pipe_args, dryrun=False, single_sample=True):
    """
    run qual_classifier
    If single_sample, we'll be making the command output to inside the biograph,
    Otherwise, there will be a named file in pipe_args
    """
    if single_sample:
        cmd = (f"biograph qual_classifier --grm {GRMDF} -d {CLSDF} -v {COVVCF} -m {pipe_args.model} "
               f"--tmp {pipe_args.tmp} -t {pipe_args.threads} {args}") + vcf_compress(RESVCF)
    else:
        cmd = (f"biograph qual_classifier --clsf GT -d {pipe_args.tmpdf} -v {pipe_args.tmpvcf} -m {pipe_args.model} "
               f"--tmp {pipe_args.tmp} -t {pipe_args.threads} {args}") + vcf_compress(pipe_args.output)
    if dryrun:
        return cmd
    ret = cmd_exe(cmd, cap_stderr=False, pipefail=True)
    log.debug(ret)
    return ret

def gtcls_cmd(args, pipe_args, dryrun=False):
    """
    Run gt classifier
    """
    cmd = (f"biograph gt_classifier -d {pipe_args.tmpdf} -v {pipe_args.tmpvcf} -m {pipe_args.model} "
           f"--tmp {pipe_args.tmp} -t {pipe_args.threads} {args}") + vcf_compress(pipe_args.output)
    if dryrun:
        return cmd
    ret = cmd_exe(cmd, cap_stderr=False, pipefail=True)
    log.debug(ret)
    return ret

def clean_files(steps, keep):
    """
    Remove the intermediate vcfs/dataframes, analyze steps to make sure we're not removing the last vcf
    """
    def rm(path):
        if os.path.exists(path):
            os.remove(path)
    if "coverage" in steps:
        if keep not in ["all", "vcf"]:
            rm(DISVCF)
            rm(DISVCF + ".tbi")
    if "qual_classifier" in steps:
        if keep not in ["all", "vcf"]:
            rm(COVVCF)
            rm(COVVCF + ".tbi")
        if keep not in ["all", "jl"]:
            rm(GRMDF)
            rm(CLSDF)

def full_pipe_args(args):
    """
    Parse and check full_pipeline arguments
    """
    parser = argparse.ArgumentParser(prog="full_pipeline",
                                     description="Run the standard BioGraph pipeline for a single sample: create, discovery, coverage, grm, qual_classifier",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--biograph", metavar="BG", required=True,
                        help="BioGraph file (will be created if running the create step)")
    parser.add_argument("-r", "--reference", metavar="REF", required=True,
                        help="Reference genome folder, in BioGraph reference format")
    parser.add_argument("--reads",
                        help="Input reads for BioGraph create, if run (fastq, bam, cram)")
    parser.add_argument("-m", "--model", metavar="ML",
                        help="BioGraph classifier model for qual_classifier, if run")
    parser.add_argument("--tmp", default=tempfile.gettempdir(),
                        help="Temporary directory (%(default)s)")
    parser.add_argument("-t", "--threads", default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--keep", choices=["all", "jl", "vcf"],
                        help="Keep intermediate dataframes/VCFs or all")
    parser.add_argument("--dry-run", action="store_true",
                        help="Run a preflight check, print all steps to be run, and exit")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite any existing intermediary files")

    pipear = parser.add_argument_group(title="Pipeline Arguments",
                                       description="Control which section of the pipeline is run")
    pipear.add_argument("--resume", default=None,
                        help="Step at which to resume the pipeline")
    pipear.add_argument("--stop", default=None,
                        help="Last step of the pipeline to run")

    cmdpar = parser.add_argument_group(title="Individual step arguments",
                                       description="Specify any additional parameters to be passed to steps. Must be a single \"string\"")
    cmdpar.add_argument("--create", default="",
                        help="create parameters")
    cmdpar.add_argument("--discovery", default="",
                        help="discovery parameters")
    cmdpar.add_argument("--coverage", default="",
                        help="coverage parameters")
    cmdpar.add_argument("--grm", default="",
                        help="truvari grm parameters")
    cmdpar.add_argument("--qual_classifier", default="",
                        help="qual_classifier parameters")

    args = parser.parse_args(args)

    return args


def full_pipeline(args): # pylint: disable=too-many-statements
    ''' Run the full BioGraph single-sample pipeline '''
    global LOGGINGSETUP, VCFDIR, DISVCF, COVVCF, RESVCF, GRMDF, CLSDF, TIMFN # pylint: disable=global-statement

    args = full_pipe_args(args)
    os.environ['TMPDIR'] = args.tmp

    log.setup_logging()
    VCFDIR = VCFDIR.format(bgdir=args.biograph)
    DISVCF = DISVCF.format(bgdir=args.biograph)
    COVVCF = COVVCF.format(bgdir=args.biograph)
    RESVCF = RESVCF.format(bgdir=args.biograph)
    GRMDF = GRMDF.format(bgdir=args.biograph)
    CLSDF = CLSDF.format(bgdir=args.biograph)
    TIMFN = TIMFN.format(bgdir=args.biograph)

    LOGGINGSETUP = True
    # steps value indexes
    CHKIDX = 0 # the check function
    CMDIDX = 1 # the command function
    PARIDX = 2 # the parameters
    # the command to submit if its binary otherwise full_pipe_args
    PIPIDX = 3

    steps = OrderedDict()

    # Use the force
    if args.force:
        args.create += " --force"

    steps["create"] = [check_create, create_cmd, args.create]
    steps["create"].append((f"bgbinary create --reads {args.reads} --tmp {args.tmp} --ref {args.reference} "
                            f"--out {args.biograph} --threads {args.threads} {{args}}"))

    steps["discovery"] = [check_discovery, discovery_cmd, args.discovery]
    steps["discovery"].append((f"bgbinary discovery --tmp {args.tmp} --in {args.biograph} "
                               f"--ref {args.reference} --threads {args.threads} {{args}}") + vcf_compress(DISVCF))

    steps["coverage"] = [check_coverage, coverage_cmd, args.coverage, args]
    steps["grm"] = [check_grms, grm_cmd, args.grm, args]
    steps["qual_classifier"] = [check_qual, qual_cmd, args.qual_classifier, args]

    # controlling steps
    step_keys = list(steps.keys())
    if args.resume is not None:
        step_keys = step_keys[step_keys.index(args.resume):]
    if args.stop is not None:
        step_keys = step_keys[:step_keys.index(args.stop) + 1]

    # check parameters
    par_successful = True

    if "create" in step_keys and not args.reads:
        par_successful = False
        log.error("--reads required when running the create step")

    if "qual_classifier" in step_keys and not args.model:
        par_successful = False
        log.error("--model required when running the qual_classifier step")

    for cur_key in step_keys:
        if not steps[cur_key][CHKIDX](steps[cur_key][PARIDX], args):
            par_successful = False
            log.error(f"The {cur_key} step parameters have problems")

    # tool check - ensure the environment is setup
    cmd_successful = tool_check(["bgbinary", "vcf-sort", "bgzip", "truvari", "tabix"])
    if not cmd_successful:
        log.error("At least one program not found in PATH")

    if not cmd_successful or not par_successful:
        log.error("Setup errors detected. See above. Aborting.")
        exit(1)

    # reference check - just make sure it's valid
    biograph.Reference(args.reference)

    if args.dry_run:
        log.warning("--dry-run specified. The following steps would be run:\n")
        for cur_key in step_keys:
            print(steps[cur_key][CMDIDX](steps[cur_key][PARIDX], steps[cur_key][PIPIDX], dryrun=True))
        exit(0)

    # run commands
    if os.path.exists(TIMFN):
        timings = json.load(open(TIMFN))
    else:
        timings = OrderedDict()
    for cur_key in step_keys:
        log.info(f"Starting {cur_key}")
        ret = steps[cur_key][CMDIDX](steps[cur_key][PARIDX], steps[cur_key][PIPIDX])
        timings[cur_key] = format_timedelta(ret.run_time)
        if ret.ret_code != 0:
            log.error(f"Couldn't run {cur_key}")
            log.error("stdout:", ret.stdout)
            log.error("retcode:", ret.ret_code)
            json.dump(timings, open(TIMFN, 'w'))
            exit(ret.ret_code)
        log.info(f"Finished {cur_key}")

    # post processing
    clean_files(step_keys, args.keep)
    json.dump(timings, open(TIMFN, 'w'))
    log.info("Timings", json.dumps(timings, indent=4))
    log.info("Finished full_pipeline")

def squareoff_args(args):
    """
    Parse and check squareoff arguments
    """
    parser = argparse.ArgumentParser(prog="squareoff",
                                     description="Run the standard BioGraph pipeline for square-off of a sample",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--biograph", metavar="BG", required=True,
                        help="BioGraph file (will be created if running the create step)")
    parser.add_argument("-r", "--reference", metavar="REF", required=True,
                        help="Reference genome folder, in BioGraph reference format")
    parser.add_argument("-v", "--variants", metavar="VCF", required=True,
                        help="Input variants VCF file to square-off")
    parser.add_argument("-m", "--model", metavar="ML", required=True,
                        help="BioGraph classifier model for gt_classifier")
    parser.add_argument("-o", "--output", default="output.vcf.gz",
                        help="Output squared-off VCF file name (%(default)s)")
    parser.add_argument("--tmp", default=tempfile.gettempdir(),
                        help="Temporary directory (%(default)s)")
    parser.add_argument("-t", "--threads", default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--timings", default=None, type=str,
                        help="File to write a timings.json file (off)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Run a preflight check, print all steps to be run, and exit")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite any existing intermediary files")

    cmdpar = parser.add_argument_group(title="Individual step arguments",
                                       description="Specify any additional parameters to be passed to steps. Must be a single \"string\"")
    cmdpar.add_argument("--coverage", default="",
                        help="coverage parameters")
    cmdpar.add_argument("--gt_classifier", default="",
                        help="gt_classifier parameters")

    args = parser.parse_args(args)
    if not args.output.endswith(".vcf.gz"):
        log.error(f"--output must end with '.vcf.gz'")
        exit(1)
    return args


def squareoff(args): # pylint: disable=too-many-statements
    ''' Run the full BioGraph VCF square-off pipeline '''
    global LOGGINGSETUP # pylint: disable=global-statement

    args = squareoff_args(args)
    os.environ['TMPDIR'] = args.tmp

    log.setup_logging()
    temp_dir = tempfile.mkdtemp(prefix='bgtmp_', dir=args.tmp)
    args.tmpdf = os.path.join(temp_dir, "temp.df")
    args.tmpvcf = os.path.join(temp_dir, "temp.vcf.gz")

    LOGGINGSETUP = True
    # steps value indexes
    CHKIDX = 0 # the check function
    CMDIDX = 1 # the command function
    PARIDX = 2 # the parameters
    # the command to submit if its binary otherwise full_pipe_args
    PIPIDX = 3

    steps = OrderedDict()

    steps["coverage"] = [check_coverage, coverage_cmd, args.coverage, args]
    steps["gt_classifier"] = [check_qual, gtcls_cmd, args.gt_classifier, args]

    # check parameters
    par_successful = True

    for cur_key in steps:
        if not steps[cur_key][CHKIDX](steps[cur_key][PARIDX], args):
            par_successful = False
            log.error(f"The {cur_key} step parameters have problems")

    # tool check - ensure the environment is setup
    cmd_successful = tool_check(["vcf-sort", "bgzip", "tabix"])
    if not cmd_successful or not par_successful:
        log.error("Setup errors detected. See above. Aborting.")
        exit(1)

    # reference check - just make sure it's valid
    biograph.Reference(args.reference)

    if args.dry_run:
        log.warning("--dry-run specified. The following steps would be run:\n")
        for cur_key in steps:
            print(steps[cur_key][CMDIDX](steps[cur_key][PARIDX], steps[cur_key][PIPIDX], dryrun=True, single_sample=False))
        exit(0)

    # run commands
    if args.timings and os.path.exists(args.timings):
        timings = json.load(open(args.timings))
    else:
        timings = OrderedDict()
    for cur_key in steps:
        log.info(f"Starting {cur_key}")
        ret = steps[cur_key][CMDIDX](steps[cur_key][PARIDX], steps[cur_key][PIPIDX], single_sample=False)
        timings[cur_key] = format_timedelta(ret.run_time)
        if ret.ret_code != 0:
            log.error(f"Couldn't run {cur_key}")
            log.error("stdout:", ret.stdout)
            log.error("retcode:", ret.ret_code)
            if args.timings:
                json.dump(timings, open(args.timings, 'w'))
            exit(ret.ret_code)
        log.info(f"Finished {cur_key}")

    # post processing
    try:
        shutil.rmtree(temp_dir)
    except OSError:
        pass

    if args.timings:
        json.dump(timings, open(args.timings, 'w'))
    log.info("Timings", json.dumps(timings, indent=4))
    log.info("Finished squareoff")
