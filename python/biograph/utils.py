"""
Utilities for biograph sdk

This library exports all of the functionality of biograph, and also provides
the following helper functions:

  * find_region_variants(), find_breakpoint_variants() automates the process of finding
    anchors, assembling between them, and finding the correct coverage for the resulting
    assembly.

  * Assembly is a container class for assembly results. The find_variants functions
    returns Assembly object(s).

  * visualize() prints a text visualization of Assembly objects.

  * genotyper - given all coverage and alternate coverage depths, use a bayesian
    estimation at the genotype
"""
import codecs
import gzip
import logging
import math
import multiprocessing
import os
import signal
import subprocess
import sys
import time
import traceback
import re

from collections import namedtuple, defaultdict
from datetime import timedelta, datetime
from io import BytesIO

# pylint: disable=too-few-public-methods,too-many-lines

class Assembly:

    """
        Simple container class for assembly results including:
            ref_range: a ReferenceRange for the containing region
            variants: a list of Variants found
            coverage: a list of reference coverage depths for each base in the range
    """

    def __init__(self, ref_range=None, variants=None, coverage=None):
        self.data = []

        self.ref_range = ref_range
        self.variants = variants
        self.coverage = coverage

    def to_vcf(self, header=False):
        """
        Turns this assembly into a StringIO vcf.
        if header=True, it will first write a VCF header to the return (do this on the first call)

        """
        ret = BytesIO()
        if header:
            ret.write((
                '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples">\n'
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV">\n'
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural Variant Type">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">\n'
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample Depth">\n'
                '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'))

        for var in self.variants:
            if var.left_contig != var.right_contig:
                raise ValueError("Assembly variant to vcf does not currently support translocations")
            chrom = var.left_contig
            start = var.left_position
            end = var.right_position
            ref_seq = str(self.ref_range.sequence[start - self.ref_range.start + 1: end - self.ref_range.start])
            alt_seq = str(var.assembly_sequence[var.assembly_begin:var.assembly_end])
            if len(ref_seq) != len(alt_seq):
                anchor_base = self.ref_range.sequence[start - self.ref_range.start]
            else:
                anchor_base = ""
            alt_depth = sum(var.depths) / len(var.depths)
            rcov = self.coverage[start - self.ref_range.start + 1: end - self.ref_range.start]
            ref_depth = sum(rcov) / len(rcov)
            genotype, genoqual = genotyper(alt_depth + ref_depth, alt_depth)
            svtype = ""
            svlen = ""
            if var.is_structural:
                svtype = "SVTYPE=DEL;" if len(ref_seq) > len(alt_seq) else "SVTYPE=INS;"
                svlen = "SVLEN=%d;" % (len(ref_seq) - len(alt_seq))
            ret.write("{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tNS=1;{svtype}{svlen}\tGT:GQ:DP:AD\t{gt}:{gq:.2f}:{dp}:{rd},{ad}\n".format(
                chrom=chrom, pos=start, ref=anchor_base + ref_seq, alt=anchor_base + alt_seq, svtype=svtype, svlen=svlen, gt=genotype,
                gq=genoqual, dp=alt_depth + ref_depth, rd=ref_depth, ad=alt_depth))
        ret.seek(0)
        return ret

# pylint: disable=too-many-locals


def find_breakpoint_variants(my_bg, ref, supercontig, start, end,
                             min_overlap=70, max_anchors=10000, max_steps=100000,
                             skip_ambiguous=False, buf_len=300):
    """
    Find anchors and return variants + true reference coverage for the
    given region. This includes coverage for reads that exactly match
    reference, as well as the reference component of variant assemblies.
    Results are returned as an Assembly object

    The method works by finding reads that do not perfectly match the
    reference on the 5' side (a.k.a. start) of the variant and the 3' side
    (a.k.a. end). These reads are called 'anchors'.

    buf_len specifies how much of a region around the start/end to look for anchors.
    For example, if we look for
        supercontig = "chr1"
        start = 5000
        end = 10000
        buf_len = 500

    We assemble between all anchors found within chr1:4500-5500 and chr1:9500-10500.

    If start and end are within buf_len of each other, we choose the boundary between
    5'side anchors and 3'side anchors as the midpoint between start and end.
    For example, if we look for
        supercontig = "chr1"
        start = 8000
        end = 8100
        buf_len = 200

    We assemble between all anchors found within chr1:7800-8050 and chr1:8050-8300

    If the start/end region do not lie in the same supercontig, or the start - buf_len
    is less than zero, an exception is raised.
    See biograph.Reference.make_range for details.

    This method is faster than find_region_variants and can help finding events like
    large spanning deletions, however it cannot find inversions and requires an
    approximate knowledge of breakpoint coordinates (e.g. breakpoints will
    be within +- buf_len of the variant.
    """
    if start >= end:
        raise RuntimeError("start must be < end")

    # find_ranges would need to behave differently all_variants = []
    fwd_start, rev_start = 0, 0
    fwd_end, rev_end = 0, 0
    if start + buf_len >= end - buf_len:
        mid = int((end + start) / 2)
        fwd_start = start - buf_len
        fwd_end = mid
        rev_start = mid
        rev_end = end + buf_len
    else:
        fwd_start = start - buf_len
        fwd_end = start + buf_len
        rev_start = end - buf_len
        rev_end = end + buf_len

    ref_range = ref.make_range(supercontig, fwd_start, rev_end)

    # this could also be fwd = fwd_start, fwd_end, False and rev = rev_start,
    # rev_end, True
    from biograph.internal import find_anchors, assemble
    fwd = find_anchors(my_bg, ref.make_range(
        supercontig, rev_start, rev_end), True, min_overlap, max_anchors)
    rev = find_anchors(my_bg, ref.make_range(
        supercontig, fwd_start, fwd_end), False, min_overlap, max_anchors)

    results = assemble(
        fwd, rev, min_overlap, max_steps, skip_ambiguous, my_bg.readmap)
    ref_range = ref.make_range(supercontig, start - buf_len, end + buf_len)
    # Add in reference coverage for anchors. Start with reference coverage.
    fixed_coverage = my_bg.seq_coverage(ref_range.sequence)

    # Add hom-reference object (no variants)
    if not results:
        # all_variants.append(Assembly(ref_range=ref_range, variants=[],
        # coverage=fixed_coverage))
        return Assembly(ref_range=ref_range, variants=[], coverage=fixed_coverage)

    # Coverage entries are of the format:
    #     ['scaffold', position, [25,25,26,26,26...]]
    for cov in results[1]:
        if cov[0] != ref_range.scaffold:
            continue

        for i in range(len(cov[2])):
            mod_pos = cov[1] + i - ref_range.start
            if 0 <= mod_pos < ref_range.size:
                fixed_coverage[mod_pos] += cov[2][i]

    # all_variants.append(Assembly(ref_range=ref_range,
    # variants=sorted(results[0]), coverage=fixed_coverage))
    return Assembly(ref_range=ref_range, variants=sorted(results[0]), coverage=fixed_coverage)


def find_region_variants(my_bg, ref, supercontig, start, end,
                         min_overlap=70, max_anchors=10000, max_steps=100000,
                         skip_ambiguous=True):
    """
    Find anchors and return variants + true reference coverage for the
    given region. This includes coverage for reads that exactly match
    reference, as well as the reference component of variant assemblies.

    Returns a list of Assembly objects for each supercontig in the region.
    If no variants are found in the assembly, it'll have an empty list around
    the variants. If no assemblies are found, an empty list will be returned.

    This works differently than find_breakpoint_variants in that it does not try to divide
    the anchors. Instead, every anchor is compared to every other anchor within
    the region.
    This creates more false-positives and runs more slowly than find_breakpoint_variants, but
    can be more sensitive in some cases.
    """
    if start >= end:
        raise RuntimeError("start must be < end")

    all_variants = []

    # For all contigs within the region...
    for ref_range in ref.find_ranges(supercontig, start, end):

        # Collect the forward and reverse anchors
        from biograph.internal import find_anchors, assemble
        fwd = find_anchors(
            my_bg, ref_range, True, min_overlap, max_anchors)
        rev = find_anchors(
            my_bg, ref_range, False, min_overlap, max_anchors)

        # Assemble between anchors. Returns [variants, anchor coverage]
        results = assemble(
            fwd, rev, min_overlap, max_steps, skip_ambiguous, my_bg.readmap)

        # Add in reference coverage for anchors. Start with reference coverage.
        fixed_coverage = my_bg.seq_coverage(ref_range.sequence)

        # Add hom-reference object (no variants)
        if not results:
            all_variants.append(
                Assembly(ref_range=ref_range, variants=[], coverage=fixed_coverage))
            continue

        # Coverage entries are of the format:
        #     ['scaffold', position, [25,25,26,26,26...]]
        for cov in results[1]:
            if cov[0] != ref_range.scaffold:
                continue

            for i in range(len(cov[2])):
                mod_pos = cov[1] + i - ref_range.start
                if 0 <= mod_pos < ref_range.size:
                    fixed_coverage[mod_pos] += cov[2][i]

        all_variants.append(
            Assembly(ref_range=ref_range, variants=sorted(results[0]), coverage=fixed_coverage))

    return all_variants

# This was produced as a "quick demo" and could use some refactoring.
# pylint: disable=too-many-statements


def visualize(assembly, min_size=0, use_ascii=True):
    '''
        Print a textual variant visualization from the assemblies generated by find_variant methods
    '''
    if not assembly.variants:
        print("No variants found.")
        return

    # Filter variants by min-size
    filtered_variants = [var for var in assembly.variants if len(
        var.sequence) >= min_size or var.range.size >= min_size]

    if not filtered_variants:
        print("No variants passed filtering (out of {count} found).".format(count=len(assembly.variants)))
        return

    # Filtering complete. Display code below.

    # Save stdout (keep Jupyter notebook happy)
    oldstdout = sys.stdout

    # UTF8 support
    if use_ascii:
        p_horz = "-"
        p_hdot = ">"
        p_vert = "|"
        p_splt = "|"
        p_top = " "
        p_bot = " "
        p_svs = "|"
        p_sve = "|"
    else:
        utf8writer = codecs.getwriter('utf8')
        sys.stdout = utf8writer(sys.stdout)

        p_horz = chr(0x2500)
        p_hdot = chr(0x2504)
        p_vert = chr(0x2502)
        p_splt = chr(0x251c)
        p_top = chr(0x256e)
        p_bot = chr(0x256f)
        p_svs = chr(0x2570)
        p_sve = chr(0x256d)

    # Walk vars once to decide how much vertical room is needed for various
    # reference positions
    ref_seq = assembly.ref_range.sequence
    ref_size = assembly.ref_range.size
    ref_contig = assembly.ref_range.scaffold
    ref_start = assembly.ref_range.start
    ref_ctx = 3
    ref_vspace = [0 for _ in range(ref_size)]
    sv_list = {}
    for var in filtered_variants:
        if var.is_structural:
            left_vstart = var.left_position + \
                (+1 if var.left_forward else -1) - ref_start
            right_vstart = var.right_position + \
                (-1 if var.right_forward else +1) - ref_start
            if var.left_contig == ref_contig and 0 <= left_vstart < ref_size:
                for i in range(max(0, left_vstart - ref_ctx), min(ref_size, left_vstart + ref_ctx)):
                    ref_vspace[i] = max(ref_vspace[i], 1)
                sv_list[left_vstart] = sv_list.get(left_vstart, []) + [var]
            if var.right_contig == ref_contig and 0 <= right_vstart < ref_size:
                for i in range(max(0, right_vstart - ref_ctx), min(ref_size, right_vstart + ref_ctx)):
                    ref_vspace[i] = max(ref_vspace[i], 1)
                sv_list[right_vstart] = sv_list.get(
                    right_vstart - 1, []) + [var.flip()]
            continue
        # Compute [start, end) of reference region
        spos = var.left_position - ref_start + 1
        epos = var.right_position - ref_start
        # Get size of reference region
        vref_size = epos - spos
        # Get size of assembled region
        vasm_size = len(var.sequence)
        # Make sure we display reference across the whole variant
        for i in range(max(0, spos - ref_ctx), min(ref_size, epos + ref_ctx)):
            ref_vspace[i] = max(ref_vspace[i], 1)
        # If we need more vertical space to fit the variant, add it
        if vasm_size > vref_size:
            ref_vspace[spos - 1] = max(
                ref_vspace[spos - 1], vasm_size - vref_size + 1)

    # Map vlocs to ref-space
    vloc_to_ref = []  # vloc -> ref map initialy empty
    ref_to_vloc = [-1 for i in range(ref_size)]
                                     # Initially, no ref location has a vloc
    cur_in_ref = True  # Start in reference to make initial ... if we don't start at the top
    for i, _ in enumerate(ref_vspace):
        if ref_vspace[i] == 0:  # If this area of reference isn't needed
            if cur_in_ref:  # If we were in a reference area
                vloc_to_ref += [-2]  # Add a reference break
                cur_in_ref = False
            continue
        cur_in_ref = True  # We are now in reference
        vloc_to_ref += [i]  # Add in the offical position
        ref_to_vloc[i] = len(vloc_to_ref) - 1  # Make the reverse mapping
        for j in range(ref_vspace[i] - 1):  # Add in 'inserted' bases
            vloc_to_ref += [-1]

    # Map variants to their positions
    vsize = len(vloc_to_ref)
    v_to_hsize = [0 for i in range(vsize)]
    hv_to_base = {}
    hv_to_depth = {}
    hv_start = {}
    hv_end = {}
    maxlen = 12

    for var in filtered_variants:
        if var.is_structural:
            continue  # Ignore SV's for now
        # Compute start and end within vertical space
        spos = ref_to_vloc[var.left_position - ref_start] + 1
        epos = ref_to_vloc[var.right_position - ref_start]
        # Find an unused hpos
        hpos = 0
        for j in range(spos, epos):
            if v_to_hsize[j] > hpos:
                hpos += 1
        # Compute distance in vloc space and assembly size
        vlen = epos - spos
        alen = len(var.sequence)
        # Make sure things fit
        if alen > vlen:
            print(ref_vspace[var.left_position - ref_start])
            print(ref_vspace[var.right_position - ref_start])
            print(var.sequence)
            raise RuntimeError("Error in variant positioning")
        # Make things wider as needed
        for j in range(spos, epos):
            v_to_hsize[j] = max(v_to_hsize[j], hpos + 1)
        # Add in 'start' and 'end' with depth at split
        hv_start[(spos, hpos)] = var.depths[0]
        hv_end[(epos - 1, hpos)] = var.depths[-1]
        # Add in sequence data
        for j in range(spos, epos):
            sloc = j - spos  # Compute location in sequence
            if sloc < len(var.sequence):
                # Add actual bases
                hv_to_base[(j, hpos)] = var.sequence[sloc]
                hv_to_depth[(j, hpos)] = var.depths[sloc + 1] if var.depths[sloc + 1] >= 0 else ""
            else:
                # Add empty bases
                hv_to_base[(j, hpos)] = '.'
                hv_to_depth[(j, hpos)] = 0

    # Print everything
    column_active = {}
    for i in range(vsize):
        ref_loc = vloc_to_ref[i]
        more = u""
        prev_ref = assembly.coverage[max(0, vloc_to_ref[max(0, i - 1)])]
        if ref_loc in sv_list:
            for var in sv_list[ref_loc]:
                if not var.left_forward:
                    continue
                prev_ref -= var.depths[0]
                pline = u"{:2}  {}{}{}  {:2}  {}".format(
                    min(99, prev_ref),
                    p_splt, (3 * p_horz), p_top,
                    min(99, var.depths[0]),
                    var.sequence if len(var.sequence) <= maxlen else (
                        str(var.sequence)[:maxlen] + '...[' + str(len(var.sequence)) + ']')
                )
                rline = u"    {}   {}{}{} {}:{} {}".format(
                    p_vert, p_svs, (
                        min(len(var.sequence), maxlen) + 6) * p_horz, p_hdot,
                    var.right_contig, var.right_position, p_top if var.right_forward else p_bot
                )
                print(u"                            " + pline)
                print(u"                            " + rline)
        for j in range(v_to_hsize[i]):
            if (i, j) in hv_start:
                prev_ref -= hv_start[(i, j)]
                pline = u"{:2}  {}".format(min(99, prev_ref), p_splt)
                for k in range(3 + j * 8):
                    pline += p_horz
                pline += p_top
                pline += u"  {:2}".format(min(99, hv_start[(i, j)]),)
                for k in range(j + 1, v_to_hsize[i]):
                    if column_active.get(k, 0):
                        pline += u"   {}    ".format(p_vert)
                    else:
                        pline += u"        "
                print(u"                            " + pline)
                column_active[j] = 1
        for j in range(v_to_hsize[i]):
            if (i, j) in hv_to_base:
                if hv_to_base[(i, j)] == '.':
                    more += u"  {}     ".format(p_vert)
                else:
                    more += u"  {}{} {:2} ".format(
                        p_vert,
                        hv_to_base[(i, j)],
                        min(99, hv_to_depth[(i, j)]))
            else:
                more += u"        "
        if ref_loc >= 0:
            print(u"{:<12}:{:>12}   {:2} {:.1}{:.1} {}".format(
                ref_contig,
                ref_start + ref_loc,
                min(99, assembly.coverage[ref_loc]),
                ref_seq[ref_loc],
                p_vert,
                more))
        elif ref_loc == -1:
            print(u"                                {} {}".format(p_vert, more))
        else:
            print(u"                                .")
            print(u"                                .")
            print(u"                                .")
        # Compute coverage sum at the next point
        prev_ref = assembly.coverage[
            max(0, vloc_to_ref[min(vsize - 1, i + 1)])]
        for j in range(v_to_hsize[i]):
            if (i, j) in hv_end:
                prev_ref -= hv_end[(i, j)]
        if ref_loc in sv_list:
            for var in sv_list[ref_loc]:
                if var.left_forward:
                    continue
                prev_ref -= var.depths[0]
        # End SV's
        if ref_loc in sv_list:
            for var in sv_list[ref_loc]:
                if var.left_forward:
                    continue
                pline = u"{:2}  {}{}{}  {:2} {}".format(
                    min(99, prev_ref) if prev_ref >= 0 else "",
                    p_splt, (3 * p_horz), p_bot,
                    min(99, var.depths[0]),
                    var.sequence.rev_comp() if len(var.sequence) <= maxlen else (
                        str(var.sequence.rev_comp())[:maxlen] + '...[' + str(len(var.sequence)) + ']')
                )
                rline = u"    {}   {}{}{} {}:{} {}".format(
                    p_vert, p_sve, (
                        min(len(var.sequence), maxlen) + 6) * p_horz, p_hdot,
                    var.right_contig, var.right_position, p_top if var.right_forward else p_bot)
                prev_ref += var.depths[0]
                print(u"                            " + rline)
                print(u"                            " + pline)
        # End everything else
        for j in range(v_to_hsize[i]):
            if (i, j) in hv_end:
                pline = u"{:2}  {}".format(min(99, prev_ref) if prev_ref >= 0 else "", p_splt)
                for k in range(3 + j * 8):
                    pline += p_horz
                pline += p_bot
                pline += u"  {:2}".format(min(99, hv_end[(i, j)]),)
                for k in range(j + 1, v_to_hsize[i]):
                    if column_active.get(k, 0):
                        pline += u"   {}    ".format(p_vert)
                    else:
                        pline += u"        "
                print(u"                            " + pline)
                prev_ref += hv_end[(i, j)]
                del column_active[j]

    # Restore stdout (keep Jupyter notebook happy)
    sys.stdout = oldstdout

#{{{ http://code.activestate.com/recipes/511478/ (r1)
# pylint: disable=invalid-name


def genotyper(totCov, altCov, priors=None):
    """
    Given total coverage and altCoverage, try to calculate how many copies
    of the alt are at the position (ref/het/hom) and a quality score
    returns two lists
    - probabilities of 0, 1, or two copies of the allele at the location
    - phred-scaled quality scores of those probs
    """
    #We have no information.. should give up
    if totCov == 0:
        return None

    # previously had avgCov
    if priors is None:
        priors = [0.05, 0.5, 0.95]

    # if len(priors) != 3: # raise exception?

    def log_choose(n, k):
        """ swap for efficiency if k is more than half of n """
        r = 0.0
        if k * 2 > n:
            k = n - k

        for d in range(1, k + 1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1

        return r

    total = totCov  # refCoverage + altCoverage if avgCov is None else avgCov
    alt = altCov  # int(spot.tags["szCount"])
    non_alt = total - alt

    gtList = []

    comb = log_choose(total, alt)
    for p_alt in priors:
        gtList.append(comb + alt * math.log(p_alt, 10) + non_alt * math.log(1 - p_alt, 10))

    return gtList

def load_regions(region_bed, region_strings):
    """
    Loads/parses region arguments in place
    """
    ret = []
    if region_bed is not None:
        raise RuntimeError("region bed no longer supported.")
    if region_strings is not None:
        bed_entry = namedtuple("Region", "chrom start end")
        for x in region_strings:
            try:
                chrom, coord = x.split(':')
                start, end = [int(y) for y in coord.split('-')]
                ret.append(bed_entry(chrom, start, end))
            except ValueError:
                sys.stderr.write("Malformed region string - %s\n" % x)
                sys.exit(1)
    return ret if ret != [] else None

class Alarm(Exception):

    """ Alarm Class for command timeouts """
    pass # pylint: disable=unnecessary-pass

def alarm_handler(signum, frame=None):
    """ Alarm handler for command timeouts """
    raise Alarm

def format_timedelta(td):
    """
    format timedelta to a string of hh:mm:ss
    """
    tot_seconds = (3600*24) * td.days + td.seconds
    hours, remainder = divmod(tot_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return '{:02}:{:02}:{:02}'.format(int(hours), int(minutes), int(seconds))

def cmd_exe(cmd, timeout=-1, cap_stderr=True, pipefail=False):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    returns namedtuple(ret_code, stdout, stderr, run_time)
    where ret_code is the exit code for the command executed
    stdout/err is the Standard Output Error from the command
    and runtime is hh:mm:ss of the execution time

    cap_stderr will capture the stderr and return it as part of the
    returned cmd_result. Otherwise, stderr will be streamed through.
    set pipefail=True if you're using pipes
    """
    cmd_result = namedtuple("cmd_result", "ret_code stdout stderr run_time")
    t_start = time.time()
    stderr = subprocess.PIPE if cap_stderr else None
    if pipefail:
        cmd = f"set -o pipefail; {cmd}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stdin=sys.stdin, stderr=stderr, close_fds=True,
                            start_new_session=True, executable="/bin/bash")
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout * 60))
    try:
        stdoutVal, stderrVal = proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d"), timeout)
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return cmd_result(214, None, None, timedelta(seconds=time.time() - t_start))
    except KeyboardInterrupt:
        logging.error("KeyboardInterrupt on cmd %s", cmd)
        os.killpg(proc.pid, signal.SIGKILL)
        proc.kill()
        try:
            sys.exit(214)
        except SystemExit:
            os._exit(214) # pylint: disable=protected-access

    stdoutVal = bytes.decode(stdoutVal)
    retCode = proc.returncode
    ret = cmd_result(retCode, stdoutVal, stderrVal, timedelta(seconds=time.time() - t_start))
    return ret

def get_opener(filename):
    """
    Auto gzip decoder. Returns an open filehandle in binary mode.
    """

    # Can't rewind stdin. Assume plain text.
    if filename == '/dev/stdin':
        return open(filename, 'rb')

    with open(filename, 'rb') as f:
        if f.read(2) == b'\x1f\x8b':
            return gzip.open(filename, 'rb')

    return open(filename, 'rb')

# File for handling pedigree information
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype - See more at:
# http://gatkforums.broadinstitute.org/gatk/discussion/7696/pedigree-ped-files

class Pedigree(dict):

    """
    Parses a pedigree file and allows different views
    """

    def __init__(self, file_name):
        super(Pedigree, self).__init__()
        self.samples = self.keys
        self.families = defaultdict(list)

        with open(file_name, 'r') as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                data = line.strip().split('\t')
                # check the data length
                n_ped = _PedSample(*data[:5], phenotype=data[5:])
                self.families[n_ped.fam_id].append(n_ped)
                if n_ped.ind_id in self:
                    raise KeyError("Duplicate Individual Id %s" % n_ped.ind_id)
                self[n_ped.ind_id] = n_ped

            # Give parents a presence in the ped, even if they didn't have a line
            for ind in self.values():
                if ind.pat_id not in self and ind.pat_id != "0":
                    self[ind.pat_id] = _PedSample(ind.fam_id, ind.pat_id, "0", "0", "1", "0")
                if ind.mat_id not in self and ind.mat_id != "0":
                    self[ind.mat_id] = _PedSample(ind.fam_id, ind.mat_id, "0", "0", "2", "0")
            # Set parent's offspring
            for n_ped in self.values():
                if n_ped.pat_id in self:
                    self[n_ped.pat_id].offspring.append(n_ped)
                    n_ped.father = self[n_ped.pat_id]
                if n_ped.mat_id in self:
                    self[n_ped.mat_id].offspring.append(n_ped)
                    n_ped.mother = self[n_ped.mat_id]

    def filter(self, inc_fam=None, exc_fam=None, inc_indiv=None, exc_indiv=None):
        """
        Exclude anything that's exc in the pedigree.
        Include only anything that's inc in the pedigree.
        """
        if inc_fam is not None:
            for i in self.keys():
                if self[i].fam_id not in inc_fam:
                    del self[i]
        if inc_indiv is not None:
            for i in self.keys():
                if self[i].ind_id not in inc_indiv:
                    del self[i]
        if exc_fam is not None:
            for i in self.keys():
                if self[i].fam_id in exc_fam:
                    del self[i]
        if exc_indiv is not None:
            for i in self.keys():
                if self[i].ind_id in exc_indiv:
                    del self[i]

    def all_male(self):
        """
        Returns all male individuals
        """
        for i in self:
            if self[i].sex == "1":
                yield self[i]

    def all_female(self):
        """
        Returns all female individuals
        """
        for i in self:
            if self[i].sex == "2":
                yield self[i]

    def all_affected(self):
        """
        Returns all affected individuals
        """
        for i in self:
            if self[i].phenotype == "2":
                yield self[i]

    def all_unaffected(self):
        """
        Returns all unaffected individuals
        """
        for i in self:
            if self[i].phenotype == "1":
                yield self[i]

    def get_siblings(self, indiv):
        """
        Returns the siblings of an individual
        """
        for i in self:
            if self[i].pat_id == self[indiv].pat_id or self[i].mat_id == self[indiv].mat_id:
                yield self[i]

    def get_trio_probands(self):
        """
        Yields _PedSample probands that are part of a trio i.e. niether parent is 0
        """
        for indiv in self.values():
            if indiv.mat_id != '0'and indiv.pat_id != '0':
                yield indiv

    def get_quad_probands(self):
        """
        Yields _PedSample proband tuples that are part of an exact quad.
        """
        for fam in self.families:
            already_yielded = {}
            for indiv in self.families[fam]:
                if indiv.ind_id in already_yielded:
                    continue
                if indiv.mat_id != "0" and indiv.pat_id != "0":
                    siblings = set(self[indiv.mat_id].offspring).intersection(set(self[indiv.pat_id].offspring))
                    if len(siblings) == 2:
                        yield list(siblings)
                    for sib in siblings:
                        if indiv != sib:
                            already_yielded[sib.ind_id] = 1
                            yield (indiv, sib)


class _PedSample():

    """
    An individual in a pedigree
    Family ID
    Individual ID
    Paternal ID
    Maternal ID
    Sex (1=male; 2=female; other=unknown)
    Phenotype
    """
    def __init__(self, fam_id, ind_id, pat_id, mat_id, sex, phenotype):
        self.fam_id = fam_id
        self.ind_id = ind_id
        self.pat_id = pat_id
        self.mat_id = mat_id
        self.sex = sex
        self.phenotype = phenotype
        self.father = None
        self.mother = None
        self.offspring = []

    def __hash__(self):
        return hash(self.ind_id)

    def __repr__(self):
        return "PedigreeSample<%s:%s %s>" % (self.fam_id, self.ind_id, self.sex)

    def __str__(self):
        return "\t".join([self.fam_id, self.ind_id, self.pat_id,
                          self.mat_id, self.pat_id, self.sex,
                          "\t".join(self.phenotype)])


# Convience objects to start a multi processed consumer/producer strategy
class Consumer(multiprocessing.Process):

    """
    Basic Consumer. Follow the two queues with your *args and **kwargs that should be sent
    to the task when __call__ 'd

    NOTE! args can't hold anything that isn't pickle-able for the subprocess
    task_queue, result_queue

    Should add timeout functionality - I know I have it somewhere
    """

    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        try:
            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    # Poison pill means shutdown
                    self.task_queue.task_done()
                    break
                try:
                    next_task()
                except Exception as e:  # pylint: disable=broad-except
                    logging.error("Exception raised in task - %s", str(e))

                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    logging.error("Dumping Traceback:")
                    traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stderr)

                    next_task.failed = True
                    next_task.errMessage = str(e)

                self.result_queue.put(next_task)
                self.task_queue.task_done()

            return
        except Exception as e:  # pylint: disable=broad-except
            logging.error("Consumer %s Died\nERROR: %s", self.name, e)
            return


class ConsumerPool():
    """
    A resource for making a pool of consumer multiprocesses

    The tasks passed in via put must be callable (__call__)
    finished tasks are then yielded back.
    """

    def __init__(self, threads=1):
        """
        Does all the work
        """
        self.threads = threads
        self.input_queue = multiprocessing.JoinableQueue()
        self.output_queue = multiprocessing.Queue()
        self.processes = [Consumer(self.input_queue, self.output_queue)
                          for i in range(self.threads)]
        self.task_count = 0

    def start_pool(self):
        """
        run start on all processes
        """
        for proc in self.processes:
            proc.start()

    def put_task(self, task):
        """
        Add a callable task to the input_queue
        """
        self.task_count += 1
        self.input_queue.put(task)

    def put_poison(self):
        """
        For each process, add a poison pill so that it will close
        once the input_queue is depleted
        """
        for i in range(self.threads):
            logging.debug("Putting poison %d", i)
            self.input_queue.put(None)

    def get_tasks(self):
        """
        Yields the finished tasks
        """
        remaining = self.task_count
        while remaining:
            ret_task = self.output_queue.get()
            remaining -= 1
            yield ret_task

        if not self.input_queue.empty():
            raise RuntimeError("A worker thread quit unexpectedly, aborting.")

def package_versions():
    """Returns a list of tuples (packagename, version) of packages that are currently loaded
into this python instance.  Tries to guess which are system and internal packages, and skips
those."""
    version_table = {}
    for modname, mod in sys.modules.items():
        try:
            version = mod.__version__
            if not re.search(r"\d", str(version)):
                continue
        except AttributeError:
            # Not a versioned package
            continue
        try:
            path = mod.__path__
        except AttributeError:
            path = []
        try:
            path.append(mod.__file__)
        except AttributeError:
            pass
        try:
            package = mod.__package__
            if package and package != modname and not modname.startswith(package):
                # Not sure what the real name of this package is; include both
                # package name and module name.
                modname = f"{package}?{modname}"
        except AttributeError:
            pass
        # Skip system packages
        if any(p.startswith("/usr/lib/python") for p in path):
            continue
        # Skip internal packages
        if "._" in modname or modname[0] == "_":
            continue

        version_table[modname] = version

    # Skip modules whose versions are the same as their parent packages.
    versions = []
    for pkg in sorted(version_table.keys()):
        version = version_table[pkg]
        parts = pkg.rsplit(".", 1)
        if len(parts) > 1 and parts[0] in version_table:
            parent_version = version_table[parts[0]]
            if parent_version == version:
                continue

        versions.append((pkg, version))
    return versions

def timestamp(datestr=None):
    ''' Return a properly formatted timestamp. If datestr is supplied, convert it to datetime '''
    datefmt = '%Y-%m-%d %H:%M:%S.%f'
    if datestr:
        return datetime.strptime(datestr, datefmt)

    return datetime.now().strftime(datefmt)

def typed(dtype, value):
    ''' Return the value converted to the specified type. Raises if the type name is not valid. '''
    types = {
        'str': str,
        'int': int,
        'float': float,
        'timestamp': timestamp
    }
    if dtype in types:
        return types[dtype](value)

    raise SystemExit(f"Invalid type: {dtype}")

def plural(num):
    ''' Return 's' if num is anything other than 1 '''
    if num == 1:
        return ''
    return 's'

def confirm(prompt="y/N: ", default=False):
    '''
    If connected to a terminal, confirm with the user.
    Return the default if stdin is not a terminal, or if the user hits enter.
    Return False if ^C is detected.
    '''
    if not sys.stdin.isatty():
        return default

    try:
        reply = input(prompt).lower()
    except KeyboardInterrupt:
        return False

    if not reply or reply in ('y', 'yes'):
        return True

    return False

def chunked(chunks, chunk_size):
    '''
    Break a list of chunks into a list of chunk_size lists
    '''
    return [(chunks[i:i + chunk_size]) for i in range(0, len(chunks), chunk_size)]
