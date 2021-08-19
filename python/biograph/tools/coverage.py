"""
Genotype events in a VCF

VCF must have one entry per chrom:start-end. Use `bcftools norm -m+any`
By default, check the BioGraph file for the sample to annotate.
If BG sample name is not in the vcf's sample fields, `--sample` must be specified.

Outputs all VCF records and the single SAMPLE column  (i.e. other samples are removed)
"""
import os
import re
import sys
import gzip
import json
import signal
import argparse
import multiprocessing
import threading
import queue
import traceback
from collections import OrderedDict, defaultdict
from time import time
from setproctitle import setproctitle  # pylint: disable=no-name-in-module

import tabix
import vcf
import joblib
import pandas as pd

import biograph
import biograph.variants as bgexvar
import biograph.coverage as bganno
import biograph.tools.log as log

from biograph.tools.refhash import refhash

import numpy as np

# Set to 'True' to run work in separate processes.  Set to 'False' to
# run work in separate threads.
USE_MULTIPROCESSING = True

if USE_MULTIPROCESSING:
    THREAD_TYPE = multiprocessing.Process
    QUEUE_TYPE = multiprocessing.JoinableQueue
else:
    THREAD_TYPE = threading.Thread
    QUEUE_TYPE = queue.Queue

def build_table_header(): # pylint: disable=too-many-statements
    """
    Carefully build the column_name, dtypes OrderedDict of a parsed VCF
    """
    meta_header = OrderedDict()
    meta_header["key"] = str
    meta_header['sample'] = str
    meta_header['chrom'] = str
    meta_header['start'] = np.uint32
    meta_header['end'] = np.uint32
    meta_header['var_type'] = str

    info_header = OrderedDict()
    info_header["POP"] = bool
    # lies
    info_header["NUMASM"] = np.uint16
    info_header["SCORE"] = np.uint32
    info_header["VARLEN"] = np.uint32
    info_header["REFSPAN"] = np.uint16
    info_header["ASMLEN"] = np.uint16
    info_header["LANCH"] = np.uint8
    info_header["RANCH"] = np.uint8
    info_header["REFGC"] = np.float16
    info_header["ALTGC"] = np.float16

    fmt_header = OrderedDict()
    # Categorical
    # need to categorize these... ./. 0/0 0/1 1/1 # only possibilities
    fmt_header["GT"] = np.uint8
    # fmt_header["PG"] = np.int32 # don't use..
    fmt_header["GQ"] = np.uint8
    # fmt_header["PI"] = don't use
    fmt_header["OV"] = np.uint8
    fmt_header["DP"] = np.uint16
    #split where _r is ref-allele and _a is alt-allele
    fmt_header["AD_r"] = np.uint16
    fmt_header["AD_a"] = np.uint16
    fmt_header["PDP"] = np.uint16
    fmt_header["PAD_r"] = np.uint16
    fmt_header["PAD_a"] = np.uint16
    fmt_header["US_r"] = np.uint16
    fmt_header["US_a"] = np.uint16
    fmt_header["DS_r"] = np.uint16
    fmt_header["DS_a"] = np.uint16
    fmt_header["UC_r"] = np.uint16
    fmt_header["UC_a"] = np.uint16
    fmt_header["DC_r"] = np.uint16
    fmt_header["DC_a"] = np.uint16
    fmt_header["UDC_r"] = np.uint16
    fmt_header["UDC_a"] = np.uint16
    fmt_header["UCC_r"] = np.uint16
    fmt_header["UCC_a"] = np.uint16
    fmt_header["DDC_r"] = np.uint16
    fmt_header["DDC_a"] = np.uint16
    fmt_header["DCC_r"] = np.uint16
    fmt_header["DCC_a"] = np.uint16
    fmt_header["UMO_r"] = np.uint16
    fmt_header["UMO_a"] = np.uint16
    fmt_header["DMO_r"] = np.uint16
    fmt_header["DMO_a"] = np.uint16
    fmt_header["UXO_r"] = np.uint16
    fmt_header["UXO_a"] = np.uint16
    fmt_header["DXO_r"] = np.uint16
    fmt_header["DXO_a"] = np.uint16
    fmt_header["NR_r"] = np.uint16
    fmt_header["NR_a"] = np.uint16
    fmt_header["MO_r"] = np.uint16
    fmt_header["MO_a"] = np.uint16
    fmt_header["XO_r"] = np.uint16
    fmt_header["XO_a"] = np.uint16
    fmt_header["XC_r"] = np.uint16
    fmt_header["XC_a"] = np.uint16
    fmt_header["AC_r"] = np.uint16
    fmt_header["AC_a"] = np.uint16
    fmt_header["MC_r"] = np.uint16
    fmt_header["MC_a"] = np.uint16
    fmt_header["EC_r"] = np.uint16
    fmt_header["EC_a"] = np.uint16
    fmt_header["PL_ref"] = np.uint8
    fmt_header["PL_het"] = np.uint8
    fmt_header["PL_hom"] = np.uint8
    fmt_header["RC"] = np.uint16
    fmt_header["AC_LR"] = np.uint32
    fmt_header["AC_AB"] = np.uint32
    fmt_header["AC_TA"] = np.uint32

    ret_header = OrderedDict()
    ret_header.update(meta_header)
    ret_header.update(info_header)
    ret_header.update(fmt_header)
    return ret_header

def get_type_lens(entry, mREF=None, mALT=None):
    """
    Parse an entry and return its sv_type and its sv_len
    """
    if entry is not None:
        mREF = entry.ref
        mALT = entry.alts[0]
    sv_types = None
    sv_lens = None
    # Get type for counting - MYTYPES
    if len(mREF) == len(mALT):
        sv_types = "REPL"
        sv_lens = len(mREF)
    elif len(mREF) == 1:
        sv_types = "INS"
        sv_lens = len(mALT) - 1
    elif len(mALT) == 1:
        sv_types = "DEL"
        sv_lens = len(mREF) - 1
    elif len(mREF) > len(mALT):
        sv_types = "SUBSDEL"
        sv_lens = len(mREF) - len(mALT)
    elif len(mALT) > len(mREF):
        sv_types = "SUBSINS"
        sv_lens = len(mALT) - len(mREF)
    else:
        log.error(str(entry))
        log.error("shouldn't have some new crazy type\n")
        exit(1)

    # Cheating and assuming it's single variant per-line vcf
    return sv_types, sv_lens

def info_to_dict(infos):
    """
    Convert INFO field to dict
    """
    ret = {}
    for info in infos:
        s = info.split("=", 1)
        if len(s) == 1:
            ret[info] = None
        else:
            ret[s[0]] = s[1]
    return ret

def parse_info(cur_data, v):
    """
    Parse the ML features for the longest discovered assembly from info into the cur_data
    """
    if "POP" in v:
        cur_data["POP"] = v["POP"]
    else:
        cur_data["POP"] = False

def parse_format(cur_data, var_sample): #pylint: disable=too-many-statements
    """
    Parse the format information from a sample into cur_data
    """
    # Need to see what all these could be...

    if var_sample["GT"] == (0, 0) or var_sample["GT"] == "0/0":
        cur_data["GT"] = 0
    elif var_sample["GT"] == (0, 1) or var_sample["GT"] == "0/1":
        cur_data["GT"] = 1
    elif var_sample["GT"] == (1, 1) or var_sample["GT"] == "1/1":
        cur_data["GT"] = 2
    else:
        cur_data["GT"] = 3
    # BioGraph Specific
    # need to get these if they exist
    cur_data["NUMASM"] = var_sample["NUMASM"] if "NUMASM" in var_sample else None
    cur_data["SCORE"] = var_sample["LASCORE"] if "LASCORE" in var_sample else None
    cur_data["REFSPAN"] = var_sample["LAREFSPAN"] if "LAREFSPAN" in var_sample else None
    cur_data["ASMLEN"] = var_sample["LAALTSEQLEN"] if "LAALTSEQLEN" in var_sample else None
    cur_data["LANCH"] = var_sample["LALANCH"] if "LALANCH" in var_sample else None
    cur_data["RANCH"] = var_sample["LARANCH"] if "LARANCH" in var_sample else None
    cur_data["REFGC"] = var_sample["LAREFGC"] if "LAREFGC" in var_sample else None
    cur_data["ALTGC"] = var_sample["LAALTGC"] if "LAALTGC" in var_sample else None
    cur_data["OV"] = var_sample["OV"] if "OV" in var_sample else None
    cur_data["PDP"] = var_sample["PDP"] if "PDP" in var_sample else None
    cur_data["PAD_r"] = var_sample["PAD"][0] if "PAD" in var_sample else None
    cur_data["PAD_a"] = var_sample["PAD"][1] if "PAD" in var_sample else None

    # pcmp specific
    cur_data["GQ"] = var_sample["GQ"] if var_sample["GQ"] is not None else 0
    cur_data["DP"] = var_sample["DP"] # be careful these aren't '.'
    #split where _r is ref-allele and _a is alt-allele
    cur_data["AD_r"] = var_sample["AD"][0]
    cur_data["AD_a"] = var_sample["AD"][1]
    cur_data["US_r"] = var_sample["US"][0]
    cur_data["US_a"] = var_sample["US"][1]
    cur_data["DS_r"] = var_sample["DS"][0]
    cur_data["DS_a"] = var_sample["DS"][1]
    cur_data["UC_r"] = var_sample["UC"][0]
    cur_data["UC_a"] = var_sample["UC"][1]
    cur_data["DC_r"] = var_sample["DC"][0]
    cur_data["DC_a"] = var_sample["DC"][1]
    cur_data["UDC_r"] = var_sample["UDC"][0]
    cur_data["UDC_a"] = var_sample["UDC"][1]
    cur_data["UCC_r"] = var_sample["UCC"][0]
    cur_data["UCC_a"] = var_sample["UCC"][1]
    cur_data["DDC_r"] = var_sample["DDC"][0]
    cur_data["DDC_a"] = var_sample["DDC"][1]
    cur_data["DCC_r"] = var_sample["DCC"][0]
    cur_data["DCC_a"] = var_sample["DCC"][1]
    cur_data["UMO_r"] = var_sample["UMO"][0]
    cur_data["UMO_a"] = var_sample["UMO"][1]
    cur_data["DMO_r"] = var_sample["DMO"][0]
    cur_data["DMO_a"] = var_sample["DMO"][1]
    cur_data["UXO_r"] = var_sample["UXO"][0]
    cur_data["UXO_a"] = var_sample["UXO"][1]
    cur_data["DXO_r"] = var_sample["DXO"][0]
    cur_data["DXO_a"] = var_sample["DXO"][1]
    cur_data["NR_r"] = var_sample["NR"][0]
    cur_data["NR_a"] = var_sample["NR"][1]
    cur_data["MO_r"] = var_sample["MO"][0]
    cur_data["MO_a"] = var_sample["MO"][1]
    cur_data["XO_r"] = var_sample["XO"][0]
    cur_data["XO_a"] = var_sample["XO"][1]
    cur_data["XC_r"] = var_sample["XC"][0]
    cur_data["XC_a"] = var_sample["XC"][1]
    cur_data["AC_r"] = var_sample["AC"][0]
    cur_data["AC_a"] = var_sample["AC"][1]
    cur_data["MC_r"] = var_sample["MC"][0]
    cur_data["MC_a"] = var_sample["MC"][1]
    cur_data["EC_r"] = var_sample["EC"][0]
    cur_data["EC_a"] = var_sample["EC"][1]
    if "PL" in var_sample and var_sample["PL"] != '.':
        cur_data["PL_ref"] = var_sample["PL"][0] if var_sample["PL"][0] is not None else 0
        cur_data["PL_het"] = var_sample["PL"][1] if var_sample["PL"][0] is not None else 0
        cur_data["PL_hom"] = var_sample["PL"][2] if var_sample["PL"][0] is not None else 0
    else:
        cur_data["PL_ref"] = 0
        cur_data["PL_het"] = 0
        cur_data["PL_hom"] = 0

    cur_data["RC"] = var_sample["RC"]

def df_row_from_str(entry, samp_data):
    """
    Convert a vcf entry to a dataframe row
    """
    cur_data = defaultdict()
    # This is a total hack to stop errors
    for key in samp_data:
        if isinstance(samp_data[key], str) and ',' in samp_data[key]:
            samp_data[key] = samp_data[key].split(',')
        elif samp_data["GQ"] == '.':
            samp_data["GQ"] = None
        elif samp_data["RC"] == ".":
            samp_data["RC"] = None
    start = int(entry[1]) - 1
    stop = start + len(entry[3])
    cur_data["key"] = "%s:%d-%d.%s" % (entry[0], start, stop, entry[4])
    cur_data["chrom"] = entry[0]
    cur_data["start"] = start
    cur_data["end"] = stop
    var_type, var_len = get_type_lens(None, entry[3], entry[4])
    cur_data["var_type"] = var_type
    cur_data["VARLEN"] = var_len
    var_info = info_to_dict(entry[7])
    try:
        parse_info(cur_data, var_info)
        parse_format(cur_data, samp_data)
    except KeyError as e:
        raise RuntimeError('Could not find required VCF fields. Input VCF must be run through pcmp first. %s' % (str(e)))
    return cur_data

class PCMPWriter(THREAD_TYPE):  # pylint: disable=too-many-instance-attributes

    """
    Thread dedicated to writing the output
    """
    # Format fields BioGraph Coverage populates
    __cov_formats = ["GT", "GQ", "PL", "DP", "AD", "AC_LR", "AC_AB", "AC_TA"]

    def __init__(self, vcf_output_name, df_output_name, template_vcf, sample, annotations, strip_fmt=False):
        """
        Using a template_vcf, edit the header appropriately
        and create the output_name
        annotations are objects with
        """
        super().__init__()
        self.queue = QUEUE_TYPE()
        self.template_vcf = template_vcf
        self.output_name = vcf_output_name
        self.df_name = df_output_name
        self.out_vcf = None
        self.sample = sample
        self.sample_column = None
        self.annotations = annotations
        self.new_headers = []
        self.new_tags = []
        # could put safety checks here to make sure we don't replicate IDs..?
        self.format_defaults = {}

        # These are the fields coverage populates in the VCF
        self.vcf_formats = list(self.__cov_formats)
        # Optionally preserve the remaining existing formats
        if not strip_fmt:
            m_vcf = vcf.VCFReader(filename=self.template_vcf)
            for fmt in m_vcf.formats.keys():
                if fmt not in self.vcf_formats:
                    self.vcf_formats.append(fmt)
        self.vcf_format_str = ":".join(self.vcf_formats)

        for anno in self.annotations:
            # Hack to prevent CovAnno from being placed in VCF header
            if not isinstance(anno, bganno.CovAnno):
                self.new_headers.extend(anno.get_header())
            self.new_tags.extend(anno.get_format_tags())
            self.format_defaults.update(anno.get_blank_format())
        self.header = self.make_header()

    def make_header(self):
        """
        Puts in the new annotation header lines
        """
        header = ""
        tag_extractor = re.compile("##FORMAT=<ID=([^,]+)")
        with gzip.GzipFile(self.template_vcf) as fh:
            for line in fh:
                line = line.decode()
                if line.startswith("##"):
                    match = tag_extractor.match(line)
                    if not (match and match.groups()[0] in self.new_tags):
                        header += line
                elif line.startswith("#"):
                    # add header stuff
                    h_dat = line.strip().split('\t')
                    try:
                        self.sample_column = h_dat.index(self.sample)
                        h_dat[9] = h_dat[self.sample_column]
                    except ValueError:
                        h_dat[9] = self.sample
                    header += "\n".join(self.new_headers) + '\n'
                    header += ('##source="Spiral Genetics BioGraph Genotype",version="%s",' %
                               biograph.version())
                    header += ('description="build-revision=\'%s\',command-line=\'%s\'"\n' %
                               (biograph.build_revision(), " ".join(sys.argv)))
                    if '##refhash=' not in header:
                        rh = refhash(self.template_vcf)
                        header += f"##{rh.verbose_name()}\n"
                    header += "\t".join(h_dat[:10]) + '\n'
                    break
        return header

    def format_entry(self, vcf_entry, sample_field, new_fmt):
        """
        Edit vcf_entry with fmt_dict {FORMAT:Value} before writing
        """
        fmt = vcf_entry[8].split(':')
        if not sample_field or sample_field == '.':
            sample_field = ":".join(["."] * len(fmt))

        samp_data = {x:y for x, y in zip(fmt, sample_field.split(':'))}

        for i in self.vcf_formats:
            if i not in samp_data:
                samp_data[i] = "."

        for tag in self.format_defaults:
            samp_data[tag] = self.format_defaults[tag]

        samp_data.update(new_fmt)

        # I need to subset this to only the fields I want to output
        vcf_entry[8] = self.vcf_format_str
        sample_field = []
        for key in self.vcf_formats:
            if key not in samp_data:
                continue
            val = samp_data[key]
            if isinstance(val, float):
                val = '{:.2f}'.format(val)
            sample_field.append(str(val))

        df_data = df_row_from_str(vcf_entry, samp_data)

        return ('\t'.join(vcf_entry[:9]) + '\t' + ":".join(sample_field)+ '\n').encode(), df_data

    def run(self):
        """
        Continually get from the queue and write the record
        """
        # Open in binary mode so we don't have to worry about having a
        # bottleneck here dealing with unicode.
        self.out_vcf = open(self.output_name, 'wb')
        self.out_vcf.write(self.header.encode())
        df_data = []
        while True:
            # queue holds tuples of vcf_data str and df_data list
            data = self.queue.get()
            if data is None:
                break
            elif isinstance(data, Exception):
                raise data
            else:
                self.out_vcf.write(data[0])
                df_data.extend(data[1])
        log.info("Making DataFrame")
        data = pd.DataFrame(df_data)
        data["sample"] = "sample"
        dtypes = build_table_header()
        for col in data.columns:
            try:
                data[col] = data[col].astype(dtypes[col])
            except TypeError as e:
                log.debug(f"Unable to cast {col} ({str(e)})")

        #data = data.astype(build_table_header())
        joblib.dump(data, self.df_name)
        self.out_vcf.close()


class PCMP(THREAD_TYPE):  # pylint: disable=too-many-instance-attributes

    """
    Fetches vcf entries in a region and calculates coverage to enotype
    """

    # pylint: disable=too-many-arguments
    def __init__(self, bg_file, rmap, ref_file, vcf_file, inputq, outputter, annotations,
                 ref_pad_bases, passonly=None, phasing=False):
        """
        Runs over regions
        """
        super().__init__()
        self.bg_file = bg_file
        self.rmap = rmap
        self.ref_file = ref_file
        self.vcf_file = vcf_file
        self.inputq = inputq
        self.outputter = outputter
        self.annotations = annotations
        if outputter.sample_column:
            self.sample_index = outputter.sample_column - 9
        else:
            self.sample_index = None
        self.passonly = passonly
        self.ref_pad_bases = ref_pad_bases
        self.error = RuntimeError("Unknown error")
        self.formatted = []
        self.phasing = phasing

    def process_region(self, chrom, start, end):
        """
        work on all the entries in a region
        # this needs to fetch the vcf entries
        # Feed the entries to annotators
        # annotators return the vcf_entry edited with new Annotations.get_fmt() dict
        #               can eventually expand to get_info..?
        # and here we forward the information along to write(vcf_entry)
        """
        # 3-way flag. None is all
        passonly = None if not self.passonly else True
        entries = bganno.VcfEntryInfo.parse_region(self.vcf_file, chrom, start, end,
                                                   sample_index=self.sample_index,
                                                   passonly=passonly, phasing=self.phasing)
        entries = bgexvar.trim_ref(self.ref_file, chrom, entries)
        entries = bgexvar.add_ref_assemblies(self.ref_file, chrom, entries, self.ref_pad_bases)

        # Chain them
        for anno in self.annotations:
            entries = anno.parse(self.rmap, entries)

        # Output them
        count = 0
        for count, variant in enumerate(entries):
            self.output_entry(variant)

        # Still need to output the non-PASS
        if self.passonly:
            for variant in bganno.VcfEntryInfo.parse_region(self.vcf_file, chrom, start, end, False):
                self.output_entry(variant)
        self.flush_formatted()

        return count + 1

    def output_entry(self, variant):
        """
        Buffers the given variant for output to the vcf writer
        """
        formatted = self.outputter.format_entry(list(variant.vcf_entry_info.vcf_entry),
                                                variant.vcf_entry_info.sample_field,
                                                dict(variant.vcf_entry_info.new_fmt))
        self.formatted.append(formatted)
        if len(self.formatted) > 100:
            self.flush_formatted()

    def flush_formatted(self):
        """
        Sends the pending buffered vcf entries to the vcf writer
        """
        if not self.formatted:
            return
        vcfdat = b"".join([_[0] for _ in self.formatted])
        dfdat = [_[1] for _ in self.formatted]
        self.outputter.queue.put((vcfdat, dfdat))
        self.formatted = []

    def abort(self, signum=None, frame=None):  # pylint: disable=unused-argument
        """
        Log an error, put the exception on the queue, and exit 1.
        Can be called directly or as a signal handler.
        """
        if signum:
            log.info(f"Aborting on signal {signum}: {self.error}")
        else:
            log.info(f"Aborting: {self.error}")
        self.outputter.queue.put_nowait(self.error)
        sys.exit(1)

    def run(self):
        """
        Continually polls the queue for something else to work on
        """
        is_multiprocess = isinstance(self, multiprocessing.Process)
        if is_multiprocess:
            # Abort on INT or TERM
            self.error = RuntimeError('Worker terminated by signal')
            signal.signal(signal.SIGINT, self.abort)
            signal.signal(signal.SIGTERM, self.abort)

        try:
            while True:
                item = self.inputq.get()
                if item is None:
                    break
                chrom, start, end = item
                msg = f"{chrom}:{start}-{end}"
                start_time = time()
                log.debug(f"Starting region {msg}, estimated input queue size = {self.inputq.qsize()}")
                if is_multiprocess:
                    setproctitle(f"biograph coverage worker: {msg}")
                count = self.process_region(chrom, start, end)
                if is_multiprocess:
                    setproctitle("biograph coverage worker: idle")
                self.inputq.task_done()
                delta = float(time() - start_time)
                log.debug(f"{msg} done in {delta:0.2f}s ({count} v, {count / delta:0.0f} v/s, {(end - start) / delta:0.0f} b/s)")

        except Exception as e:  # pylint: disable=broad-except
            log.error(f"When processing {msg}: {traceback.format_exc()}")
            self.error = e
            self.abort()

class MICROCONTIGS(THREAD_TYPE):
    """ Break contigs into microcontigs """
    def __init__(self, vcf_file, contigq, workerq, min_clearance, min_contig_size):
        """
        Runs over regions
        """
        super().__init__()
        self.vcf_file = vcf_file
        self.contigq = contigq
        self.workerq = workerq
        self.min_clearance = min_clearance
        self.min_contig_size = min_contig_size
        self.error = RuntimeError("Unknown error")
        self.total_mcs = 0

    @staticmethod
    def find(region, vcf_file, min_clearance, min_contig_size, workerq):
        """
        Break regions into smaller microcontigs separated by a minimum clearance.
        This function is static to allow easier testing.
        """
        (chrom, chrom_start, chrom_end) = region
        chrom_start = chrom_start - 1
        last_start = 0
        last_end = 0
        padding = min_clearance // 2
        total_mcs = 0
        t = tabix.open(vcf_file)
        try:
            for v in t.query(chrom, chrom_start, chrom_end):
                pos = int(v[1]) - 1
                #     pos + len(ref) - 1 (since VCF always has at least 1 ref base)
                end = pos + len(v[3]) - 1

                if last_start == 0:
                    last_start = pos
                    last_end = end
                    continue

                # tabix returns any variant that overlaps a region. Skip
                # variants that don't start within a region
                if not (chrom_start <= pos < chrom_end):
                    continue

                # Is there room to break this contig?
                # And the microcontig is above the minimum size
                if pos - last_end > min_clearance and last_end - last_start + min_clearance > min_contig_size:
                    workerq.put((chrom, max(0, last_start - padding), min(last_end + padding, chrom_end)))
                    last_start = pos
                    last_end = end
                else: # otherwise, continue expanding
                    last_end = max(end, last_end)
                total_mcs += 1
            if last_start > 0:
                workerq.put((chrom, max(0, last_start - padding), min(end + padding, chrom_end)))

        # No variants in this region? No problem.
        except tabix.TabixError:
            pass
        return total_mcs

    def run(self):
        """
        Continually polls the queue for something else to work on
        """
        is_multiprocess = isinstance(self, multiprocessing.Process)
        while True:
            item = self.contigq.get()
            if item is None:
                break

            if is_multiprocess:
                setproctitle(f"biograph coverage worker: find microcontigs {item[0]}:{item[1]}-{item[2]}")

            self.total_mcs = self.find(item, self.vcf_file, self.min_clearance, self.min_contig_size, self.workerq)
            self.contigq.task_done()


def parse_args(clargs):
    """ Command line UI """
    parser = argparse.ArgumentParser(prog="coverage", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--biograph", metavar="BG", required=True,
                        help="Merged BioGraph file containing individuals")
    parser.add_argument("-v", "--variants", metavar="VCF", required=True,
                        help="The vcf containing all samples")
    parser.add_argument("-r", "--reference", metavar="REF", required=True,
                        help="Reference genome folder")
    parser.add_argument("-R", "--regions", "--bed", default=None,
                        help="Bed file of regions to process")
    parser.add_argument("-s", "--sample", default=None,
                        help="Sample in merged BioGraph to use (%(default)s)")
    parser.add_argument("-o", "--output", metavar="OUT", default="/dev/stdout",
                        help="Annotated vcf file output (%(default)s)")
    parser.add_argument("-d", "--dataframe", metavar="DF", default="coverage.df",
                        help="Output coverage DataFrame (%(default)s)")
    parser.add_argument("--strip-fmt", action="store_true",
                        help="Strip any existing FORMAT fields")
    parser.add_argument("-m", "--min-insert", default=200, type=int,
                        help="Minimum insert size to consider paired (%(default)s)")
    parser.add_argument("-M", "--max-insert", default=1000, type=int,
                        help="Maximum insert size to consider paired (%(default)s)")
    parser.add_argument("--min-clearance", default=2000, type=int,
                        help="Find microcontig boundaries separated by at least this many reference bases (0 to disable, default = %(default)s)")
    parser.add_argument("--min-contig-size", default=3000, type=int,
                        help="Minimum contig size to make when computing microcontigs (default = %(default)s)")
    parser.add_argument("-p", "--passonly", action="store_true",
                        help="Only attempt to genotype PASS variants")
    parser.add_argument("--cache", default=False, action="store_true",
                        help="Attempt to cache as much as possible in RAM")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("--max-reads-per-entry", type=int, default=0,
                        help="If nonzero, maximum number of reads per seqset entry to use; additional reads are ignored (%(default)s)")
    parser.add_argument("--max-coverage-paths", type=int, default=3000,
                        help="If nonzero, limit number of parallel coverage paths to trace (%(default)s)")
    parser.add_argument("--phasing", action="store_true", help="Skip paths counterindicated by phasing information")
    parser.add_argument("--debug", action="store_true", help="Verbose logging")
    parser.add_argument("--max-coverage-alleles", type=int, default=50,
                        help="Maximum number of alleles through reference to calculate coverage on, including the reference path.  Only used when phasing is enabled.")
    parser.add_argument("--filter-dup-align", action="store_true", help="Count each read only once per local area")
    parser.add_argument("--ideal-insert", type=int, default=None, help="Place paired reads only once, optimizing for closest to the given insert size.")
    parser.add_argument("--placer-max-ambig", type=int, default=15, help="Maximum number of ambiguous choices before the pair placer will discard a read")

    args = parser.parse_args(clargs)
    log.setup_logging(args.debug)
    log.debug(f"Params:\n{json.dumps(vars(args), indent=4)}")

    if not args.variants.endswith(".vcf.gz") or not os.path.exists(args.variants + ".tbi"):
        log.error("Variants file must be compressed and indexed")
        exit(1)
    return args

def get_regions(regions, ref):
    """
    parse the regions file or use ref scaffolds
    """
    ret = []
    if regions is not None:
        with open(regions, 'r') as fh:
            for line in fh:
                data = line.strip().split('\t')
                ret.append((data[0], int(data[1]), int(data[2])))
    else:
        for ctg in ref.scaffolds:
            ret.append((ctg, 0, int(ref.scaffold_lens[ctg])))

    return ret

def main(clargs): # pylint: disable=too-many-locals, too-many-statements
    ''' Calculate coverage for VCF entries '''
    # NOTE: ^^^ this ^^^ docstring is used for the command description in __main__.py

    args = parse_args(clargs)

    # Ensure the vcf is indexed
    if not os.path.exists(args.variants + ".tbi"):
        log.error("Input vcf file needs to be indexed (.tbi)")
        exit(1)

    if args.min_clearance > 0 and args.min_clearance < (2 * args.max_insert):
        log.warning(f"Increasing min-clearance to 2x max-insert ({2 * args.max_insert})")
        args.min_clearance = 2 * args.max_insert

    setproctitle("biograph coverage main")

    m_bg = biograph.BioGraph(args.biograph,
                             biograph.CacheStrategy.RAM if args.cache
                             else biograph.CacheStrategy.MMAPCACHE)
    m_ref = biograph.Reference(args.reference)

    # Figure out the sample
    if args.sample is None:
        if len(m_bg.metadata.samples) == 1:
            args.sample = next(iter(m_bg.metadata.samples))
    elif args.sample not in m_bg.metadata.samples:
        raise KeyError("Sample %s not present in BioGraph" % args.sample)

    t_vcf = vcf.Reader(filename=args.variants)
    if args.sample not in t_vcf.samples:
        # raise KeyError("Sample %s not present in VCF" % args.sample)
        log.warning(f"Sample {args.sample} not present in VCF")

    rmap = m_bg.open_readmap(args.sample)

    # Set your annotations here...
    annotations = []
    if args.phasing:
        annotations.append(bganno.PhaseConflictResolver())
    annotations.append(bganno.CovAnno(
        args.min_insert, args.max_insert,
        max_reads_per_entry=args.max_reads_per_entry,
        max_coverage_paths=args.max_coverage_paths,
        phasing=args.phasing,
        max_coverage_alleles=args.max_coverage_alleles,
        filter_dup_align=args.filter_dup_align,
        ideal_insert=args.ideal_insert,
        placer_max_ambig=args.placer_max_ambig
    ))
    annotations.append(bganno.GTAnno())
    annotations.append(bganno.ACAnno())
    out_writer = PCMPWriter(args.output, args.dataframe, args.variants, args.sample, annotations, args.strip_fmt)
    out_writer.start()

    mc_thread_count = max(1, args.threads // 2)
    pcmp_thread_count = max(1, args.threads - mc_thread_count)

    pcmp_threads = []
    mc_threads = []

    contigq = QUEUE_TYPE()
    workerq = QUEUE_TYPE()

    # Start workers
    log.info(f"Processing regions")

    # If building microcontigs, put regions on the contigq
    if args.min_clearance > 0:
        for region in get_regions(args.regions, m_ref):
            contigq.put(region)
    else:
        mc_thread_count = 0
        for region in get_regions(args.regions, m_ref):
            workerq.put(region)

    for _ in range(mc_thread_count):
        contigq.put(None)

    # spin up microcontig workers
    for _ in range(mc_thread_count):
        new_mc = MICROCONTIGS(
            vcf_file=args.variants,
            contigq=contigq,
            workerq=workerq,
            min_clearance=args.min_clearance,
            min_contig_size=args.min_contig_size
        )
        new_mc.start()
        mc_threads.append(new_mc)

    # spin up PCMP workers
    ref_pad_bases = args.max_insert + m_bg.seqset.max_sequence_length()

    for _ in range(pcmp_thread_count):
        new_pcmp = PCMP(
            bg_file=m_bg,
            rmap=rmap,
            ref_file=m_ref,
            vcf_file=args.variants,
            inputq=workerq,
            outputter=out_writer,
            annotations=annotations,
            ref_pad_bases=ref_pad_bases,
            passonly=args.passonly,
            phasing=args.phasing
        )

        new_pcmp.start()
        pcmp_threads.append(new_pcmp)
    total_mcs = 0
    for cur_thread in mc_threads:
        cur_thread.join()
        total_mcs += cur_thread.total_mcs
        if len(pcmp_threads) < args.threads:
            # mc graduates to pcmp
            new_pcmp = PCMP(
                bg_file=m_bg,
                rmap=rmap,
                ref_file=m_ref,
                vcf_file=args.variants,
                inputq=workerq,
                outputter=out_writer,
                annotations=annotations,
                ref_pad_bases=ref_pad_bases,
                passonly=args.passonly,
                phasing=args.phasing
            )
            new_pcmp.start()
            pcmp_threads.append(new_pcmp)
    if total_mcs:
        log.info(f"{total_mcs} regions to parse")
    for _ in pcmp_threads:
        workerq.put(None)

    cur_qsize = workerq.qsize()
    # Logging progress
    log.info(f"Awaiting ~{cur_qsize} regions")
    tot_done = 0
    INTERVAL = 300
    start_time = time()
    while pcmp_threads:
        # Aim to process through the thread pool for INTERVAL total time, spending per_thread_timeout on each thread.
        per_thread_timeout = INTERVAL / len(pcmp_threads)
        to_remove = []
        for cur_thread in pcmp_threads:
            cur_thread.join(timeout=per_thread_timeout)
            if not cur_thread.is_alive():
                if USE_MULTIPROCESSING and cur_thread.exitcode != 0:
                    e = RuntimeError(f"Worker {cur_thread} executed with exit code {cur_thread.exitcode}")
                    out_writer.queue.put(e)
                    raise e
                to_remove.append(cur_thread)

        for i in to_remove:
            pcmp_threads.remove(i)

        new_qsize = workerq.qsize()
        if new_qsize == 0:
            log.info(f"Awaiting final regions; {len(pcmp_threads)} workers remaining")
            INTERVAL = 60
            continue
        tot_time = time() - start_time
        tot_done += cur_qsize - new_qsize
        cur_qsize = new_qsize
        if not tot_done or tot_time < 1:
            continue
        rate = tot_done / tot_time
        remain = cur_qsize / rate
        m, s = divmod(remain, 60)
        h, m = divmod(m, 60)
        remain = '{:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s))
        log.info(f"Estimated {cur_qsize} regions ({remain}) remaining")

    # Finish up the writer
    out_writer.queue.put(None)
    out_writer.join()

    # There should be no more work to do
    if not workerq.empty():
        while True:
            data = workerq.get()
            if data is None:
                break
            print(data)
        raise RuntimeError("Worker died unexpectedly, aborting.")

    if out_writer.exitcode:
        raise RuntimeError("PCMPWriter died unexpectedly, aborting.")

    log.info("Finished")

if __name__ == '__main__':
    main(sys.argv[1:])
