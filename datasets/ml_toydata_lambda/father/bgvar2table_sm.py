"""
Converts a BioGraph discovery assembly.csv and variants.vcf into a joblib datastructure containing
 {"data": flat pandas DataFrame,
  "time": time it was build,
  "user": user who build this,
  "computer": computer that built this
  "args": command line arguments}
"""
# TODO: Make a util to auto build the time/user/computer/args metadata
import io
import os
import sys
import glob
import gzip
import joblib
import logging
import argparse

from collections import OrderedDict

import numpy as np
import pandas as pd
import pysam

from acebinf import setup_logging

def parse_args(args):
    """
    Parse the args
    """
    parser = argparse.ArgumentParser(prog="meanno", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--asm", type=str, required=True,
                        help="Assembly.csv file to parse")
    parser.add_argument("--mode", default="full", choices=["bench", "full"],
                        help="Write a 'full' table or a 'bench' table (%(default)s)")
    parser.add_argument("-v", "--vcf", type=str,
                        help="VCF file to parse")
    parser.add_argument("-s", "--sample", type=str, default="unknown",
                        help="Sample identifier (%(default)s)")
    parser.add_argument("-o", "--out", type=str, required=True,
                        help="Output joblib file to write")
    parser.add_argument("--trudir", default=None,
                        help="Truvari directory")
    parser.add_argument("--rtgdir", default=None,
                        help="RTG directory")

    args = parser.parse_args(args)
    setup_logging()
    if args.mode == 'bench' and (args.trudir is None or args.rtgdir is None):
        logging.error("Bench mode requires --trudir and --rtgdir")
        exit(1)
    elif args.mode == 'plain' and args.var is None:
        logging.error("Plain mode requires --vcf")
        exit(1)
    return args

def load_asm(asm_fn):
    """
    Read assembly.csv
    return {"asm_id": ["SCORE" # score of the assembly
                     , "REFSPAN" # length of reference-sequence
                     , "ASMLEN" # Length of the assembly-sequence
                     , "LANCH" # anchor length
                     , "RANCH" # anchor length
                     , "REFGC" # gcpct of ref sequence
                     , "ALTGC" # gcpct of alt sequence
                     , "ALTSEQ" # kept for optional mapping later ]}
    """
    logging.info("Loading assemblies")
    ret_header = OrderedDict()
    ret_header["asm_id"] = str
    ret_header["SCORE"] = np.int32
    ret_header["REFSPAN"] = np.int16
    ret_header["ASMLEN"] = np.int16
    ret_header["LANCH"] = np.int8 # probably could be 8s
    ret_header["RANCH"] = np.int8 #
    ret_header["REFGC"] = np.float16
    ret_header["ALTGC"] = np.float16
    ret_header["ALTSEQ"] = str
    if asm_fn.endswith(".gz"):
        fh = io.TextIOWrapper(gzip.open(asm_fn))
    else:
        fh = open(asm_fn)
    
    header = fh.readline().strip().split(',')
    iaid = header.index("aid")
    iscore = header.index("score")
    iref_seq = header.index("ref_seq")
    ialt_seq = header.index("seq")
    ileft_anchor = header.index("left_anchor_len")
    iright_anchor = header.index("right_anchor_len")
    ret = {}
    for line in fh:
        data = line.strip().split(',')
        ref_len = len(data[iref_seq])
        alt_len = len(data[ialt_seq])
        if ref_len != 0:
            ref_gc = (data[iref_seq].count("G") + data[iref_seq].count("C")) / ref_len
        else:
            ref_gc = 0
        if alt_len != 0:
            alt_gc = (data[ialt_seq].count("G") + data[ialt_seq].count("C")) / alt_len
        else:
            alt_gc = 0
        ret[int(data[iaid])] = [int(data[iaid]),
                                int(data[iscore]),
                                ref_len,
                                alt_len,
                                int(data[ileft_anchor]),
                                int(data[iright_anchor]),
                                ref_gc,
                                alt_gc,
                                data[ialt_seq]]
    logging.info("%d assemblies loaded", len(ret))
    return ret_header, ret
 
def get_type_lens(entry):
    """
    Parse an entry and return it's sv_type and it's sv_len
    """
    mREF = entry.ref
    # TODO - should get the longest?
    mALTs = entry.alts[:1]
    sv_types = []
    sv_lens = []
    # Get type for counting - MYTYPES
    for mALT in mALTs:
        if len(mREF) == len(mALT):
            sv_types.append("REPL")
            sv_lens.append(len(mREF))
        elif len(mREF) == 1:
            sv_types.append("INS")
            sv_lens.append(len(mALT) - 1)
        elif len(mALT) == 1:
            sv_types.append("DEL")
            sv_lens.append(len(mREF) - 1)
        elif len(mREF) > len(mALT):
            sv_types.append("SUBSDEL")
            sv_lens.append(len(mREF) - len(mALT))
        elif len(mALT) > len(mREF):
            sv_types.append("SUBSINS")
            sv_lens.append(len(mALT) - len(mREF))
        else:
            logging.error(str(entry))
            logging.error("shouldn't have some new crazy type\n")
            exit()

    # MYSIZES
    ret_lens = []
    for sv_len in sv_lens:
        if sv_len < 10:
            ret_lens.append("1-9")
        elif sv_len < 50:
            ret_lens.append("10-49")
        elif sv_len < 100:
            ret_lens.append("50-99")
        elif sv_len < 300:
            ret_lens.append("100-299")
        elif sv_len < 1000:
            ret_lens.append("300-999")
        else:
            ret_lens.append("gt1000")
    # Cheating and assuming it's single variant per-line vcf
    return sv_types[0], sv_lens[0]


def parse_format(var_sample):
    """
    Parse the format information from a sample
    """
    # ugh
    ret = []
    # Parsing format information
    # Need to see what all these could be...
    if None in var_sample["GT"]:
        ret.append(3)
    elif var_sample["GT"] == (0, 0):
        ret.append(0)
    elif var_sample["GT"] ==  (0, 1):
        ret.append(1)
    elif var_sample["GT"] == (1, 1):
        ret.append(2)
        
    ret.extend([var_sample["GQ"] if var_sample["GQ"] is not None else 0,
               var_sample["OV"],
               var_sample["DP"], # be careful these aren't '.'
               #split where _r is ref-allele and _a is alt-allele
               var_sample["AD"][0],
               var_sample["AD"][1],
               var_sample["PDP"],
               var_sample["PAD"][0],
               var_sample["PAD"][1],
               var_sample["US"][0],
               var_sample["US"][1],
               var_sample["DS"][0],
               var_sample["DS"][1],
               var_sample["UC"][0],
               var_sample["UC"][1],
               var_sample["DC"][0],
               var_sample["DC"][1],
               var_sample["UDC"][0],
               var_sample["UDC"][1],
               var_sample["UCC"][0],
               var_sample["UCC"][1],
               var_sample["DDC"][0],
               var_sample["DDC"][1],
               var_sample["DCC"][0],
               var_sample["DCC"][1],
               var_sample["UMO"][0],
               var_sample["UMO"][1],
               var_sample["DMO"][0],
               var_sample["DMO"][1],
               var_sample["UXO"][0],
               var_sample["UXO"][1],
               var_sample["DXO"][0],
               var_sample["DXO"][1],
               var_sample["NR"][0],
               var_sample["NR"][1],
               var_sample["MO"][0],
               var_sample["MO"][1],
               var_sample["XO"][0],
               var_sample["XO"][1],
               var_sample["XC"][0],
               var_sample["XC"][1],
               var_sample["AC"][0],
               var_sample["AC"][1],
               var_sample["MC"][0],
               var_sample["MC"][1],
               var_sample["EC"][0],
               var_sample["EC"][1],
               var_sample["PL"][0] if var_sample["PL"][0] is not None else 0,
               var_sample["PL"][1] if var_sample["PL"][0] is not None else 0,
               var_sample["PL"][2] if var_sample["PL"][0] is not None else 0])
    return ret
    #END


def build_var_header():
    """
    Carefully build the column_name, dtypes OrderedDict of a parsed VCF
    """
    # var_id will only be populated by... ew.. this one is difficult...
    # "var_id" # - need to generate these...
    meta_header = OrderedDict()
    meta_header['sample'] = str
    meta_header['chrom'] =  str
    meta_header['start'] =  np.int16
    meta_header['end'] =  np.int16
    meta_header['var_type'] = str
    meta_header['state'] = str

    # Target - made from parsing truvari/rtg information
    # "state" #

    info_header = OrderedDict()
    info_header["POP"] = bool
    info_header["VARLEN"] = np.int32
    info_header["NUMASM"] = np.int16
    
    fmt_header = OrderedDict() 
    # Categorical
    # need to categorize these... ./. 0/0 0/1 1/1 # only possibilities
    fmt_header["GT"] = np.int8 
    # fmt_header["PG"] = np.int32 # don't use..
    fmt_header["GQ"] = np.int8
    # fmt_header["PI"] = don't use
    fmt_header["OV"] = np.int8
    fmt_header["DP"] = np.int16
    #split where _r is ref-allele and _a is alt-allele
    fmt_header["AD_r"] = np.int16
    fmt_header["AD_a"] = np.int16
    fmt_header["PDP"] = np.int16
    fmt_header["PAD_r"] = np.int16
    fmt_header["PAD_a"] = np.int16
    fmt_header["US_r"] = np.int16
    fmt_header["US_a"] = np.int16
    fmt_header["DS_r"] = np.int16
    fmt_header["DS_a"] = np.int16
    fmt_header["UC_r"] = np.int16
    fmt_header["UC_a"] = np.int16
    fmt_header["DC_r"] = np.int16
    fmt_header["DC_a"] = np.int16
    fmt_header["UDC_r"] = np.int16
    fmt_header["UDC_a"] = np.int16
    fmt_header["UCC_r"] = np.int16
    fmt_header["UCC_a"] = np.int16
    fmt_header["DDC_r"] = np.int16
    fmt_header["DDC_a"] = np.int16
    fmt_header["DCC_r"] = np.int16
    fmt_header["DCC_a"] = np.int16
    fmt_header["UMO_r"] = np.int16
    fmt_header["UMO_a"] = np.int16
    fmt_header["DMO_r"] = np.int16
    fmt_header["DMO_a"] = np.int16
    fmt_header["UXO_r"] = np.int16
    fmt_header["UXO_a"] = np.int16
    fmt_header["DXO_r"] = np.int16
    fmt_header["DXO_a"] = np.int16
    fmt_header["NR_r"] = np.int16
    fmt_header["NR_a"] = np.int16
    fmt_header["MO_r"] = np.int16
    fmt_header["MO_a"] = np.int16
    fmt_header["XO_r"] = np.int16
    fmt_header["XO_a"] = np.int16
    fmt_header["XC_r"] = np.int16
    fmt_header["XC_a"] = np.int16
    fmt_header["AC_r"] = np.int16
    fmt_header["AC_a"] = np.int16
    fmt_header["MC_r"] = np.int16
    fmt_header["MC_a"] = np.int16
    fmt_header["EC_r"] = np.int16
    fmt_header["EC_a"] = np.int16
    fmt_header["PL_ref"] = np.int8
    fmt_header["PL_het"] = np.int8
    fmt_header["PL_hom"] = np.int8

    ret_header = OrderedDict()
    ret_header.update(meta_header)
    ret_header.update(info_header)
    ret_header.update(fmt_header)
    return ret_header
    
def load_vcf(vcf_fn, asm_dat, state="UNK", sample='unknown'):
    """
    load the vcf, unting with the asm as we go. asm_data is what's returned from load_asm
    returns DataFrame of stuff
    """
    logging.info("Loading Variants from %s", vcf_fn)
    asm_header, asms = asm_dat
    ret_header = build_var_header()
    ret_header.update(asm_header)
    asmheadidx = list(asm_header.keys())
    ret = []
    with pysam.VariantFile(vcf_fn) as fh:
        for var in fh:
            cur_data = [sample]
            cur_data.append(var.chrom)
            cur_data.append(var.start)
            cur_data.append(var.stop)

            var_type, var_len = get_type_lens(var)
            
            cur_data.append(var_type)
            cur_data.append(state)
            
            baid = None # best aid
            blen = 0 # best aid length
            best = None # pulled data
            num_asms = 0
            for aid in var.info["AID"]:
                if aid not in asms:
                    continue
                num_asms += 1
                dat = asms[aid]
                alen = dat[asmheadidx.index("ASMLEN")]
                if alen > blen:
                    baid = aid
                    blen = alen
                    best = dat
            
            cur_data.append(var.info["POP"])
            cur_data.append(var_len)
            cur_data.append(num_asms)

            cur_data.extend(parse_format(var.samples[0]))
            
            if baid is not None:
                cur_data.extend(best)
            else:
                cur_data.extend(([0] * (len(asmheadidx) - 1)) + [""])


            ret.append(cur_data)
    logging.info("Loaded %d variants", len(ret))
    return pd.DataFrame(ret, columns = ret_header.keys())

def save_table(data, out_file):
    """
    Save to a joblib
    """
    logging.info("Saving table")
    #header, data = data
    #out = pd.DataFrame(data=data, columns = header.keys())
    joblib.dump(data, out_file)


def main(args):
    """
    Main runner
    """
    def glob_get(base_dir, fn):
        ret = glob.glob(os.path.join(base_dir, fn))
        if not len(ret) == 1:
            logging.error("Couldn't find exactly one %s/%s*", base_dir, fn)
            exit(1)
        return ret[0]

    args = parse_args(args)
    # python bgvcfasm_to_table.py
    asms = load_asm(args.asm)
    # Do this differently for truthset parsing
    variants = None
    if args.mode == 'full':
        variants = load_vcf(args.vcf, asms, "UNK", args.sample)#, sys.argv[3], sys.argv[4])
    elif args.mode == 'bench':
        tp_small = load_vcf(glob_get(args.rtgdir, "tp.vcf.gz"), asms, "TP", args.sample)
        fp_small = load_vcf(glob_get(args.rtgdir, "fp.vcf.gz"), asms, "FP", args.sample)
        tp_large = load_vcf(glob_get(args.trudir, "tp-call.vcf.gz"), asms, "TP", args.sample)
        fp_large = load_vcf(glob_get(args.trudir, "fp.vcf.gz"), asms, "FP", args.sample)
        variants = pd.DataFrame()
        variants.append(tp_small)
        for i in [tp_small, fp_small, tp_large, fp_large]:
            if not  i.empty :
                variants = variants.append(i)

        #variants = pd.concat([tp_small, fp_small, tp_large, fp_large])
        #variants = tp_large.append(tp_small)
        #variants = variants.append(fp_small)
        variants.reset_index()
    
    save_table(variants, args.out)
    logging.info("Finished")

if __name__ == '__main__':
    main(sys.argv[1:])
