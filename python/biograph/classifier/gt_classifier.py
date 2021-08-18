#!/usr/bin/env python3
"""
Classify Genotypes
"""
import os
import re
import sys
import types
import logging
import argparse
import tempfile
import multiprocessing

from io import StringIO
from collections import OrderedDict

import tabix
import joblib
import pandas as pd
import numpy as np

from biograph.utils import get_opener
from truvari import setup_logging
#pd.set_option('mode.chained_assignment', None)

qcls_shared = types.SimpleNamespace()

def parse_args(clargs):
    """ Command line UI """
    parser = argparse.ArgumentParser(prog="gt_classifier", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--vcf", required=True,
                        help="VCF to parse")
    parser.add_argument("-d", "--dataframe", required=True,
                        help="Coverage DataFrame frame")
    parser.add_argument("-m", "--model", required=True,
                        help="Model to apply to data")
    parser.add_argument("-o", "--out", default="/dev/stdout",
                        help="VCF to output")
    parser.add_argument("--sample", type=str, default="",
                        help="Sample identifier (only required for multi-sample VCFs)")
    parser.add_argument("--tmp", type=str, default=tempfile.gettempdir(),
                        help="Temporary directory (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    # We could put heristics for setting a GT Filter

    args = parser.parse_args(clargs)
    setup_logging()
    return args

def transform_dp(data):
    """
    Preprocessing the coverage metrics for classification
    Mainly the _r / _a depth numbers into all ratios of one another
    Returns new dataframe
    """
    features = ["VARLEN", "GT", "GQ", "DP", "AD_r", "AD_a", "US_r", "US_a",
                "DS_r", "DS_a", "UC_r", "UC_a", "DC_r", "DC_a", "UDC_r",
                "UDC_a", "UCC_r", "UCC_a", "DDC_r", "DDC_a", "DCC_r", "DCC_a",
                "UMO_r", "UMO_a", "DMO_r", "DMO_a", "NR_r", "NR_a", "MO_r", "MO_a",
                "XC_r", "XC_a", "AC_r", "AC_a", "MC_r", "MC_a", "EC_r", "EC_a",
                "PL_ref", "PL_het", "PL_hom", "UXO_r", "UXO_a", "DXO_r", "DXO_a",
                "XO_r", "XO_a"]
    trans_features = ["AD_r", "AD_a", "UC_r", "UC_a", "DC_a", "DC_r", "UDC_r", "UDC_a",
                      "UCC_r", "UCC_a", "DDC_r", "DDC_a", "DCC_r", "DCC_a", "XC_r",
                      "XC_a", "AC_r", "AC_a", "MC_r", "MC_a", "EC_r", "EC_a"]
    keep_features = [_ for _ in features if _ not in trans_features] + ["RC"]

    ret = pd.DataFrame()
    data['VARLEN'] = data.apply(set_size, axis=1)
    union = set()
    for i in trans_features:
        if i.endswith("_r") or i.endswith("_a"):
            union.add(i.split("_")[0])
        else:
            ret[i] = data[i]

    for i in sorted(list(union)):
        denom = data[i + "_a"] + data[i + "_r"]
        alt = data[i + "_a"] / denom
        ref = data[i + "_r"] / denom
        ret[i] = alt - ref

    for i in keep_features:
        ret[i] = data[i]
    # No coverage is assumed to be -1
    to_return = norm_rl(ret.fillna(-1))
    return to_return

def norm_rl(df):
    """
    Normalizes read length dependent features
    """
    rl_feat = ["US_r", "US_a", "DS_a", "DS_r", "UXO_r", "UXO_a", "DXO_r",
               "DXO_a", "UMO_r", "UMO_a", "DMO_r", "DMO_a", "MO_r", "MO_a",
               "XO_r", "XO_a"]
    rl = df['MO_r'].max()
    df[rl_feat] = df[rl_feat] / rl
    return df

def set_size(x):
    """
    Set VARLEN to the appropriate sign
    """
    return x["VARLEN"] if x["var_type"] == "INS" else -x["VARLEN"]

def run_model(data, model):
    """
    Run the gt model on the data
    """
    np.seterr(all="ignore")
    t_data = transform_dp(data)
    gts = pd.Series(model.predict(t_data), name="predict_GT", index=data.index)
    pls = pd.DataFrame(model.predict_proba(t_data), columns=["REF", "HET", "HOM"], index=data.index)
    def rounder(d):
        r = -10 * np.log10(d)
        r = r.replace(-0, 0).replace([np.inf], 99).round()
        return r
    pls['PL_ref'] = rounder(pls["REF"])
    pls['PL_het'] = rounder(pls['HET'])
    pls['PL_hom'] = rounder(pls['HOM'])
    pls["predict_GQ"] = get_gq(pls)
    ret = pls.join(gts)
    return ret

def get_gq(df):
    """
    Calculate GQ
    """
    pl_arr = df[['PL_ref', 'PL_het', 'PL_hom']].values
    pl_arr.sort(axis=1)
    gq = [(row[1] - row[0]) for row in pl_arr]
    return gq

def fmt_to_dict(fmt, samp):
    """
    zip format and sample information to make a dictionary
    """
    ret = OrderedDict()
    for key, val in zip(fmt.split(':'), samp.split(':')):
        ret[key] = val
    return ret

def dict_to_fmt(samp):
    """
    Return tuple for VCF's FORMAT\tSAMPLE columns
    """
    return ":".join(samp.keys()), ":".join(samp.values())

def make_stats():
    """
    Counter stats box
    """
    stats = {"none_cnt": 0,
             "filt_cnt": 0,
             "lowq_cnt": 0,
             "pass_cnt": 0,
             "het_cnt": 0,
             "hom_cnt": 0,
             "hom_ref": 0}
    return stats

def edit_vcf(chunk):
    """
    Edits/filters vcf based of classifications
    """
    # indexes of chunk
    # ref_name, chunk_start, chunk_stop, in_vcf_name, filter_thresh, lowqual_sv, lowqual_ao, clsf = chunk
    fdat = qcls_shared.data
    stats = make_stats()
    out_vcf = StringIO()
    in_vcf = tabix.open(chunk[3])
    try:
        region = in_vcf.query(chunk[0], chunk[1], chunk[2])
    except tabix.TabixError:
        logging.warning("No Fetch over %s", str(chunk))
        out_vcf.seek(0)
        return out_vcf, stats
    for entry in region:
        start = int(entry[1]) - 1
        stop = start + len(entry[3])
        sample_data = fmt_to_dict(entry[8], entry[9])
        key = "%s:%d-%d.%s" % (entry[0], start, stop, entry[4])

        if key not in fdat.index:
            stats["none_cnt"] += 1
            out_vcf.write("\t".join(entry) + '\n')
            continue
        predict = fdat.loc[key]
        if predict["predict_GT"] == 1:
            stats["het_cnt"] += 1
            sample_data["GT"] = "0/1"
        elif predict["predict_GT"] == 2:
            stats["hom_cnt"] += 1
            sample_data["GT"] = "1/1"
        else:
            stats["hom_ref"] += 1
            sample_data["GT"] = "0/0"

        sample_data["GQ"] = str(int(predict["predict_GQ"]))
        sample_data["PL"] = ",".join([str(int(_)) for _ in predict[["PL_ref", "PL_het", "PL_hom"]].to_list()])
        entry[8], entry[9] = dict_to_fmt(sample_data)
        out_vcf.write("\t".join(entry) + '\n')
    out_vcf.seek(0)
    return out_vcf, stats

def ref_ranges(contigs, chunk_size, in_vcf_name):
    """
    Chunk reference into pieces
    """
    for ref_name, final_stop in contigs:
        start = 0
        stop = start + chunk_size
        while stop < final_stop:
            yield ref_name, start, stop, in_vcf_name
            start = stop
            stop += chunk_size
        yield ref_name, start, final_stop, in_vcf_name

def edit_header(in_vcf):
    """
    Create an edited header
    return the header as well as contig information from vcf header
    """
    header = []
    ctgre = re.compile(r"##contig=<ID=(?P<name>.*),length=(?P<length>\d+)>$")
    contigs = []
    while True:
        line = in_vcf.readline().decode()
        header.append(line)
        mat = ctgre.match(line)
        if mat:
            mat = mat.groupdict()
            contigs.append((mat["name"], int(mat["length"])))
        if line.startswith("#CHROM"):
            num_cols = line.strip().split('\t')
            if len(num_cols) != 10:
                logging.error(f"Input VCF doesn't have exactly 10 columns (found {num_cols})")
                logging.error("Input VCF must only have a single sample")
                exit(1)
            break
    return header, contigs

def write_vcf(results, input_vcf, output_vcf, threads=1):
    """
    Write results to output vcf
    """
    # edit header
    in_vcf = get_opener(input_vcf)
    header, contigs = edit_header(in_vcf)

    if not contigs:
        logging.error("No contigs found in VCF header")
        exit(1)

    qcls_shared.data = results

    out_vcf = open(output_vcf, 'w')
    out_vcf.write("".join(header))
    logging.info("Writing Results")
    stats = make_stats()
    with multiprocessing.Pool(threads) as pool:
        iters = ref_ranges(contigs, chunk_size=100000000, in_vcf_name=input_vcf)
        for out, chunk_stats in pool.map(edit_vcf, iters):
            out_vcf.write(out.read())
            for key in chunk_stats:
                stats[key] += chunk_stats[key]
    out_vcf.close()
    ratio = stats["het_cnt"] / stats["hom_cnt"] if stats["hom_cnt"] != 0 else 0
    logging.info(f'Reference {stats["hom_ref"]}; Homozygous {stats["hom_cnt"]}; Heterozygous {stats["het_cnt"]}; Ratio {ratio:.2f}')

def main(args):
    """
    Classify genotypes by coverage metrics
    """
    args = parse_args(args)
    if not os.path.exists(args.vcf):
        logging.error(f"VCF {args.vcf} is missing")
        exit(1)
    logging.info(f"Loading model {args.model}")
    model = joblib.load(args.model)['gtcls_model']
    logging.info(f"Loading {args.dataframe}")
    # Load the model
    cov_data = joblib.load(args.dataframe)
    # Load the coverage info
    cov_data = cov_data[cov_data["VARLEN"] >= 50]
    cov_data = cov_data.sort_values('key')
    cov_data = cov_data.drop_duplicates('key')
    cov_data = cov_data.set_index('key')
    # Run the model
    results = run_model(cov_data, model)
    # Write the VCF
    write_vcf(results, args.vcf, args.out, args.threads)
    logging.info("Finished")

if __name__ == '__main__':
    main(sys.argv[1:])
