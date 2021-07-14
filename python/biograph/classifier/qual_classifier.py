#!/usr/bin/env python3
"""
Classify VCF variants
"""
import re
import sys
import types
import argparse
import tempfile
import multiprocessing

from io import StringIO
from collections import OrderedDict

import numpy
import tabix
import joblib
import pandas as pd

import biograph.tools.log as log
from biograph.utils import get_opener
pd.set_option('mode.chained_assignment', None)

qcls_shared = types.SimpleNamespace()

CLSCHOICES = {"GT": 1,
              "Qual": 2,
              "Both": 3}

def parse_args(clargs):
    """ Command line UI """
    parser = argparse.ArgumentParser(prog="qual_classifier", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--vcf", required=True,
                        help="VCF to parse")
    parser.add_argument("-d", "--dataframe", required=True,
                        help="Coverage DataFrame frame")
    parser.add_argument("-m", "--model", required=True,
                        help="Model to apply to data")
    parser.add_argument("-o", "--out", default="/dev/stdout",
                        help="VCF to output")
    parser.add_argument("-x", "--grm", default=None,
                        help="DataFrame conaining grm features from truvari")
    parser.add_argument("-f", "--filter", default=0.10, type=float,
                        help="Maximum threshold of calls to filter (%(default)s)")
    parser.add_argument("-s", "--lowqual_sv", default=0.352, type=float,
                        help="Maximum threshold for calls to mark as lowqual_sv (%(default)s)")
    parser.add_argument("-a", "--lowqual_ao", default=0.22, type=float,
                        help="Maximum threshold for calls to mark as lowqual_ao (%(default)s)")
    parser.add_argument("--sample", type=str, default="",
                        help="Sample identifier (only required for multi-sample VCFs)")
    parser.add_argument("--tmp", type=str, default=tempfile.gettempdir(),
                        help="Temporary directory (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of threads to use (%(default)s)")
    parser.add_argument("-g", "--thresh_gt", default=0.75, type=float,
                        help="threshold for GT")
    parser.add_argument("-c", "--clsf", default="Both", type=str, choices=CLSCHOICES,
                        help="Flag for which classifiers to run (%(default)s)")

    args = parser.parse_args(clargs)
    args.clsf = CLSCHOICES[args.clsf]
    log.setup_logging()
    if args.clsf & 2 and args.grm is None:
        log.error("Expecting --grm dataframe for this classifier run. Exiting")
        exit(1)

    return args

# @profile
def prepare_input(df):
    """
    Encodes categorical values and splits SVs and AO

    NOTE: Modifies df!
    """
    df.drop(['sample'], axis=1, inplace=True)

    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, float('inf')]
    labels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    df['DP_Bins'] = pd.cut(df['DP'], bins, labels=labels)

    POP_d = {False:0, True:1}
    df['POP_True'] = df['POP'].map(POP_d).fillna(df['POP'])

    var_type_d = {'REPL':0, 'SUBSINS':1, 'INS':2, 'SUBSDEL':3, 'DEL':4}
    df['VAR_TYPE'] = df['var_type'].map(var_type_d).fillna(0)

    type_d = {'DEL':-1, 'INS':1, 'SUBSDEL':-1, 'SUBSINS':1}
    df['VARTYPE'] = df['var_type'].map(type_d).fillna(0)
    df['VAR_LEN_TYPE'] = df['VARLEN']*df['VARTYPE']

    ret_sv = df[abs(df["VARLEN"]) >= 20]
    ret_sv.reset_index(inplace=True, drop=True)
    #ret_sv = norm_rl(ret_sv)

    ret_ao = df[abs(df["VARLEN"]) < 20]
    ret_ao.reset_index(inplace=True, drop=True)

    return ret_sv, ret_ao

def norm_rl(df):
    """
    Normalizes read length dependent features
    """
    rl_feat = ["US_r", "US_a", "DS_a", "DS_r", "UXO_r", "UXO_a", "DXO_r", "DXO_a", "UMO_r", "UMO_a", "DMO_r", "DMO_a", "MO_r", "MO_a", "XO_r", "XO_a"]
    rl = df['MO_r'].max()
    df[rl_feat] = df[rl_feat]/rl
    return df

# @profile
def run_classifier(classifier, data, features, typ):
    """
    Run a classfier on data with appropriate features
    adds a "probs" column
    """
    to_classify = data.dropna(subset=features)
    probs = pd.DataFrame(classifier.predict_proba(to_classify[features]),
                         columns=["fp", "tp"], index=to_classify.index)
    probs["type"] = typ
    return pd.concat([data, probs], axis=1)

def run_gt(classifier, data, features=None, thresh=0.50):
    """
    Run gt classification - drops then replaces "GT" column
    """
    data = data.drop(["GT"], axis=1)
    to_classify = data.dropna(subset=features)
    gt_prob = classifier.predict_proba(to_classify[features])
    gt_pred = pd.Series([1 if i > thresh else 2 for i in gt_prob[:, 1]],
                        name="GT", index=to_classify.index)

    return pd.concat([data, gt_pred], axis=1)

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
             "hom_cnt": 0}
    return stats

# @profile
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
        out_vcf.seek(0)
        return out_vcf, stats
    filter_thresh, lowqual_sv, lowqual_ao, clsf = chunk[4:]
    for entry in region:
        start = int(entry[1]) - 1
        stop = start + len(entry[3])
        qual = entry[5]
        filt = entry[6]
        sample_data = fmt_to_dict(entry[8], entry[9])
        key = "%s:%d-%d.%s" % (entry[0], start, stop, entry[4])

        if key not in fdat.index:
            entry[6] = "lowq"
            stats["none_cnt"] += 1
            out_vcf.write("\t".join(entry) + '\n')
            continue

        if clsf & 1:
            if fdat["GT"][key] == 1:
                stats["het_cnt"] += 1
                sample_data["GT"] = "0/1"
            else:
                stats["hom_cnt"] += 1
                sample_data["GT"] = "1/1"

        if clsf & 2:
            score = fdat["tp"][key]
            if fdat["type"][key] == "sv":
                lowqual = lowqual_sv
                filter_thresh = filter_thresh
            else:
                lowqual = lowqual_ao
                filter_thresh = lowqual_ao - 0.05

            if score < filter_thresh:
                stats["filt_cnt"] += 1
                continue
            qual = str(int(score * 100) if not numpy.isnan(score) else 0)
            if score < lowqual:
                filt = "lowq"
                stats["lowq_cnt"] += 1
            else:
                stats["pass_cnt"] += 1
        entry[5] = qual
        entry[6] = filt
        entry[8], entry[9] = dict_to_fmt(sample_data)
        out_vcf.write("\t".join(entry) + '\n')
    out_vcf.seek(0)
    return out_vcf, stats

def ref_ranges(contigs, chunk_size, in_vcf_name, filter_thresh, lowqual_sv, lowqual_ao, clsf):
    """
    Chunk reference into pieces
    """
    for ref_name, final_stop in contigs:
        start = 0
        stop = start + chunk_size
        while stop < final_stop:
            yield ref_name, start, stop, in_vcf_name, filter_thresh, lowqual_sv, lowqual_ao, clsf
            start = stop
            stop += chunk_size
        yield ref_name, start, final_stop, in_vcf_name, filter_thresh, lowqual_sv, lowqual_ao, clsf

def edit_header(in_vcf):
    """
    Create an edited header
    return the header as well as contig information from vcf header
    """
    header = []
    has_lowq = False
    lowq_line = '##FILTER=<ID=lowq,Description="Low confidence variants">\n'
    ctgre = re.compile(r"##contig=<ID=(?P<name>.*),length=(?P<length>\d+)>$")
    contigs = []
    while True:
        line = in_vcf.readline().decode()
        header.append(line)
        if not has_lowq and line == lowq_line:
            has_lowq = True
        mat = ctgre.match(line)
        if mat:
            mat = mat.groupdict()
            contigs.append((mat["name"], int(mat["length"])))
        if line.startswith("#CHROM"):
            num_cols = line.strip().split('\t')
            if len(num_cols) != 10:
                log.error(f"Input VCF doesn't have exactly 10 columns (found {num_cols})")
                log.error("Input VCF must only have a single sample")
                exit(1)
            if not has_lowq:
                header.insert(-1, lowq_line)
            break
    return header, contigs

# @profile
def main(args):
    ''' Assign quality scores and filter variants '''
    # NOTE: ^^^ this ^^^ docstring is used for the command description in __main__.py
    args = parse_args(args)

    log.info(f"Loading {args.dataframe}")
    data = joblib.load(args.dataframe)
    # Backwards compatibility with pre 6.0.2dev.
    if isinstance(data, dict):
        data = data["data"]
    model = joblib.load(args.model)

    data_sv, data_ao = prepare_input(data)

    # edit header
    in_vcf = get_opener(args.vcf)
    header, contigs = edit_header(in_vcf)

    if not contigs:
        log.error("No contigs found in VCF header")
        exit(1)

    # Add GT
    if args.clsf & 1:
        log.info("Classifying GT")
        data_sv = run_gt(model["model_gt"], data_sv, model["features_gt"], thresh=args.thresh_gt)
        if args.clsf == 1:
            data_sv["tp"] = 1

    if args.clsf & 2:
        # Extend dataframe with grms
        data_grm = joblib.load(args.grm)
        data_sv = norm_rl(data_sv)
        data_sv = pd.merge(data_sv, data_grm, on='key', how='inner')
        del(data_grm)

        log.info("Classifying SV, SNPs & INDELS")
        data_sv = run_classifier(model["model_sv"], data_sv, model["features_sv"], "sv")

        data_ao = run_classifier(model["model_ao"], data_ao, model["features_ao"], "ao")
    data = pd.concat([data_sv, data_ao])
    data = data.sort_values(["VARLEN", "tp"], ascending=False).drop_duplicates('key').sort_index()
    qcls_shared.data = data.set_index('key')

    out_vcf = open(args.out, 'w')
    out_vcf.write("".join(header))
    stats = make_stats()
    log.info("Writing Results")
    with multiprocessing.Pool(args.threads) as pool:
        iters = ref_ranges(contigs, chunk_size=100000000, in_vcf_name=args.vcf, filter_thresh=args.filter, \
             lowqual_sv=args.lowqual_sv, lowqual_ao=args.lowqual_ao, clsf=args.clsf)
        for out, chunk_stats in pool.map(edit_vcf, iters):
            out_vcf.write(out.read())
            for key in chunk_stats:
                stats[key] += chunk_stats[key]
    out_vcf.close()
    if args.clsf & 1:
        ratio = stats["het_cnt"] / stats["hom_cnt"] if stats["hom_cnt"] != 0 else 0
        log.info(f'Homozygous {stats["hom_cnt"]}; Heterozygous {stats["het_cnt"]}; Ratio {ratio:.2f}')

    if args.clsf & 2:
        log.info(f'None {stats["none_cnt"]}; Removed {stats["filt_cnt"]}; LowQ {stats["lowq_cnt"]}; PASS {stats["pass_cnt"]}')
    log.info("Finished")

if __name__ == '__main__':
    main(sys.argv[1:])
