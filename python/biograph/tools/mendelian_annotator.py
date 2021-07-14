"""
Annotates a trio vcf's INFO field with a 'szcat'egory and Mendelian Error Families (MEF)
Outputs a json and txt summary of these numbers
"""
import os
import re
import sys
import json
import gzip
import argparse
from collections import OrderedDict
from biograph.utils import Pedigree
import biograph.tools.log as log

# TODOs - pedigree files are broken when parents are present because
# I automatically add them and then it does a weird double check


class MendelReport:

    """
    Holds all the stats from counting mendelian errors
    takes a Pedigree object and a dict with {sample_id: column_in_vcf}
    Writes a final report
    """

    def __init__(self, ped, samples):
        self.ped = ped
        self.samples = samples
        self.__report = {}
        self.tot_count = 0
        self.__depth_cats = ["All", 0, 1, 10, 15]
        self.__sv_lens = ["All", "1-9", "10-49", "50-99", "100-299", "300-999", "gt1000"]
        self.__sv_types = ["DEL", "INS", "SUBSDEL", "SUBSINS", "REPL"]
        self.trio_columns = []

        self.__make_trio_vcf_columns()
        self.__make_blank_report()

    def __make_trio_vcf_columns(self):
        """
        create a list of tuples of all families' trios tied to the vcf samples' columns
        tuples fields : fam_id, proband_id_column, father_id_column, mother_id_column
        """
        for indiv in self.ped.get_trio_probands():
            if indiv.ind_id not in self.samples:
                log.warning("Individual %s in pedigree, but not in vcf", indiv.ind_id)
                continue
            if indiv.pat_id not in self.samples:
                log.warning("Individual %s is missing father %s in vcf", indiv.ind_id, indiv.pat_id)
                continue
            if indiv.mat_id not in self.samples:
                log.warning("Individual %s is missing mother %s in vcf", indiv.ind_id, indiv.mat_id)
                continue
            self.trio_columns.append((indiv.fam_id, self.samples[indiv.ind_id],
                                      self.samples[indiv.pat_id], self.samples[indiv.mat_id]))

    def __make_blank_report(self):
        """
        Setup my blank report
        """

        def make_sub_report():
            """
            Single collection of a report
            ret[a grouping]      # a single collection of numbers e.g. overall or a single family
                ["var_count"]     # number of variants looked at
                ["type_counts"]   # 2d of sv_type x tot|err|pct
                    [sv_type] [tot | err | pct]
                ["size_type_counts"] # 3d of size x sv_type x tot|err|pct
                    [size] [sv_type] [ tot | err | pct]
                ["cov_size_me"] # 3d of size x sv_type x tot|err| pct
                    [depth_cat] [sv_len] [ tot | err | pct ]
            """
            type_counts = OrderedDict()
            for sv_type in self.__sv_types:
                type_counts[sv_type] = {"tot": 0, "err": 0, "pct": 0.0}

            size_type_counts = OrderedDict()
            for sv_len in self.__sv_lens:
                size_type_counts[sv_len] = OrderedDict()
                for sv_type in self.__sv_types:
                    size_type_counts[sv_len][sv_type] = {"tot": 0, "err": 0, "pct": 0.0}

            cov_size_me = OrderedDict()
            for depth_cat in self.__depth_cats:
                cov_size_me[depth_cat] = OrderedDict()
                for sv_len in self.__sv_lens:
                    cov_size_me[depth_cat][sv_len] = {"tot": 0, "err": 0, "pct": 0.0}

            return {
                "type_counts": type_counts,
                "size_type_counts": size_type_counts,
                "cov_size_me": cov_size_me,
                "var_count": 0
            }

        # End make_sub_report
        if "overall" in self.ped.families:
            log.error("Cannot use Family id 'overall'. Rename them in pedigree file")
            exit(1)
        for i in self.trio_columns:
            self.__report[i[0]] = make_sub_report()
        self.__report["overall"] = make_sub_report()

    def increment(self, group, sv_type, sv_len, depth_cat, is_err):
        """
        increment all the values based on the properties
        """
        self.__report[group]["type_counts"][sv_type]["tot"] += 1
        self.__report[group]["size_type_counts"]["All"][sv_type]["tot"] += 1
        self.__report[group]["size_type_counts"][sv_len][sv_type]["tot"] += 1
        # depth_cat = None will be overall where I know to just report the average
        # over all families
        self.__report[group]["cov_size_me"][depth_cat][sv_len]["tot"] += 1
        self.__report[group]["cov_size_me"][depth_cat]["All"]["tot"] += 1
        self.__report[group]["cov_size_me"]["All"][sv_len]["tot"] += 1
        self.__report[group]["cov_size_me"]["All"]["All"]["tot"] += 1

        if is_err:
            self.__report[group]["type_counts"][sv_type]["err"] += 1
            self.__report[group]["size_type_counts"]["All"][sv_type]["err"] += 1
            self.__report[group]["size_type_counts"][sv_len][sv_type]["err"] += 1

            self.__report[group]["cov_size_me"][depth_cat][sv_len]["err"] += 1
            self.__report[group]["cov_size_me"][depth_cat]["All"]["err"] += 1
            self.__report[group]["cov_size_me"]["All"][sv_len]["err"] += 1
            self.__report[group]["cov_size_me"]["All"]["All"]["err"] += 1

    def consolidate_report(self):
        """
        Given a filled out reports with counts, calculate the merr pcts
        """

        def pctify(table_row):
            """
            populate the pct of a category
            """
            if table_row["tot"] != 0:
                table_row["pct"] = table_row["err"] / float(table_row["tot"])

        # special case for the overall - only keep the All
        self.__report["overall"]["cov_size_me"] = {"All": self.__report["overall"]["cov_size_me"]["All"]}
        for group in self.__report.values():
            for type_col in group["type_counts"].values():
                pctify(type_col)

            for size_x in group["size_type_counts"].values():
                for type_y in size_x.values():
                    pctify(type_y)

            for depth_x in group["cov_size_me"].values():
                for size_y in depth_x.values():
                    pctify(size_y)

    def write_json(self, out_name=None):
        """
        Dump my report to a json
        """
        if out_name is None:
            out = sys.stdout
        else:
            out = open(out_name, 'w')
        json.dump(self.__report, out)
        out.close()

    def write(self, out_name=None):
        """
        given a json file of a mendel_anno report,
        convert to a simplier, tab-delimited file
        """
        if out_name is None:
            out = sys.stdout
        else:
            out = open(out_name, 'w')

        def write_group(grouping):
            """
            output all the stats of an individual grouping
            """
            out.write("Type\tTot\tErr\tPct\n")
            for key in self.__sv_types:
                out.write(
                    "%s\t%d\t%d\t%.2f\n" % (key, grouping["type_counts"][key]['tot'],
                                            grouping["type_counts"][key]['err'], grouping["type_counts"][key]['pct']))

            out.write("\nCounts By Size and Type\n")
            group_cat = grouping["size_type_counts"]
            for title, fmt in [("tot", "\t%d"), ("err", "\t%d"), ("pct", "\t%.2f")]:
                out.write("%s\t%s\n" % (title.upper(), \
                  "\t".join(next(iter(group_cat.values())).keys())))
                for size in group_cat.keys():
                    out.write(size)
                    for sv_type in group_cat[size]:
                        out.write(fmt % group_cat[size][sv_type][title])
                    out.write("\n")
                out.write("\n")
            out.write("\nMendelian Error Table\n")
            out.write("Cov\tSize\tTotal\tErrors\tMendelianError\n")
            group_cat = grouping["cov_size_me"]
            for row1 in group_cat.keys():
                for row2 in group_cat[row1]:
                    out.write("{0}\t{1}\t{tot}\t{err}\t{pct:.4}\n".format(row1, row2,
                                                                          **grouping["cov_size_me"][row1][row2]))

        # End write_group
        out.write("var_count\t%d\n" % (self.tot_count))
        out.write("========\nOverall\n========\n")
        write_group(self.__report["overall"])
        for group in self.__report:
            if group != "overall":
                out.write("\n%s\n%s\n%s\n" % ("=" * len(group), group, "=" * len(group)))
                write_group(self.__report[group])
        out.close()


def setup_output_vcf(outname, t_vcf):
    """
    Create an output vcf.Writer given the input vcf file as a templte
    writes the full header and
    Adds info fields:
        sizeCat
        MEF
    Returns a file handler and a dict with {individual_id: column in vcf}
    """
    out = open(outname, 'w')
    line = t_vcf.readline()
    samp_columns = {}
    while not line.startswith("#CHROM"):
        out.write(line)
        line = t_vcf.readline()

    # edit the header
    out.write('##INFO=<ID=sizeCat,Number=A,Type=String,Description="Size category of variant">\n')
    out.write('##INFO=<ID=MEF,Number=.,Type=String,Description="Names of families that contain mendelian error">\n')
    out.write(line)
    for pos, iid in enumerate(line.strip().split('\t')[9:]):
        samp_columns[iid] = pos + 9
    return out, samp_columns


def get_type_lens(entry):
    """
    Parse an entry and return it's sv_type and it's sv_len
    """
    mREF = entry[3]
    # TODO - should get the longest?
    mALTs = entry[4].split(',')
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
            log.error(str(entry))
            log.error("shouldn't have some new crazy type\n")
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

    return sv_types, ret_lens, sv_lens


def get_coverage_category(format_keys, indivs):
    """
    Parse a vcf entry and set it's coverage category based on the min depth
    of all samples
    """
    coverages = [0, 1, 10, 15]
    depth_cat = None
    idx = format_keys.index("DP")
    get_cov = lambda x: None if x[idx] == '.' else int(x[idx])

    for min_coverage in coverages:
        min_covered = True
        for sample in indivs:
            data = sample.split(':')
            cov = get_cov(data)
            if cov is None:
                continue
            if cov < min_coverage:
                min_covered = False
        if min_covered:
            depth_cat = min_coverage
    return depth_cat


def me_check(pr, fa, ma, ref_count=False):
    """
    Simplest possible way to check ME
    -1 for skipped sites ('.' or hom-ref in all samples
    0 for consistent
    1 for inconsistent
    """
    pr = pr.split(':')[0].split('/')
    fa = fa.split(':')[0].split('/')
    ma = ma.split(':')[0].split('/')
    if '.' in pr + fa + ma:
        return -1
    null = ["0", "0"]
    if not ref_count and pr == null and fa == null and ma == null:
        return -1
    if pr[0] in fa and pr[1] in ma:
        return 0
    if pr[1] in fa and pr[0] in ma:
        return 0
    return 1


def check_family_mendel(entry, sv_type, sv_len, ret_report, ref_count):
    """
    For an entry, iterate over every full family in a pedigree and
    check their mendelian consistency.
    Also add to the report per-family AND overall ME stats
    return a list of family ids that have a mendelian error
    """
    format_keys = entry[8].split(':')

    error_families = []
    # fam_id, pro_id_col, fat_id_col, mat_id_col
    for family in ret_report.trio_columns:
        fam = family[0]
        # p_dat, f_dat, m_dat
        samp_dat = [entry[family[1]],
                    entry[family[2]],
                    entry[family[3]]]
        m_error = me_check(*samp_dat, ref_count=ref_count)
        if m_error == -1:  # This does happen with the ref_count now. Need to count skipped?
            continue

        depth_cat = get_coverage_category(format_keys, samp_dat)
        ret_report.increment(fam, sv_type, sv_len, depth_cat=depth_cat, is_err=m_error)
        if m_error:
            error_families.append(fam)

    return error_families


def annotate_vcf(m_vcf, output_vcf, ped, ref_count):
    """
    Separate merged vcf by size and then calculate mendelian error
    Create the annotated vcf output
    and returns the json report
    """
    log.info("Annotating VCF")

    cur_vcf = None
    # pylint has a false positive here
    # pylint: disable=useless-suppression
    if m_vcf == '-':
        cur_vcf = sys.stdin
    elif m_vcf.endswith(".gz"):
        cur_vcf = gzip.GzipFile(m_vcf, 'r')
    else:
        cur_vcf = open(m_vcf, 'r')

    out_vcf, samples = setup_output_vcf(output_vcf, cur_vcf)

    ret_report = MendelReport(ped, samples)

    for entry in cur_vcf:
        entry = entry.strip().split('\t')
        sv_info = get_type_lens(entry)

        ret_report.tot_count += 1
        types = ",".join(sv_info[0])
        entry[7] += ";sizeCat=%s" % (",".join(sv_info[1]))
        sizes = ",".join([str(x) for x in sv_info[2]])

        has_sv = [x for x in sv_info[2] if x >= 25]
        if "SVTYPE" in entry[7]:
            entry[7] = re.sub(r"SVTYPE=\d+", "SVTYPE=%s" % types, entry[7])
        elif has_sv:
            entry[7] += ";SVTYPE=%s" % (types)

        if "SVLEN" in entry[7]:
            entry[7] = re.sub(r"SVLEN=\d+", "SVLEN=%s" % sizes, entry[7])
        elif has_sv:
            entry[7] += ";SVLEN=%s" % (sizes)

        err_fams = check_family_mendel(entry, sv_info[0][0], sv_info[1][0], ret_report, ref_count)
        ret_report.increment("overall", sv_info[0][0], sv_info[1][0], depth_cat=0, is_err=len(err_fams))

        if err_fams:
            entry[7] += ";MEF=%s" % (",".join(err_fams))

        out_vcf.write("\t".join(entry) + '\n')

    ret_report.consolidate_report()
    return ret_report


def parse_args(args):
    """
    Parse the args
    """
    parser = argparse.ArgumentParser(prog="meanno", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--vcf", type=str, required=True,
                        help="VCF file to annotate (use '-' for uncompressed stdin)")
    parser.add_argument("-p", "--pedigree", type=str, required=True,
                        help="Pedigree file")
    parser.add_argument("-r", "--ref-count", action="store_true", default=False,
                        help="Count homozygous reference in trios as mendelian consistent (%(default)s)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output vcf file to write")
    args = parser.parse_args(args)
    log.setup_logging()

    # Get absolute paths of everything
    if args.vcf != "-":
        args.vcf = os.path.abspath(args.vcf)
    args.output = os.path.abspath(args.output)
    return args


def main(args):
    """
    Run this as a tool
    """
    args = parse_args(args)
    ped = Pedigree(args.pedigree)

    m_report = annotate_vcf(args.vcf, args.output, ped, args.ref_count)
    m_report.write_json(args.output.rstrip(".vcf.gz") + "_stats.json")
    m_report.write(args.output.rstrip(".vcf.gz") + "_stats.txt")

    log.info("Finished")


if __name__ == '__main__':
    main(sys.argv[1:])
