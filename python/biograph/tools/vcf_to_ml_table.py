"""
Translates a VCF into a simple table with structure

AlleleId    AF  MAF AC  AC_Het  AC_Hom  AC_Hemi HWE ExecHet Sample_Depth    Sample_Depth    ...

Where AlleleId is a unique descriptor of the the variant and depth is the allele's coverage
AlleleId has structure:
    bg_{chrom}:{pos}.{alleleNumber}

To transpose the output, run the following awk command:
    awk '
    {
        for (i=1; i<=NF; i++)  {
            a[NR,i] = $i
        }
    }
    NF>p { p = NF }
    END {
        for(j=1; j<=p; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){
                str=str"\t"a[i,j];
            }
            print str
        }
    }' output.txt
"""
import sys
import gzip
import argparse

from collections import defaultdict

def parse_args(args):
    """
    Parse all arguments
    """
    parser = argparse.ArgumentParser(
        prog="vcf_to_ml_table", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', metavar='VCF', type=str, help='VCF files to tabulate')
    args = parser.parse_args(args)
    return args


def make_key(entry):
    """
    Turn vcf entry into its allele keys
    """
    yield "bg_%s:%s.%d" % (entry[0], entry[1], 0), 0
    for pos in range(len(entry[4].split(','))):
        yield "bg_%s:%s.%d" % (entry[0], entry[1], pos + 1), pos + 1


def info_to_dict(entry):
    """
    Translate the info INFO of an entry to a dict for lookup
    """
    ret = defaultdict(lambda: '.')
    for i in entry[7].split(';'):
        if i.count('='):
            key, value = i.split("=")
        else:
            key = i
            value = True
        ret[key] = value
    return ret


def main():
    """
    Main
    """

    info_keys = ["AF", "MAF", "AC", "AC_Het", "AC_Hom", "AC_Hemi", "HWE", "ExcHet"]
    with gzip.GzipFile(sys.argv[1], 'r') as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                sys.stdout.write("pos\t%s\t%s\n" % ("\t".join(info_keys), "\t".join(line.strip().split('\t')[9:])))
                continue

            entry = line.strip().split('\t')
            # Should check AD (and others) are present
            fmt = entry[8].split(':')
            my_depth = []
            for i in entry[9:]:
                my_depth.append(i.split(':')[fmt.index("AD")].split(','))

            # This is per allele, not per-sample informaiton
            my_info = info_to_dict(entry)
            # For each alternate allele, hold each of the keys below
            ex_info = []  # extracted info keys for writing
            for key in info_keys:
                ex_info.append(my_info[key].split(','))

            #AF=0.590717,0.0295359;MAF=0.409283,0.0295359;AC=140,7;AC_Het=90,0;AC_Hom=50,0;AC_Hemi=0,7;HWE=4.60983e-14,1;ExcHet=4.05062e-14,1
            for key, position in make_key(entry):
                sys.stdout.write(key)

                # This is a problem for unannotated files...
                if position == 0:
                    sys.stdout.write("\t." * len(info_keys))
                else:
                    for info in ex_info:
                        sys.stdout.write("\t%s" % (info[position - 1]))
                for samp in my_depth:
                    if len(samp) > position:
                        sys.stdout.write("\t%s" % (samp[position]))
                    else:
                        sys.stdout.write("\t.")
                sys.stdout.write("\n")

if __name__ == '__main__':
    main()
