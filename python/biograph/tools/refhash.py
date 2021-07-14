#!/usr/bin/env python3
''' sha256 the sorted contig:len from VCF, FASTA, SAM, GFF, or karyotypes.json '''
import argparse
import hashlib
import json
import os
import re
import sys

from biograph.utils import get_opener
from biograph.internal.refhashes import (
    refhashes, refstyle
)

__all__ = [
    'refhash', 'refstyle'
]

class refhash():
    """
    Container for refhash operations. Initialize it with a filename, a data=list/tuple of strings,
    or a lookup of hash/common name.
    """
    def __init__(self, filename=None, data=None, lookup=None):
        """ Someone set us up the refhash """
        if filename and data:
            raise RuntimeError('Specify either filename or data, not both.')

        if filename and os.path.isdir(filename):
            if os.path.isfile(f'{filename}/karyotype.json'):
                self.filename = f'{filename}/karyotype.json'
            else:
                raise SystemExit('The specified directory is not a BioGraph reference.')
        else:
            self.filename = filename

        self.data = data
        self.vcf_refhash = None
        self.contigs = []
        self.the_digest = None

        # specify refname to skip computation and look up the hash instead
        if lookup:
            if lookup in refhashes:
                self.the_digest = lookup
            else:
                # Try to look it up by common name, or raise on fail
                self.the_digest = refhash.get_hash(lookup)

        elif filename or data:
            self.the_digest = self.compute_hash()

    def digest(self):
        """ Accessor function for symmetry """
        return self.the_digest

    def common_name(self, the_hash=None):
        ''' Return the common name for any digest we've seen before, otherwise unknown-[...] '''
        digest = the_hash or self.the_digest
        if digest in refhashes:
            return refhashes[digest]['common']
        return f'unknown-{digest[0:8]}'

    def verbose_name(self, the_hash=None):
        ''' Return a vcf-like identifier '''
        digest = the_hash or self.the_digest
        return f"refhash={digest},name={self.common_name(digest)}"

    def build(self, the_hash=None):
        ''' Return the build '''
        digest = the_hash or self.the_digest

        if digest in refhashes:
            return refhashes[digest]['build']
        return 'unknown'

    def style(self, the_hash=None):
        ''' Return the style '''
        digest = the_hash or self.the_digest
        if digest in refhashes:
            return refhashes[digest]['style']
        return refstyle.UNKNOWN

    def url(self, the_hash=None):
        ''' Return the url '''
        digest = the_hash or self.the_digest
        if digest in refhashes:
            return refhashes[digest]['url']
        return 'unknown'

    def info(self, the_hash=None):
        ''' Return the info '''
        digest = the_hash or self.the_digest
        if digest in refhashes:
            return refhashes[digest]['info']
        return 'unknown'

    def md5(self, the_hash=None):
        ''' Return the md5 '''
        digest = the_hash or self.the_digest
        if digest in refhashes:
            return refhashes[digest]['md5']
        return 'unknown'

    def sha1(self, the_hash=None):
        ''' Return the sha1 '''
        digest = the_hash or self.the_digest
        if digest in refhashes:
            return refhashes[digest]['sha1']
        return 'unknown'

    def full_info(self, the_hash=None):
        ''' Return paragraph of info text '''
        digest = the_hash or self.the_digest

        reply = [f"refhash: {digest}"]

        if digest in refhashes:
            for entry in refhashes[digest]:
                if entry == 'style':
                    reply.append(f"{entry:>7}: {refhashes[digest][entry].name}")
                else:
                    reply.append(f"{entry:>7}: {refhashes[digest][entry]}")
        else:
            reply.append(f" common: {self.common_name(digest)}")

        return '\n'.join(reply)

    @staticmethod
    def parse_vcf_refhash(line):
        ''' Return a (digest, name) named tuple from a VCF refhash line '''
        pattern = r'^##refhash=(.*),name=(.*)$'
        match = re.search(pattern, line)
        if match:
            return {'refhash': match.group(1), 'name': match.group(2)}
        return None

    @staticmethod
    def get_hash(common):
        ''' Look up a refhash by its common name '''
        for k in refhashes:
            if refhashes[k]['common'] == common:
                return k

        raise RuntimeError(f"Could not find refhash for '{common}'")

    def get_pattern(self, line):
        '''
        Return the appropriate regex

        It should include two backreferences, where \1 is the contig name and \2 is its length.
        '''
        patterns = {
            # VCF
            '#': r'^##contig=.*ID=(.*),length=(\d+)',
            # SAM header, eg. 'samtools view -H'
            '@': r'^@SQ\tSN:(.*)\tLN:(\d+)',
            # FASTA
            '>': r'^>(\S+)\s+.*LN:(\d+)',
            # GFF
            'GFF': r'^##sequence-region\s+(\S+)\s+\d+\s+(\d+)'
        }

        if 'gff-version' in line:
            return patterns['GFF']

        if line[0] in patterns:
            return patterns[line[0]]

        raise SystemExit(f"refhash could not parse file: {self.filename}")

    def gather_contigs(self): # pylint: disable=too-many-statements
        ''' Extract a dict of {chr:len} from various formats '''
        contigs = dict()
        pattern = None
        count = 0
        current_contig = None

        # BioGraph directories cannot be streamed on STDIN.
        if self.filename and self.filename.endswith('karyotype.json'):
            with open(self.filename, 'r') as f:
                kt = json.load(f)
                for c in kt['chromosomes']:
                    contigs[c['name']] = c['len']

                self.contigs = sorted(contigs.items(), key=lambda x: x[0])
                return self.contigs

        # data can be a list or tuple of strings
        if self.data:
            lines = self.data

        # file can be optionally gzipped
        else:
            lines = get_opener(self.filename)

        for line in lines:
            if self.filename:
                try:
                    line = line.decode('utf-8')
                except UnicodeDecodeError:
                    raise SystemExit('Compressed data on STDIN is not supported. Pipe through zcat (or use -i) and try again.')

            if line.startswith('##refhash='):
                self.vcf_refhash = self.parse_vcf_refhash(line)

            if not pattern:
                pattern = self.get_pattern(line)

            # stop parsing VCF after the first non-comment
            if pattern.startswith('^#') and not line.startswith('#'):
                break

            match = re.search(pattern, line)
            if match:
                contigs[match.group(1)] = int(match.group(2))
                # close out any pending FASTA contig
                if current_contig:
                    contigs[current_contig] = count
                    current_contig = None
                    count = 0

            # FASTA may not reliably include LN at all, or even on every
            # contig line (Homo_sapiens_assembly38.fasta, I'M LOOKING AT YOU!)
            # so we revert to counting bases instead.
            elif pattern[1] == '>':
                if line[0] == '>':
                    # >anything_then_space_or_newline <OK>
                    match = re.search(r'^>(\S+)', line)
                    if not match:
                        raise SystemExit(f'Unparseable FASTA contig: {line}')

                    if current_contig:
                        contigs[current_contig] = count
                        current_contig = match[1]
                        count = 0
                        continue
                    else:
                        current_contig = match[1]
                        continue
                else:
                    if current_contig:
                        count = count + len(line.rstrip())

        # save the last contig
        if current_contig:
            contigs[current_contig] = count

        if not contigs:
            # Non-fatal warning since some VCFs may not contain contig lines
            if pattern and pattern.startswith("^##contig="):
                print(f"Warning: no contigs present in {self.filename}", file=sys.stderr)
            else:
                raise SystemExit(f"refhash could not parse file: {self.filename}")

        if self.filename:
            lines.close()

        self.contigs = sorted(contigs.items(), key=lambda x: x[0])
        return self.contigs

    def compute_hash(self):
        ''' Compute a sha256 over the contigs '''
        self.gather_contigs()
        m = hashlib.sha256()
        for k, v in self.contigs:
            # Delimit chr and len to distinguish unlikely collisions, eg.
            # chr1 249250621 vs. chr12 49250621
            # chr1:249250621\n
            m.update(f"{k}:{v}\n".encode())

        digest = m.hexdigest()

        if self.vcf_refhash:
            if digest != self.vcf_refhash['refhash']:
                raise SystemExit(f"VCF contains an invalid refhash. Computed:\n{digest}\nvs. ##refhash=\n{self.vcf_refhash['refhash']}")

        return digest

    # ref style mapping

    @staticmethod
    def to_ebi(chrom, build=None):
        ''' Convert major chromosomes to EBI style '''
        chroms = {
            'GRCh38': {
                'CM000663.2': '1',
                'CM000664.2': '2',
                'CM000665.2': '3',
                'CM000666.2': '4',
                'CM000667.2': '5',
                'CM000668.2': '6',
                'CM000669.2': '7',
                'CM000670.2': '8',
                'CM000671.2': '9',
                'CM000672.2': '10',
                'CM000673.2': '11',
                'CM000674.2': '12',
                'CM000675.2': '13',
                'CM000676.2': '14',
                'CM000677.2': '15',
                'CM000678.2': '16',
                'CM000679.2': '17',
                'CM000680.2': '18',
                'CM000681.2': '19',
                'CM000682.2': '20',
                'CM000683.2': '21',
                'CM000684.2': '22',
                'CM000685.2': 'X',
                'CM000686.2': 'Y',
                'J01415.2': 'MT',
                'NC_000001.11': '1',
                'NC_000002.12': '2',
                'NC_000003.12': '3',
                'NC_000004.12': '4',
                'NC_000005.10': '5',
                'NC_000006.12': '6',
                'NC_000007.14': '7',
                'NC_000008.11': '8',
                'NC_000009.12': '9',
                'NC_000010.11': '10',
                'NC_000011.10': '11',
                'NC_000012.12': '12',
                'NC_000013.11': '13',
                'NC_000014.9': '14',
                'NC_000015.10': '15',
                'NC_000016.10': '16',
                'NC_000017.11': '17',
                'NC_000018.10': '18',
                'NC_000019.10': '19',
                'NC_000020.11': '20',
                'NC_000021.9': '21',
                'NC_000022.11': '22',
                'NC_000023.11': 'X',
                'NC_000024.10': 'Y',
                'NC_012920.1': 'MT',
                'chr1': '1',
                'chr2': '2',
                'chr3': '3',
                'chr4': '4',
                'chr5': '5',
                'chr6': '6',
                'chr7': '7',
                'chr8': '8',
                'chr9': '9',
                'chr10': '10',
                'chr11': '11',
                'chr12': '12',
                'chr13': '13',
                'chr14': '14',
                'chr15': '15',
                'chr16': '16',
                'chr17': '17',
                'chr18': '18',
                'chr19': '19',
                'chr20': '20',
                'chr21': '21',
                'chr22': '22',
                'chrX': 'X',
                'chrY': 'Y',
                'chrM': 'MT',
                'MT': 'MT',
                'M': 'MT'
            },
            'GRCh37': {
                'CM000663.1': '1',
                'CM000664.1': '2',
                'CM000665.1': '3',
                'CM000666.1': '4',
                'CM000667.1': '5',
                'CM000668.1': '6',
                'CM000669.1': '7',
                'CM000670.1': '8',
                'CM000671.1': '9',
                'CM000672.1': '10',
                'CM000673.1': '11',
                'CM000674.1': '12',
                'CM000675.1': '13',
                'CM000676.1': '14',
                'CM000677.1': '15',
                'CM000678.1': '16',
                'CM000679.1': '17',
                'CM000680.1': '18',
                'CM000681.1': '19',
                'CM000682.1': '20',
                'CM000683.1': '21',
                'CM000684.1': '22',
                'CM000685.1': 'X',
                'CM000686.1': 'Y',
                'J01415.2': 'MT',
                'NC_000001.10': '1',
                'NC_000002.11': '2',
                'NC_000003.11': '3',
                'NC_000004.11': '4',
                'NC_000005.9': '5',
                'NC_000006.11': '6',
                'NC_000007.13': '7',
                'NC_000008.10': '8',
                'NC_000009.11': '9',
                'NC_000010.10': '10',
                'NC_000011.9': '11',
                'NC_000012.11': '12',
                'NC_000013.10': '13',
                'NC_000014.8': '14',
                'NC_000015.9': '15',
                'NC_000016.9': '16',
                'NC_000017.10': '17',
                'NC_000018.9': '18',
                'NC_000019.9': '19',
                'NC_000020.10': '20',
                'NC_000021.8': '21',
                'NC_000022.10': '22',
                'NC_000023.10': 'X',
                'NC_000024.9': 'Y',
                'NC_012920.1': 'MT',
                'chr1': '1',
                'chr2': '2',
                'chr3': '3',
                'chr4': '4',
                'chr5': '5',
                'chr6': '6',
                'chr7': '7',
                'chr8': '8',
                'chr9': '9',
                'chr10': '10',
                'chr11': '11',
                'chr12': '12',
                'chr13': '13',
                'chr14': '14',
                'chr15': '15',
                'chr16': '16',
                'chr17': '17',
                'chr18': '18',
                'chr19': '19',
                'chr20': '20',
                'chr21': '21',
                'chr22': '22',
                'chrX': 'X',
                'chrY': 'Y',
                'chrM': 'MT',
                'MT': 'MT',
                'M': 'MT'
            }
        }

        if build in chroms and chrom in chroms[build]:
            return chroms[build][chrom]

        if build is None:
            for any_build in chroms:
                if chrom in chroms[any_build]:
                    return chroms[any_build][chrom]

        return chrom

    @staticmethod
    def from_ebi(chrom, build, style):
        ''' Convert major chromosomes from EBI to any style '''

        # UCSC is the same for GRCh37 and 38
        if style == refstyle.UCSC:
            if chrom in (
                    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
                    'X', 'Y'
            ):
                return f'chr{chrom}'

            if chrom == 'MT':
                return 'chrM'

            return chrom

        # Everything else is a lookup
        chroms = {
            'GRCh38': {
                refstyle.GenBank: {
                    '1': 'CM000663.2',
                    '2': 'CM000664.2',
                    '3': 'CM000665.2',
                    '4': 'CM000666.2',
                    '5': 'CM000667.2',
                    '6': 'CM000668.2',
                    '7': 'CM000669.2',
                    '8': 'CM000670.2',
                    '9': 'CM000671.2',
                    '10': 'CM000672.2',
                    '11': 'CM000673.2',
                    '12': 'CM000674.2',
                    '13': 'CM000675.2',
                    '14': 'CM000676.2',
                    '15': 'CM000677.2',
                    '16': 'CM000678.2',
                    '17': 'CM000679.2',
                    '18': 'CM000680.2',
                    '19': 'CM000681.2',
                    '20': 'CM000682.2',
                    '21': 'CM000683.2',
                    '22': 'CM000684.2',
                    'X': 'CM000685.2',
                    'Y': 'CM000686.2',
                    'MT': 'J01415.2',
                },
                refstyle.RefSeq: {
                    '1': 'NC_000001.11',
                    '2': 'NC_000002.12',
                    '3': 'NC_000003.12',
                    '4': 'NC_000004.12',
                    '5': 'NC_000005.10',
                    '6': 'NC_000006.12',
                    '7': 'NC_000007.14',
                    '8': 'NC_000008.11',
                    '9': 'NC_000009.12',
                    '10': 'NC_000010.11',
                    '11': 'NC_000011.10',
                    '12': 'NC_000012.12',
                    '13': 'NC_000013.11',
                    '14': 'NC_000014.9',
                    '15': 'NC_000015.10',
                    '16': 'NC_000016.10',
                    '17': 'NC_000017.11',
                    '18': 'NC_000018.10',
                    '19': 'NC_000019.10',
                    '20': 'NC_000020.11',
                    '21': 'NC_000021.9',
                    '22': 'NC_000022.11',
                    'X': 'NC_000023.11',
                    'Y': 'NC_000024.10',
                    'MT': 'NC_012920.1',
                },
            },
            'GRCh37': {
                refstyle.GenBank: {
                    '1': 'CM000663.1',
                    '2': 'CM000664.1',
                    '3': 'CM000665.1',
                    '4': 'CM000666.1',
                    '5': 'CM000667.1',
                    '6': 'CM000668.1',
                    '7': 'CM000669.1',
                    '8': 'CM000670.1',
                    '9': 'CM000671.1',
                    '10': 'CM000672.1',
                    '11': 'CM000673.1',
                    '12': 'CM000674.1',
                    '13': 'CM000675.1',
                    '14': 'CM000676.1',
                    '15': 'CM000677.1',
                    '16': 'CM000678.1',
                    '17': 'CM000679.1',
                    '18': 'CM000680.1',
                    '19': 'CM000681.1',
                    '20': 'CM000682.1',
                    '21': 'CM000683.1',
                    '22': 'CM000684.1',
                    'X': 'CM000685.1',
                    'Y': 'CM000686.1',
                    'MT': 'J01415.2',
                },
                refstyle.RefSeq: {
                    '1': 'NC_000001.10',
                    '2': 'NC_000002.11',
                    '3': 'NC_000003.11',
                    '4': 'NC_000004.11',
                    '5': 'NC_000005.9',
                    '6': 'NC_000006.11',
                    '7': 'NC_000007.13',
                    '8': 'NC_000008.10',
                    '9': 'NC_000009.11',
                    '10': 'NC_000010.10',
                    '11': 'NC_000011.9',
                    '12': 'NC_000012.11',
                    '13': 'NC_000013.10',
                    '14': 'NC_000014.8',
                    '15': 'NC_000015.9',
                    '16': 'NC_000016.9',
                    '17': 'NC_000017.10',
                    '18': 'NC_000018.9',
                    '19': 'NC_000019.9',
                    '20': 'NC_000020.10',
                    '21': 'NC_000021.8',
                    '22': 'NC_000022.10',
                    'X': 'NC_000023.10',
                    'Y': 'NC_000024.9',
                    'MT': 'NC_012920.1',
                },
            }
        }

        if build in chroms and style in chroms[build] and chrom in chroms[build][style]:
            return chroms[build][style][chrom]

        return chrom

    def to_ucsc(self, chrom, build):
        ''' Convert major chromosomes to UCSC style '''
        return self.from_ebi(self.to_ebi(chrom, build), build, refstyle.UCSC)

    def to_genbank(self, chrom, build):
        ''' Convert major chromosomes to NCBI GenBank style '''
        return self.from_ebi(self.to_ebi(chrom, build), build, refstyle.GenBank)

    def to_refseq(self, chrom, build):
        ''' Convert major chromosomes to NCBI RefSeq style '''
        return self.from_ebi(self.to_ebi(chrom, build), build, refstyle.RefSeq)

    def to_native(self, chrom, build, style):
        ''' Convert major chromosomes to the specified style '''
        if style == refstyle.EBI:
            return self.to_ebi(chrom, build)

        if style == refstyle.GenBank:
            return self.to_genbank(chrom, build)

        if style == refstyle.RefSeq:
            return self.to_refseq(chrom, build)

        if style == refstyle.UCSC:
            return self.to_ucsc(chrom, build)

        return chrom

def parse_args(clargs):
    ''' biograph vdb import args '''
    parser = argparse.ArgumentParser(
        description=main.__doc__
    )

    parser.add_argument("input", type=str, default="/dev/stdin", nargs='?',
                        help="Input filename or refdir (%(default)s)")
    parser.add_argument("-c", "--common", action='store_true',
                        help="Print the common name if known, otherwise print the hash")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="Print the hash and the common name")
    parser.add_argument("-f", "--full", action='store_true',
                        help="Print all available information")
    parser.add_argument("-l", "--list", action='store_true',
                        help="List all known hashes and names")
    parser.add_argument("--debug", action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args(clargs)

    if sys.stdin.isatty() and args.input == '/dev/stdin' and not args.list:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return args

def main(clargs):
    ''' Identify the reference in a VCF, FASTA, SAM, or BioGraph refdir '''
    args = parse_args(clargs)

    if args.list:
        for k in refhashes:
            print("refhash:", k)
            for sub in refhashes[k]:
                print(f"{sub:>7}:", refhashes[k][sub])
            print()
        sys.exit(0)

    rh = refhash(args.input)

    if args.full:
        print(rh.full_info())
    elif args.verbose:
        print(rh.verbose_name())
    elif args.common:
        print(rh.common_name())
    else:
        print(rh.digest())

    if args.debug:
        for contig in rh.contigs:
            print(contig, file=sys.stderr)

if __name__ == '__main__':
    main(sys.argv[1:])
