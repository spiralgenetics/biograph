#!/usr/bin/env python3
"""
    vcf.py

    Utility functions for manipulating VCF files.

    http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
    http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
"""
# NOTE: This is a very incomplete implementation!

DEFAULT_VCF_VERSION = '4.1'


def vcf_format(version=DEFAULT_VCF_VERSION):
    """
        The canonical definition of the VCF format.

        'name' is the name of the field.
        'type' is a function that should be applied to convert it to the expected type when parsing.
    """

    if (version not in ['4.0', '4.1']):
        raise Exception('Unsupported VCF version: %s' % version)

    # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample

    return [
        {
            'name': 'CHROM',
            'type': str
        },
        {
            'name': 'POS',
            'type': int
        },
        {
            'name': 'ID',
            'type': str
        },
        {
            'name': 'REF',
            'type': str
        },
        {
            'name': 'ALT',
            'type': str
        },
        {
            'name': 'QUAL',
            'type': float
        },
        {
            'name': 'FILTER',
            'type': str
        },
        {
            'name': 'INFO',
            'type': str
        },
        {
            'name': 'FORMAT',
            'type': str
        },
        {
            'name': 'Sample',
            'type': str
        },
    ]


def vcf_column_header(version=DEFAULT_VCF_VERSION):
    """
        The VCF column names, in order.
    """
    return [field['name'] for field in vcf_format(version)]


def variant_size(variant):
    """ Can only call this function on a vcf variant dict """
    return abs(len(variant['ALT']) - len(variant['REF']))


def variant_type(variant):
    """ does what it says on the tin """
    return call_vcf_variant(variant)


def classify_variant(variant):
    """ does what it says on the tin """
    var_size = abs(len(variant['ALT']) - len(variant['REF']))

    alt_len = len(variant['ALT'])
    ref_len = len(variant['REF'])

    # SNPs are size zero
    if var_size == 0:
        return 'snp', var_size

    # An ALT smaller than REF is a delete
    if alt_len < ref_len:
        return 'del', var_size

    # Anything else is an insert
    return 'ins', var_size


def parse_vcf_line(line, version=DEFAULT_VCF_VERSION):
    """
        Parse one line from an VCF file. Return the variant type and all columns.

        Raises an exception on parse failure.
    """
    if (version not in ['4.0', '4.1']):
        raise Exception('Unsupported VCF version: %s' % version)

    variant = {}
    props = vcf_column_header()

    entries = line.split('\t')

    if len(props) != len(entries):
        raise Exception('Expected %d fields, got %d\n%s' % (len(props), len(entries), line))

    for i, field in enumerate(vcf_format()):
        variant[field['name']] = field['type'](entries[i]) if entries[i] != '.' else field['type']()

    return variant


def call_vcf_variant(variant):
    """
        Given a variant object, return the variant type.
    """
    alt_len = len(variant['ALT'])
    ref_len = len(variant['REF'])

    var_size = variant_size(variant)

    # Unfortunately, GATK doesn't use variant calling, even though it claims v4.1 support. Sigh.

    # SNPs are size zero
    if var_size == 0:
        return 'snp'

    # An ALT smaller than REF is a delete
    if alt_len < ref_len:
        return 'del'

    # Anything else is an insert
    return 'ins'

    # Inversions and repeats need to be detected with a post-processing pass, since
    # multiple lines must be considered when making that determination

    # TODO: VCF v4.1 variant calling support (#ALT=...) and imprecise calling (<DEL>)
    # See http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41

    # TODO: What about mixed VCF records (ALT=ACT,ACG,ATA...)?

# def make_vcf_line(variant, version=DEFAULT_VCF_VERSION):
#     """
#         Given a variant object, generate an equivalent VCF line.

#         Inverse of parse_vcf_line()
#     """

#     raise Exception('Not implemented yet.')

    # props = svs_header()

    # if len(variant.keys()) != (len(props) + 2):
    #     raise Exception('Wrong number of columns (%d should be %d)' % (len(variant.keys()), len(props)))

    # variants = []
    # for header in props:
    #     variants.append(str(variant[header]))

    # return '\t'.join(variants)

# def count_vcf_variations(file_name):
#     """
#         Count all of the variations in a VCF file.

#         Returns an object indicating the number of records found for each variation type.
#         Raises an exception if an invalid vcf file is detected.
#     """

#     raise Exception('Not implemented yet.')

    # count = {
    #     'total': 0,
    #     'snp': 0,
    #     'del': 0,
    #     'ins': 0,
    #     'trn': 0,
    #     'rpt': 0
    # }
    # with open(file_name, 'r') as f:
    #     for line in f:
    # if line[0] == '#':
    #             continue
    #         variant = parse_svs_line(line)
    #         count[variant['var_type']] = count[variant['var_type']] + 1
    #         count['total'] = count['total'] + 1

    # return count


def enumerate_vcf_variations(file_name):
    """
        Return a list of all unique variants in file_name.
    """
    if not file_name:
        raise Exception('enumerate_variations: file_name required.')

    variants = []

    # TODO: Parse VCF header and grab version, etc.

    with open(file_name, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            variant = parse_vcf_line(line)

            if variant not in variants:
                variants.append(variant)

    return variants

# def compare_vcf_variations(file_one, file_two, position_window=0):
#     """
#         Return a list of variation differences between two files.

#         Fields listed in do_not_compare[] are ignored on comparison, but are
#         returned in the variant list.
#     """

#     raise Exception('Not implemented yet.')


def complement(sequence, reverse=False):
    """
        Compute the complement of a DNA sequence.

        If reverse is True, reverse it too.
    """
    flip = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    revcomp = []
    for i in list(sequence):
        revcomp.append(flip[i])

    if reverse:
        return ''.join(revcomp[::-1])
    return ''.join(revcomp)
