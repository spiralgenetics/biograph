"""
    fasta.py

    Utility functions for manipulating fasta files.

    NOTE: This is largely incomplete; only .fai support is currently implemented.
"""
from python.functest.utils.fileops import opener

def fai_format():
    """
        The canonical definition of the fai format.

        'name' is the name of the field.
        'type' is a function that should be applied to convert it to the expected type when parsing.
    """

    # http://samtools.sourceforge.net/samtools.shtml
    # http://www.biostars.org/p/11523/#11524

    return [
        {
            'name': 'name',
            'type': str
        },
        {
            'name': 'length',
            'type': int
        },
        {
            'name': 'offset',
            'type': int
        },
        {
            'name': 'bases_per_line',
            'type': int
        },
        {
            'name': 'bytes_per_line',
            'type': int
        },
    ]

def fai_column_header():
    """
        The .fai column names, in order.
    """
    return [field['name'] for field in fai_format()]

def fasta_index(fai_file):
    """
        Return a lookup table for a fasta index.
    """
    if not fai_file:
        raise Exception('fai file required.')

    faindex = {}
    with opener(fai_file, 'r') as fafile:
        for line in fafile:
            entry = parse_fai_line(line)
            faindex[entry['name']] = entry

    return faindex

def parse_fai_line(line):
    """
        Parse a single fasta index line. Return a dict containing the index entry.
    """
    entry = {}
    props = fai_column_header()
    entries = line.split('\t')

    if len(props) != len(entries):
        raise Exception('Expected %d fields, got %d\n%s' % (len(props), len(entries), line))

    for i, field in enumerate(fai_format()):
        entry[field['name']] = field['type'](entries[i])

    if ':' in entry['name']:
        entry['chromosome'], entry['contig'] = entry['name'].split(':')
        entry['contig'] = int(entry['contig'])

    return entry
