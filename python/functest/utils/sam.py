#!/usr/bin/env python3
"""
    sam.py

    Utility functions for parsing SAM files.
"""


def sam_format():
    """
        The canonical definition of the SAM format.

        http://samtools.github.io/hts-specs/SAMv1.pdf

        FLAG definitions:
            0x1 template having multiple segments in sequencing
            0x2 each segment properly aligned according to the aligner
            0x4 segment unmapped
            0x8 next segment in the template unmapped
            0x10 SEQ being reverse complemented
            0x20 SEQ of the next segment in the template being reversed
            0x40 the first segment in the template
            0x80 the last segment in the template
            0x100 secondary alignment
            0x200 not passing quality controls
            0x400 PCR or optical duplicate
            0x800 supplementary alignment

        CIGAR string:
            Op BAM Description
            M 0 alignment match (can be a sequence match or mismatch)
            I 1 insertion to the reference
            D 2 deletion from the reference
            N 3 skipped region from the reference
            S 4 soft clipping (clipped sequences present in SEQ)
            H 5 hard clipping (clipped sequences NOT present in SEQ)
            P 6 padding (silent deletion from padded reference)
            = 7 sequence match
            X 8 sequence mismatch

        'name' is the name of the field.
        'type' is a function that should be applied to convert it to the expected type when parsing.
    """
    all_fields = [
        {
            # Query template NAME
            'name': 'QNAME',
            'type': str
        },
        {
            # bitwise FLAG
            'name': 'FLAG',
            'type': int
        },
        {
            # Reference sequence NAME
            'name': 'RNAME',
            'type': str
        },
        {
            # 1-based leftmost mapping POSition
            'name': 'POS',
            'type': int
        },
        {
            # MAPping Quality
            'name': 'MAPQ',
            'type': int
        },
        {
            # CIGAR string
            'name': 'CIGAR',
            'type': str
        },
        {
            # Ref. name of the mate/next read
            'name': 'RNEXT',
            'type': str
        },
        {
            # Position of the mate/next read
            'name': 'PNEXT',
            'type': int
        },
        {
            # observed Template LENgth
            'name': 'TLEN',
            'type': int
        },
        {
            # segment SEQuence
            'name': 'SEQ',
            'type': str
        },
        {
            # ASCII of Phred-scaled base QUALity+33
            'name': 'QUAL',
            'type': str
        }
    ]

    return all_fields


def parse_sam_line(samline):
    """
        Parse one line from a SAM file. Return a single alignment.

        Raises an exception on parse failure.
    """
    props = [field['name'] for field in sam_format()]

    alignment = dict(zip(props, samline.rstrip().split('\t')))

    if len(props) != len(alignment):
        raise Exception('Expected %d fields, got %d\n%s' % (len(props), len(alignment), samline))

    for field in sam_format():
        alignment[field['name']] = field['type'](alignment[field['name']])

    return alignment


if __name__ == '__main__':
    import fileinput

    for line in fileinput.input():
        # Skip tags
        if line[0] == '@':
            continue

        assembly = parse_sam_line(line)

        print(assembly)
