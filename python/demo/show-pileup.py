#!/usr/bin/python3

import biograph
import argparse
from random import randint

parser = argparse.ArgumentParser(description="Show pileup for a sequence")
parser.add_argument("biograph", help="Biograph to use")
parser.add_argument(
    "sequence", help="Sequence to show pileup for")
parser.add_argument("--show-mates", help="Show mates of each read encountered",
                    action="store_true")
parser.add_argument(
    "--show-seqset", help="Show overlapping seqset entries whether or not they have an associated read in any of the given readmaps",
                    action="store_true")
readmap_group = parser.add_mutually_exclusive_group()
readmap_group.add_argument(
    "--all-readmaps", help="Open all readmaps", action="store_true")
readmap_group.add_argument(
    "--readmap", dest="readmaps", nargs="*", help="Open the given readmap (may be specified multiple times)")

args = parser.parse_args()

bg = biograph.BioGraph(args.biograph)

readmaps = dict()
if args.all_readmaps:
    for accession_id, _ in bg.metadata.samples.iteritems():
        readmaps[accession_id] = bg.open_readmap(accession_id)
elif args.readmaps:
    for accession_id in args.readmaps:
        readmaps[accession_id] = bg.open_readmap(accession_id)
else:
    readmaps["default"] = bg.open_readmap("")

SEQUENCE_HEADER = "Sequence"
accession_id_len = max([len(name)
                       for name in readmaps.iterkeys()] + [len(SEQUENCE_HEADER)])

print("Searching through %d readmaps: %s" %
      (len(readmaps), ",".join(readmaps.keys())))

# Display the main sequence up top
print("%-*s %s" % (accession_id_len, SEQUENCE_HEADER, args.sequence))

# Start from the right to find all sequences present
entry = bg.seqset.empty_entry()
offset = len(args.sequence)
for base in reversed(args.sequence):
    new_entry = entry.push_front_drop(base)
    new_offset = offset - 1

    if args.show_seqset and len(new_entry) < len(entry):
        seqset_len = len(entry)
        seqset_seq = args.sequence[offset:offset + seqset_len]
        print("%-*s %*s%s" %
              (accession_id_len, "(seqset)", offset, "", seqset_seq))

    entry = new_entry
    offset = new_offset
    # See what readmaps have this sequence present
    for accession_id, rm in readmaps.iteritems():
        for read in rm.get_prefix_reads(entry):
            read_len = len(read)
            read_seq = args.sequence[offset:offset + read_len]
            print("%-*s %*s%s" %
                  (accession_id_len, accession_id, offset, "", read_seq))
