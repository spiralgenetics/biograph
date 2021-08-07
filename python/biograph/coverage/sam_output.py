"""
Generates SAM/BAM output from biograph ReadCoverages
"""

import logging
import pysam
import biograph.variants as bgexvar

class SamMapper:
    """SamMapper provides a facility for aligned reads to be written to a sam output file.

A single instance of SamMapper supports a single supercontig at a time."""

    def __init__(self, readmap, output_queue, chrom, emit_all_mates=False):
        self.readmap = readmap
        self.output_queue = output_queue
        self.chrom = chrom
        self.queued_lines = []
        self.emit_all_mates = emit_all_mates
        # For emit_all_mates in aligned mode
        self.read1_output = bgexvar.BigReadIdSet()
        self.read2_output = bgexvar.BigReadIdSet()
        # For unaligned mode:
        self.all_read_ids = bgexvar.BigReadIdSet()

    def aligned_to_sam(self, read_id, aligned):
        "Converts an AlignedRead to a line in SAM format"
        flag = 0
        read = self.readmap.get_read_by_id(read_id)
        if read.is_original_orientation():
            fwd_read = read
            fwd_read_id = read_id
        else:
            flag |= pysam.FREVERSE
            fwd_read = read.get_rev_comp()
            fwd_read_id = fwd_read.get_read_id()
        if fwd_read.has_mate():
            flag |= pysam.FPAIRED
            canon_read_id = min(fwd_read_id, fwd_read.get_mate().get_read_id())
            if canon_read_id == fwd_read_id:
                flag |= pysam.FREAD1
                if self.emit_all_mates:
                    self.read1_output.add(canon_read_id)
            else:
                flag |= pysam.FREAD2
                if self.emit_all_mates:
                    self.read2_output.add(canon_read_id)
        else:
            flag |= pysam.FREAD1
            canon_read_id = fwd_read_id
        line = b"%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            canon_read_id, # QNAME
            flag, # FLAG
            self.chrom.encode(), # RNAME
            aligned.left_offset + 1, # POS
            b"255", # MAPQ
            aligned.cigar.encode(), # CIGAR
            b"*", # RNEXT
            b"0", # PNEXT
            b"0", # TLEN
            str(aligned.seq).encode(), # SEQ
            b"*") # QUAL
        return line

    @staticmethod
    def unaligned_read_to_sam(read_id, seq, flag):
        "Generates a SAM format for a read which hasn't been aligned."
        line = b"%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            read_id, # QNAME
            flag, # FLAG
            b"*", # RNAME
            b"0", # POS
            b"0", # MAPQnnnn
            b"*", # CIGAR
            b"*", # RNEXT
            b"0", # PNEXT
            b"0", # TLEN
            str(seq).encode(), # SEQ
            b"*") # QUAL
        return line

    def get_sequences(self, canon_read_ids, mate):
        "Lookups up sequences, caching per seqset entry"
        last_seqset_entry = None
        last_seq = None
        if mate:
            rds = ((canon_read_id, self.readmap.get_read_by_id(canon_read_id))
                   for canon_read_id in canon_read_ids)
            rds = ((rd.get_read_id(), canon_read_id, rd)
                   for canon_read_id, rd in rds)
            # resort by mate_id
            rds = sorted(rds)
            rds = ((canon_read_id, seq) for _, canon_read_id, seq in rds)
        else:
            rds = ((canon_read_id, self.readmap.get_read_by_id(canon_read_id))
                   for canon_read_id in canon_read_ids)
        for canon_read_id, rd in rds:
            seqset_entry = rd.get_seqset_entry()
            if seqset_entry != last_seqset_entry:
                last_seqset_entry = seqset_entry
                last_seq = last_seqset_entry.sequence()
            yield canon_read_id, last_seq

    def emit_mates(self):
        "For reads which we haven't emitted the mate, emit the mates."
        read1_output = self.read1_output.to_read_id_set()
        self.read1_output = bgexvar.BigReadIdSet()
        read2_output = self.read2_output.to_read_id_set()
        self.read2_output = bgexvar.BigReadIdSet()
        read1_needed = read2_output - read1_output
        read2_needed = read1_output - read2_output

        for read_id, seq in self.get_sequences(read1_needed, mate=False):
            sam_line = self.unaligned_read_to_sam(read_id, seq,
                                                  pysam.FPAIRED | pysam.FREAD1 | pysam.FUNMAP)
            self.queue_sam(sam_line)
        for read_id, seq in self.get_sequences(read2_needed, mate=True):
            sam_line = self.unaligned_read_to_sam(read_id, seq,
                                                  pysam.FPAIRED | pysam.FREAD2 | pysam.FUNMAP)
            self.queue_sam(sam_line)

    def on_aligned(self, read_ids, aligned):
        "Writes the given AlignedRead to the output SAM stream"
        for read_id in read_ids:
            sam_line = self.aligned_to_sam(read_id, aligned)
            self.queue_sam(sam_line)

    def on_unaligned(self, read_ids):
        "Saves the given read IDs to be written to output."
        self.all_read_ids |= read_ids

    def emit_unaligned(self):
        "Writes out unaligned entries for all reads in all_read_ids."
        # Catagorize the reads
        read1s = bgexvar.BigReadIdSet()
        read2s = bgexvar.BigReadIdSet()
        unpaired = bgexvar.BigReadIdSet()

        logging.info(f"emit_unaligned {self.chrom}: collating read ids")
        for read_id in self.all_read_ids:
            read = self.readmap.get_read_by_id(read_id)
            if not read.is_original_orientation():
                read = read.get_rev_comp()
                read_id = read.get_read_id()
            if read.has_mate():
                mate_id = read.get_mate().get_read_id()
                canon_read_id = min(read_id, mate_id)
                if read_id == canon_read_id:
                    read1s.add(read_id)
                    if self.emit_all_mates:
                        read2s.add(mate_id)
                else:
                    read2s.add(read_id)
                    if self.emit_all_mates:
                        read1s.add(mate_id)
            else:
                unpaired.add(read_id)
        self.all_read_ids = bgexvar.BigReadIdSet()

        logging.info(f"flush_unaligned {self.chrom}: outputting read1s")
        for read_id, seq in self.get_sequences(read1s, mate=False):
            sam_line = self.unaligned_read_to_sam(read_id, seq,
                                                  pysam.FPAIRED | pysam.FREAD1 | pysam.FUNMAP)
            self.queue_sam(sam_line)

        logging.info(f"flush_unaligned {self.chrom}: outputting read2s")
        for read_id, seq in self.get_sequences(read2s, mate=False):
            canon_read_id = self.readmap.get_read_by_id(read_id).get_mate().get_read_id()
            sam_line = self.unaligned_read_to_sam(canon_read_id, seq,
                                                  pysam.FPAIRED | pysam.FREAD2 | pysam.FUNMAP)
            self.queue_sam(sam_line)

        logging.info(f"flush_unaligned {self.chrom}: outputting unpaired")
        for read_id, seq in self.get_sequences(unpaired, mate=False):
            sam_line = self.unaligned_read_to_sam(read_id, seq,
                                                  pysam.FREAD1 | pysam.FUNMAP)
            self.queue_sam(sam_line)
        logging.info(f"flush_unaligned {self.chrom}: done")

    def queue_sam(self, sam_line):
        "Enqueues the given sam line to be output"
        self.queued_lines.append(sam_line)
        if len(self.queued_lines) >= 100:
            self.flush_queued()

    def flush_queued(self):
        "Flush queued output lines to the SAM writer thread"
        self.output_queue.put(b"".join(self.queued_lines))
        self.queued_lines.clear()

class SamOutput:
    "SamOutput supports aligning reads and writing their output to a SAM file"

    def __init__(self, readmap, *, refskip_anchor=False, emit_all_mates=False):
        self.reference_id_by_name = {}
        self.header = []
        self.readmap = readmap
        self.refskip_anchor = refskip_anchor
        self.emit_all_mates = emit_all_mates

    def add_reference(self, ref):
        "Adds all supercontigs from the given biograph reference to SAM header"
        for scaffold_name, scaffold_len in ref.scaffold_lens.items():
            self.add_scaffold(scaffold_name, scaffold_len)

    def get_header(self):
        "Returns the SAM header in a format suitable for pysam.AlignmentFile"
        return self.header

    def add_scaffold(self, scaffold_name, scaffold_len):
        "Adds the given supercontig to the header"
        if scaffold_name in self.reference_id_by_name:
            # Already have this reference
            return
        reference_id = len(self.reference_id_by_name)
        self.reference_id_by_name[scaffold_name] = reference_id
        self.header.append(b"@SQ\tSN:%s\tLN:%d\n" % (
            scaffold_name.encode(), scaffold_len))

    def output_unaligned(self, output_queue, scaffold_name, entries):
        "Outputs unaligned reads present in any of the given entries"
        mapper = SamMapper(self.readmap, output_queue, scaffold_name,
                           emit_all_mates=self.emit_all_mates)
        for a in entries:
            mapper.on_unaligned(a.read_coverage.all_read_ids())
            yield a
        mapper.emit_unaligned()
        mapper.flush_queued()

    def output_aligned(self, output_queue, scaffold_name, entries):
        """Aligns the reads in the given assemblies and outputs them to the given sam stream.

Returns a stream of the input assemblies which must be consumed."""
        mapper = SamMapper(self.readmap, output_queue, scaffold_name,
                           emit_all_mates=self.emit_all_mates)
        yield from bgexvar.align_reads(entries,
                                       on_aligned=mapper.on_aligned,
                                       refskip_anchor=self.refskip_anchor)
        if self.emit_all_mates:
            mapper.emit_mates()
        mapper.flush_queued()
