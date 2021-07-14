"""Converts variant data to parquet

Most parsing is done with 'bytes' instead of 'string' for efficiency reasons.

"""

import collections
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile

import logging
from logging import debug

from pathlib import Path

import pandas as pd

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.ipc as ipc

from biograph.tools.refhash import refhash

class ConverterBase:
    """Base class for converting variant file formats to parquet."""

    # Number of rows per chunk when reading input into chunks.
    CHUNK_ROWS = 500_000

    # Current instance of ConverterBase, since multiprocessing.pool
    instance = None

    # Default span size to use for additional "span" column for efficient range joins.
    # Large variants will use approximately (varend - pos) / SPAN_BASES
    # additional rows to associate the variants with multiple regions.
    DEFAULT_SPAN_BASES = 1000

    def __init__(self, build, tmpdir, outfile, schema, verbose=False,
                 span_bases=DEFAULT_SPAN_BASES, nthreads=1):
        self.build = build
        self.tmpdir = tmpdir
        self.schema = schema
        self.chunks = []
        self.nrows = 0
        self.outfile = outfile
        self.verbose = verbose
        self.span_bases = span_bases
        self.nthreads = nthreads
        self.header_lines = []

        logging.basicConfig(
            stream=sys.stderr,
            level=logging.DEBUG if self.verbose else logging.WARNING,
            format='%(message)s'
        )

        first_fields = [field.name for field in list(schema)[0:5]]
        if first_fields != ["spans", "reflen", "chrom", "pos", "varend"]:
            raise RuntimeError(f"First four schema fields must be spans, reflen, chrom, pos, and varend; got {first_fields}")

    def parse_chunk(self, fh):  # pylint: disable=no-self-use
        """
        Subclasses must implement parse_chunk.

        This method should read lines from fh and yield a tuple with two
        elements.  The first element is a tuple containing any additional keys
        that should be used to categorize into row groups.  The second element
        is a tuple that should match the schema starting at the 'chrom' field;
        everything before that is added automatically.
        """
        raise RuntimeError("Subclasses must implement parse_chunk")

    def make_chunks_from(self, in_fh):
        "Reads variant data from in_fh, splits into chunks, and saves the chunks to the temporary directory."

        debug("Making chunks")

        chunk_dir = Path(self.tmpdir) / "chunks"
        chunk_dir.mkdir()

        line = in_fh.readline()
        while line.startswith(b"#"):
            stripped = line.strip().strip(b"#")
            if not stripped:
                # Some files have random blank comment lines throughout the file;
                # don't save these.
                continue
            self.header_lines.append(line.decode(encoding="UTF8").rstrip())
            line = in_fh.readline()

        debug(f"{len(self.header_lines)} header lines saved")

        # ~128MB chunks seems optimal, but the value isn't critical.
        # Too big + too many threads == out of memory later.
        # -C breaks on \n boundaries
        psplit = subprocess.Popen(
            [
                "/usr/bin/split",
                "-C", "128000000",
                "-d", "-",
                f"{chunk_dir}/x"
            ],
            stdin=subprocess.PIPE
        )
        # save the first non-header line
        psplit.stdin.write(line)
        shutil.copyfileobj(in_fh, psplit.stdin)
        psplit.stdin.close()
        psplit.wait()

        self.chunks = [(str(f.absolute()), i) for i, f in enumerate(chunk_dir.glob('*'))]
        debug(f"Generated {len(self.chunks)} input chunks")

    def write_table_chunk(self, rows):
        "Writes out the given rows to a temporary file in apache arrow format"
        df = pd.DataFrame(
            columns=[f.name for f in self.schema],
            data=rows,
            dtype=object).sort_values("pos")
        batch = pa.RecordBatch.from_pandas(df, schema=self.schema)
        with tempfile.NamedTemporaryFile(dir=self.tmpdir, delete=False) as out_fh:
            writer = ipc.new_file(out_fh, self.schema)
            writer.write(batch)
            writer.close()
            return out_fh.name

    @staticmethod
    def process_chunk_in_instance(spec):
        "Helper method so that multiprocessing can call process_chunk"
        return ConverterBase.instance.process_chunk(spec)

    def process_chunk(self, spec):
        "Converts a given chunk specification from variant format to apache arrow format, grouped by chrom"

        filename, chunk_num = spec
        debug(f"Processing chunk #{chunk_num}: {filename}")
        out_chunks = collections.defaultdict(list)
        out_table_chunks = []
        last_chrom = None
        ebi_chrom = None
        span_bases = self.span_bases
        with open(filename, "rb") as fh:
            for _, line in self.parse_chunk(fh):
                chrom, start, end = line[0:3]
                if chrom != last_chrom:
                    ebi_chrom = refhash.to_ebi(chrom.decode(encoding='UTF8'), build=self.build)
                    last_chrom = chrom
                reflen = end - start
                line = (range(start // span_bases, (end // span_bases) + 1), # spans
                        reflen, ebi_chrom, *line[1:])

                out_chunk = out_chunks[ebi_chrom]
                out_chunk.append(line)
                if len(out_chunk) >= self.CHUNK_ROWS:
                    fn = self.write_table_chunk(out_chunk)
                    out_table_chunks.append(fn)
                    out_chunks[ebi_chrom].clear()

            for ebi_chrom, out_chunk in out_chunks.items():
                fn = self.write_table_chunk(out_chunk)
                out_table_chunks.append(fn)
                out_chunks[ebi_chrom].clear()

        os.unlink(filename)
        return out_table_chunks

    def generate(self):
        ''' Converts chunks saved from make_chunks_from to apache arrow format, and outputs a parquet file. '''
        debug(f"Generating chunks on {self.nthreads} threads")

        all_chunk_tables = []
        if self.nthreads > 1:
            ConverterBase.instance = self
            with multiprocessing.Pool(self.nthreads) as p:
                out_chunks = p.imap_unordered(ConverterBase.process_chunk_in_instance, self.chunks)
                p.close()
                for chunk_tables in out_chunks:
                    all_chunk_tables.extend(chunk_tables)
                p.join()
        else:
            for chunk_tables in map(self.process_chunk, self.chunks):
                all_chunk_tables.extend(chunk_tables)

        debug(f"Generated {len(all_chunk_tables)} grouped chunks from {len(self.chunks)} input chunks")

        with multiprocessing.Pool(self.nthreads) as p:
            p.imap_unordered(self.save_parquet, all_chunk_tables)
            p.close()
            p.join()

    def save_parquet(self, fn):
        "Writes tables grouped by chunks to the output parquet dataset"
        debug(f"Writing parquet from {fn}")
        chunk_batches = []
        reader = ipc.open_file(pa.OSFile(fn))
        chunk_batches.append(reader.get_batch(0))
        os.unlink(fn)
        pq.write_to_dataset(
            pa.Table.from_batches(chunk_batches),
            root_path=self.outfile
        )
        reader.close()

def parse_vcf_info(info):
    "Parses a VCF INFO field into a list of key, value pairs"
    if info == b".":
        return []

    fields = info.split(b";")
    partitioned = (p.partition(b"=") for p in fields)
    return [(key, value) for key, _, value in partitioned]

def parse_vcf_sample(fmt, sample):
    "Parses VCF FORMAT and sample field into a list of key, value pairs"
    if sample == b".":
        return []
    return list(zip(fmt.split(b":"), sample.split(b":")))
