"""

Convert variants files (e.g. VCF) to parquet format.

NOTE: make sure that any parquet schema changes here are reflected in
sql_tables.py

"""
import tempfile

import pyarrow as pa
import numpy as np

from .parquet import ConverterBase, parse_vcf_info, parse_vcf_sample

VARIANTS_ARROW_SCHEMA = pa.schema([
    pa.field("spans", pa.list_(pa.int64()), nullable=False),
    pa.field("reflen", pa.int64(), nullable=False),
    pa.field("chrom", pa.string(), nullable=False),
    pa.field("pos", pa.int64(), nullable=False),
    pa.field("varend", pa.int64(), nullable=False),
    pa.field("varid", pa.string()),
    pa.field("ref", pa.string(), nullable=False),
    pa.field("alt", pa.string(), nullable=False),
    pa.field("qual", pa.float32(), nullable=True),
    pa.field("filt", pa.string(), nullable=True),
    pa.field("info", pa.map_(pa.string(), pa.string())),
    pa.field("sample", pa.map_(pa.string(), pa.string())),
])

class VcfConverter(ConverterBase):
    """
    Converter to convert VCF to a parquet file.
    """
    def __init__(self, build, tmpdir, outfile, **kwargs):
        super().__init__(build, tmpdir, outfile, VARIANTS_ARROW_SCHEMA, **kwargs)

    def parse_chunk(self, fh):
        "Reads lines from the given chunk and parses into the schema"
        for line in fh:
            fields = line.rstrip().split(b"\t")
            if len(fields) > 10:
                raise RuntimeError(f"Only single sample VCFs are currently supported (this VCF contains {len(fields) - 9})")
            if len(fields) != 10:
                raise RuntimeError(f"Expecting 10 fields in vcf line: {repr(fields)}")
            chrom, pos, varid, ref, alt, qual, filt, info, t_format, sample = fields
            if b',' in alt:
                raise RuntimeError("Multiallelic VCFs are not supported. Normalize with 'bcftools norm -m-any'")
            start = int(pos) - 1 # Convert from 1-based to 0-based
            end = start + len(ref)
            yield ((), # extra row group keys
                   (chrom,
                    start,
                    end,
                    np.nan if varid == b"." else varid,
                    ref,
                    alt,
                    np.nan if qual == b"." else float(qual),
                    np.nan if filt == b"." else filt,
                    parse_vcf_info(info),
                    parse_vcf_sample(t_format, sample),
                   ))

def vcf_to_parquet(build, fh, outfile, *, tmpdir=None, **kwargs):
    "Converts the given VCF data to parquet format.  Returns a list of lines in the VCF header"
    with tempfile.TemporaryDirectory(dir=tmpdir) as td:
        v = VcfConverter(build, td, outfile, **kwargs)
        v.make_chunks_from(fh)
        v.generate()
        return v.header_lines
