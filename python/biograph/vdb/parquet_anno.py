"""

Convert annotation files (e.g. GFF, sample-less VCF) to parquet format.

NOTE: make sure that any parquet schema changes here are reflected in
sql_tables.py

"""
import tempfile
import urllib.parse

import pyarrow as pa
import numpy as np

from .parquet import ConverterBase, parse_vcf_info

ANNO_ARROW_SCHEMA = pa.schema([
    pa.field("spans", pa.list_(pa.int64()), nullable=False),
    pa.field("reflen", pa.int64(), nullable=False),
    pa.field("chrom", pa.string(), nullable=False),
    pa.field("pos", pa.int64(), nullable=False),
    pa.field("varend", pa.int64(), nullable=False),
    # VCF fields:
    pa.field("varid", pa.string()),
    pa.field("ref", pa.string()),
    pa.field("alt", pa.string()),
    pa.field("qual", pa.float32()),
    pa.field("filt", pa.string()),
    pa.field("info", pa.map_(pa.string(), pa.string())),
    # GFF fields:
    pa.field("source", pa.string()),
    pa.field("feature", pa.string()),
    pa.field("score", pa.float32()),
    pa.field("frame", pa.string()),
    pa.field("strand", pa.string()),
    pa.field("attributes", pa.map_(pa.string(), pa.string()))
])

class VcfAnnoConverter(ConverterBase):
    """
    Converter to convert a sample-less VCF to a parquet file.
    """
    def __init__(self, build, tmpdir, outfile, **kwargs):
        super().__init__(build, tmpdir, outfile, ANNO_ARROW_SCHEMA, **kwargs)

    def parse_chunk(self, fh):
        "Reads lines from the given chunk and parses into the schema"
        for line in fh:
            fields = line.rstrip().split(b"\t")
            if len(fields) != 8:
                raise RuntimeError(f"Expecting 8 fields in vcf annotation line: {repr(fields)}")
            chrom, pos, varid, ref, alt, qual, filt, info = fields
            if b',' in alt:
                raise RuntimeError("Multiallelic VCFs are not supported. Normalize with 'bcftools norm -m-any'")
            start = int(pos) - 1 # Convert from 1-based to 0-based
            end = start + len(ref)
            yield ((), # extra row group keys
                   (chrom,
                    start,
                    end,
                    # VCF fields:
                    np.nan if varid == b"." else varid,
                    ref,
                    alt,
                    np.nan if qual == b"." else float(qual),
                    np.nan if filt == b"." else filt,
                    parse_vcf_info(info),
                    # GFF fields:
                    np.nan, # source
                    np.nan, # feature
                    np.nan, # score
                    np.nan, # frame
                    np.nan, # strand
                    np.nan # attributes
                   ))

class GtfAnnoConverter(ConverterBase):
    """
Converter to convert a GTF (GFFv2) annotation file to a parquet file.
"""
    def __init__(self, build, tmpdir, outfile, **kwargs):
        super().__init__(
            build, tmpdir, outfile,
            ANNO_ARROW_SCHEMA, **kwargs)

    def parse_chunk(self, fh):
        "Reads lines from the given chunk and parses into the schema"
        for line in fh:
            fields = line.rstrip().split(b"\t")
            if len(fields) != 9:
                raise RuntimeError(f"Expecting 9 fields in GTF annotation line: {repr(fields)}")
            chrom, source, method, start, end, score, strand, phase, group = fields
            start = int(start) - 1 # Convert from 1-based to 0-based
            end = int(end) - 1 # Convert from 1-based to 0-based
            if end < start:
                raise RuntimeError(f"End must not be before start in annotation line: {repr(fields)}")

            varid = None
            attrs = {}
            for field in group.split(b";"):
                if not field:
                    continue
                key, _, val = field.strip().partition(b" ")
                val = urllib.parse.unquote(val.decode("UTF8"))
                attrs[key] = val.strip('"')

            # pick the first best field for varid
            for key in (b"Name", b"gene_name", b"gene_id", b"logic_name", b"ID", b"Parent"):
                if key in attrs:
                    varid = attrs[key]
                    break

            yield ((method,), # extra row group keys
                   (chrom,
                    start,
                    end,
                    # VCF fields:
                    np.nan if varid is None else varid, # varid
                    np.nan, # ref
                    np.nan, # alt
                    np.nan, # qual
                    np.nan, # filt
                    np.nan, # info
                    # GFF fields:
                    source, # source
                    method, # feature
                    np.nan if score == b"." else float(score),
                    np.nan if phase == b"." else phase, # frame
                    np.nan if strand == b"." else strand, # strand
                    np.nan if group == b"." else [(k, attrs[k]) for k in attrs], # attributes
                   ))

class GffAnnoConverter(ConverterBase):
    """
Converter to convert a GFFv3 annotation file to a parquet file.
"""
    def __init__(self, build, tmpdir, outfile, **kwargs):
        super().__init__(
            build, tmpdir, outfile,
            ANNO_ARROW_SCHEMA, **kwargs)

    def parse_chunk(self, fh):
        "Reads lines from the given chunk and parses into the schema"
        for line in fh:
            # Annoyingly, gff may contain inline comments. Just skip them.
            if line.startswith(b"#"):
                continue
            fields = line.rstrip().split(b"\t")
            if len(fields) != 9:
                raise RuntimeError(f"Expecting 9 fields in GFF annotation line: {repr(fields)}")
            chrom, source, stype, start, end, score, strand, phase, attrstr = fields
            start = int(start) - 1 # Convert from 1-based to 0-based
            end = int(end) - 1 # Convert from 1-based to 0-based
            if end < start:
                raise RuntimeError(f"End must not be before start in annotation line: {repr(fields)}")

            varid = None
            if attrstr == b".":
                attrs = np.nan
            else:
                attrs = {}
                for field in attrstr.split(b";"):
                    if not field:
                        continue
                    key, _, val = field.strip().partition(b"=")
                    val = urllib.parse.unquote(val.decode("UTF8"))
                    attrs[key] = val

            # pick the first best field for varid
            for key in (b"Name", b"gene_name", b"gene_id", b"logic_name", b"ID", b"Parent"):
                if key in attrs:
                    varid = attrs[key]
                    break

            yield ((stype,), # extra row group keys
                   (chrom,
                    start,
                    end,
                    # VCF fields:
                    np.nan if varid is None else varid, # varid
                    np.nan, # ref
                    np.nan, # alt
                    np.nan, # qual
                    np.nan, # filt
                    np.nan, # info
                    # GFF fields:
                    source, # source
                    stype, # feature
                    np.nan if score == b"." else float(score),
                    np.nan if phase == b"." else phase, # frame
                    np.nan if strand == b"." else strand, # strand
                    [(k, attrs[k]) for k in attrs], # attributes
                   ))

ANNO_CONVERTERS = {
    "vcf": VcfAnnoConverter,
    "gtf": GtfAnnoConverter,
    "gff": GffAnnoConverter,
}

def anno_to_parquet(fmt, build, fh, outfile, *, tmpdir=None, **kwargs):
    "Converts the given annotation data to parquet format.  Returns a list of lines in the VCF header"
    if fmt not in ANNO_CONVERTERS:
        raise RuntimeError(f"Unknown annotation format {fmt}; valid formats: {', '.join(ANNO_CONVERTERS.keys())}")

    with tempfile.TemporaryDirectory(dir=tmpdir) as td:
        v = ANNO_CONVERTERS[fmt](build, td, outfile, **kwargs)
        v.make_chunks_from(fh)
        v.generate()
        return v.header_lines
