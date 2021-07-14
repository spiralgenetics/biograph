#!/usr/bin/env python3
'''
VDB: VCF commands
'''
import argparse
import gzip
import inspect
import multiprocessing
import shutil
import subprocess
import sys
import tempfile
import uuid

from pathlib import Path

import pyarrow.parquet as pq
import orjson as json

import biograph.vdb as vdb
import biograph.vdb.athena as athena
from biograph.tools.log import setup_logging, log
from biograph.tools.refhash import refhash
from biograph.utils import (get_opener, timestamp)

def check_record_count(ds_path, expected):
    '''
    Verify that the parquet dataset contains the specified record count.
    '''
    ds = pq.ParquetDataset(ds_path)
    pq_count = ds.read().num_rows
    if pq_count != expected:
        raise SystemExit(f"Aborting: expected {expected} lines but found {pq_count} in parquet")

def cmd_export(clargs):
    ''' Export a VCF '''
    parser = argparse.ArgumentParser(prog=f"{CMD} export", description=inspect.getdoc(cmd_export),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--aid", required=True, type=str, help="Export the VCF with this AID")
    parser.add_argument("-o", "--output", default="/dev/stdout", type=str, help="Write output to this file (default: STDOUT)")

    args = parser.parse_args(clargs)
    db = athena.connect()

    db.download_aid(args.aid, args.output, db.path.vcf.files)

def cmd_list(clargs):
    ''' List all available VCFs '''
    parser = argparse.ArgumentParser(prog=f"{CMD} list", description=inspect.getdoc(cmd_list),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-s", "--sample", type=str, help="Only list VCFs with this sample name prefix")

    args = parser.parse_args(clargs)

    db = athena.connect()
    print(f"{'sample_name':<21}{'imported_on':<21}{'refname':<11}{'variant_count':<14}{'aid':<37}description")

    if args.sample:
        query = f"SELECT sample_name, imported_on, refname, variant_count, aid, description FROM {db.table.vcf.meta} WHERE sample_name LIKE '{args.sample}%' ORDER BY sample_name ASC, imported_on DESC;"
    else:
        query = f"SELECT sample_name, imported_on, refname, variant_count, aid, description FROM {db.table.vcf.meta} ORDER BY sample_name ASC, imported_on DESC;"

    for (sample_name, imported_on, refname, variant_count, aid, description) in db.query(query):
        print(f"{sample_name:<21}{str(imported_on)[:19]:<21}{refname:<11}{variant_count:<14}{aid:<37}{description}")

def cmd_delete(clargs):
    ''' Delete a VCF '''
    parser = argparse.ArgumentParser(prog=f"{CMD} delete", description=inspect.getdoc(cmd_delete),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("aid", nargs="+", type=str, help="Delete the VCF with this aid")

    args = parser.parse_args(clargs)
    db = athena.connect()

    db.delete_aid(args.aid)

def sort_vcf(in_vcf, out_vcf, tmp=tempfile.gettempdir(), **kwargs):
    '''
    Sort VCF from filehandle in_vcf to filehandle out_vcf using GNU sort and split.
    '''
    chunk_dir = Path(tmp) / "chunks"
    chunk_dir.mkdir()

    line = in_vcf.readline()
    while line.startswith(b"#"):
        out_vcf.write(line)
        line = in_vcf.readline()
    out_vcf.flush()

    # ~128MB chunks seems optimal, but the value isn't critical.
    # -C breaks on \n boundaries
    psplit = subprocess.Popen(
        [
            kwargs.get("gnusplit", "/usr/bin/split"),
            "-C", "128000000",
            "-d", "-",
            f"{chunk_dir}/x"
        ],
        stdin=subprocess.PIPE,
    )
    psplit.stdin.write(line)
    shutil.copyfileobj(in_vcf, psplit.stdin)
    psplit.stdin.close()
    psplit.wait()

    psort = subprocess.Popen(
        [
            kwargs.get("gnusort", "/usr/bin/sort"),
            "-k1,1V" if kwargs.get("chrom_sort", False) else "-k1,1d",
            "-k2,2n",
            "-T", tmp
        ] + [str(f) for f in chunk_dir.glob("*")],
        stdout=out_vcf,
    )
    psort.wait()

def cmd_sort(clargs):
    ''' Sort a VCF file '''
    parser = argparse.ArgumentParser(
        prog=f"{CMD} sort",
        description='''Sort a VCF file.'''
    )
    parser.add_argument("-i", "--input", default="/dev/stdin", type=str, help="Input VCF filename (%(default)s)")
    parser.add_argument("-o", "--output", default="/dev/stdout", type=str, help="Output VCF filename (%(default)s)")
    parser.add_argument("-c", "--chromosomal", action="store_true", help="Use natural order (1,2,3,10,22,X) instead of alphabetic order (1,10,2,22,3,X)")
    parser.add_argument("-t", "--tmp", type=str, default=tempfile.gettempdir(), help="Temporary directory (%(default)s)")
    parser.add_argument("--sort", default="/usr/bin/sort", type=str, help=argparse.SUPPRESS)
    parser.add_argument("--split", default="/usr/bin/split", type=str, help=argparse.SUPPRESS)

    args = parser.parse_args(clargs)

    if sys.stdin.isatty() and args.input == '/dev/stdin':
        log("Refusing to read VCF from the terminal.")
        parser.print_help(sys.stderr)
        exit(1)

    if args.input.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    with tempfile.TemporaryDirectory(prefix=f"{args.tmp}/") as tmpdir:
        with opener(args.input, "rb") as in_vcf:
            with open(args.output, "wb") as out_vcf:
                sort_vcf(in_vcf, out_vcf, tmp=tmpdir, chrom_sort=args.chromosomal, gnusort=args.sort, gnusplit=args.split)

def add_import_args(clargs):
    ''' Add import arguments and validate them '''
    parser = argparse.ArgumentParser(
        prog=f"{CMD} import",
        description='''Import a VCF file. It may be optionally gzip/bgzip compressed.'''
    )
    parser.add_argument("input", nargs='?', type=str, help="Input VCF filename")
    parser.add_argument("-s", "--sample", type=str,
                        help="Sample name for variants (default: extract from VCF file)")
    parser.add_argument("-d", "--description", type=str, default="",
                        help="Free form description of this VCF")
    parser.add_argument("-r", "--refhash", type=str,
                        help="Explicit reference name or hash (default: extract from input file)")
    parser.add_argument("--aid", type=str, default=str(uuid.uuid4()),
                        help="Unique GUID (default: autogenerate)")
    parser.add_argument("-o", "--output", type=str, default='.',
                        help="Output directory prefix (%(default)s)")
    parser.add_argument("--keep-data", action='store_true',
                        help="Keep a local copy of the converted import data (default: delete after upload)")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Maximum number of threads (%(default)s)")
    parser.add_argument("--tmp", type=str, default=tempfile.gettempdir(), help="Temporary directory (%(default)s)")

    args = parser.parse_args(clargs)

    if len(sys.argv) == 4:
        raise SystemExit(parser.print_help())

    if args.input in ('/dev/stdin', '-', 'STDIN'):
        raise SystemExit('Reading from STDIN is not supported.')

    return args

def get_count(p):
    ''' Return the number of rows in a parquet dataset '''
    return pq.ParquetDataset(str(p)).read().num_rows

def cmd_import(clargs):
    ''' Import a VCF '''
    args = add_import_args(clargs)
    db = athena.connect()
    aid = db.validate_aid(args.aid)

    local_path = Path(args.output) / aid
    try:
        local_path.mkdir(parents=True)
    except FileExistsError:
        raise SystemExit(f"Output path {local_path} already exists. Choose a different --output or --aid.")

    # Print AID to stdout for easy shell capture. All other logs are to stderr.
    print(aid)
    log(f'Importing {args.input}')

    tmp_path = local_path / "tmp"
    tmp_path.mkdir(parents=True)

    if args.input.lower().startswith("s3://"):
        in_file = db.download_s3_path(args.input, local_path)
    else:
        in_file = args.input

    # refhash tries the lookup if available, otherwise computes it from filename
    rh = refhash(lookup=args.refhash, filename=in_file)

    with get_opener(in_file) as in_vcf:
        headers = vdb.vcf_to_parquet(
            rh.build(),
            in_vcf,
            tmp_path,
            nthreads=args.threads,
            tmpdir=args.tmp
        )

    # We have to do this parquet file path dance because we might not know the
    # sample name until the headers have been processed.

    sample_name = db.validate_sample_name(args.sample or headers[-1].split()[-1])

    path_suffix = Path(f"sample_name={sample_name}") / f"build={rh.build()}" / f"aid={aid}"

    pq_path = local_path / db.path.vcf.data / path_suffix
    pq_path.parent.mkdir(parents=True)
    tmp_path.rename(pq_path)

    log("Conversion complete. Verifying parquet.")

    # These can get pretty big, and they count fairly quickly, so limit to max 8
    with multiprocessing.Pool(min(8, args.threads)) as p:
        pq_count = sum(p.map(get_count, [f for f in pq_path.glob('*')]))

    log(f"{pq_count} variants converted to parquet")

    headers_path = local_path / db.path.vcf.meta / path_suffix / 'headers.json'
    headers_path.parent.mkdir(parents=True)

    with open(f"{headers_path}", "w") as out_json:
        print(json.dumps(
            {
                "build": rh.build(),
                "refname": rh.common_name(),
                "refhash": rh.digest(),
                "header": '\n'.join(headers),
                "sample_name": sample_name,
                "aid": aid,
                "description": args.description,
                "imported_on": timestamp(),
                "variant_count": pq_count,
                "filename": Path(in_file).name
            }
        ).decode(), file=out_json)

    if pq_count == 0:
        shutil.rmtree(local_path)
        raise SystemExit(f"No variants found in {in_file}. Nothing to do.")

    log(f"Uploading to S3...")

    db.sync_to_s3(local_path / db.path.vcf.root, db.path.vcf.root)
    db.upload_to_s3(in_file, db.path.vcf.files / path_suffix / 'original.vcf.gz')
    db.add_vcf_partitions(sample_name, rh.build(), aid)

    # sanity check: remote pq matches local pq
    log(f'Validating upload...')

    vdb_count = db.get_vcf_variant_count(aid)
    if vdb_count != pq_count:
        # TODO: cleanup local fs, s3, partition
        raise SystemExit(f"Variant count in VDB ({vdb_count}) does not match local ({pq_count}), aborting.")

    log(f"VCF imported, aid={aid}")

    if args.keep_data:
        log(f'Import data retained in {local_path}/')
    else:
        shutil.rmtree(local_path)

def main(clargs):
    ''' Top level parser '''
    usage = f'''vcf [COMMAND] [options]

Import and export VCF data.

Run any command with --help for additional information.

    import    {inspect.getdoc(CMDS['import'])}
    export    {inspect.getdoc(CMDS['export'])}

    list      {inspect.getdoc(CMDS['list'])}
    delete    {inspect.getdoc(CMDS['delete'])}

    sort      {inspect.getdoc(CMDS['sort'])}

'''
    parser = argparse.ArgumentParser(prog="vcf", usage=usage,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="COMMAND", choices=CMDS.keys(), type=str, help=argparse.SUPPRESS)
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    if len(sys.argv) == 3:
        raise SystemExit(parser.print_help())

    args = parser.parse_args(clargs)

    setup_logging(debug_mode=args.debug, simple=True)

    CMDS[args.cmd](args.options)

# top level command
CMD = 'biograph vdb vcf'

# module global CMDs
CMDS = {
    'import':   cmd_import,
    'export':   cmd_export,
    'list':     cmd_list,
    'delete':   cmd_delete,
    'sort':     cmd_sort,
    # 'show':    cmd_show,
}

if __name__ == '__main__':
    main(sys.argv[1:])
