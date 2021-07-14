#!/usr/bin/env python3
'''
VDB: Annotation commands
'''
import argparse
import inspect
import multiprocessing
import shutil
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
    if pq_count == expected:
        return True

    log(f"Expected {expected} lines but found {pq_count} in parquet")
    return False

def cmd_export(clargs):
    ''' Export an annotation '''
    parser = argparse.ArgumentParser(prog=f"{CMD} export", description=inspect.getdoc(cmd_export),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--aid", required=True, type=str, help="Export the annotation with this AID")
    parser.add_argument("-o", "--output", default="/dev/stdout", type=str, help="Write output to this file (default: STDOUT)")
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args(clargs)
    db = athena.connect()

    db.download_aid(args.aid, args.output, db.path.anno.files)

def cmd_list(clargs):
    ''' List all available annotations '''
    parser = argparse.ArgumentParser(prog=f"{CMD} list", description=inspect.getdoc(cmd_list),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.parse_args(clargs)

    db = athena.connect()
    print(f"{'anno_name':<21}{'version':<13}{'imported_on':<21}{'build':<11}{'annotations':<14}{'aid':<37}description")

    query = f"SELECT anno_name, version, imported_on, build, variant_count, aid, description FROM {db.table.anno.meta} ORDER BY anno_name ASC, imported_on DESC;"

    for (anno_name, version, imported_on, refname, variant_count, aid, description) in db.query(query):
        print(f"{anno_name:<21}{version:<13}{str(imported_on)[:19]:<21}{refname:<11}{variant_count:<14}{aid:<37}{description}")

def cmd_delete(clargs):
    ''' Delete an annotation '''
    parser = argparse.ArgumentParser(prog=f"{CMD} delete", description=inspect.getdoc(cmd_delete),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("aid", nargs="+", type=str, help="Delete the annotation with this aid")
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args(clargs)
    db = athena.connect()

    db.delete_aid(args.aid, "anno")

def add_import_args(clargs):
    ''' Add import arguments and validate them '''
    parser = argparse.ArgumentParser(
        prog=f"{CMD} import",
        description='''Import variant annotation data in GFF, GTF, or VCF format.'''
    )
    parser.add_argument("anno_name", help="Name for this annotation (eg. ClinVar)")
    parser.add_argument("version", help="Version of this annotation (eg. 2020-10-03)")

    parser.add_argument("-i", "--input", default="/dev/stdin", type=str, help="Input filename (default: STDIN)")
    parser.add_argument("-f", "--format", default=None, choices=('vcf', 'gtf', 'gff'),
                        help="Input contains annotations in this format")
    parser.add_argument("-r", "--refhash", type=str,
                        help="Explicit reference name or hash (default: extract from input file if possible)")
    parser.add_argument("-d", "--description", type=str, default="",
                        help="Free form description of this annotation")
    parser.add_argument("-o", "--output", type=str, default='.',
                        help="Output directory prefix (%(default)s)")
    parser.add_argument("--aid", type=str, default=str(uuid.uuid4()),
                        help="Unique GUID (default: autogenerate)")
    parser.add_argument("--keep-data", action='store_true',
                        help="Keep a local copy of the converted import data (default: delete after upload)")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                        help="Maximum number of threads (%(default)s)")
    parser.add_argument("--tmp", type=str, default=tempfile.gettempdir(), help="Temporary directory (%(default)s)")
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args(clargs)

    if len(sys.argv) == 4:
        raise SystemExit(parser.print_help())

    if args.format is None:
        if args.input.lower().endswith(('vcf', 'vcf.gz')):
            args.format = 'vcf'
        elif args.input.lower().endswith(('gtf', 'gtf.gz')):
            args.format = 'gtf'
        elif args.input.lower().endswith(('gff', 'gff.gz', 'gff3', 'gff3.gz')):
            args.format = 'gff'
        else:
            raise SystemExit(f"Could not determine the file format. Please specify --format.")

    return args

def get_count(p):
    ''' Return the number of rows in a parquet dataset '''
    return pq.ParquetDataset(str(p)).read().num_rows

def cmd_import(clargs):
    ''' Import an annotation file '''
    args = add_import_args(clargs)
    db = athena.connect()
    aid = db.validate_aid(args.aid)

    # refhash tries the lookup if available, otherwise computes it from filename
    rh = refhash(lookup=args.refhash, filename=args.input)
    # build may be blank (no contigs present) or unknown (never seen before)
    if rh.build() in ('', 'unknown'):
        raise SystemExit(f"Could not identify the reference build from {args.input}. Please specify --refhash")

    local_path = Path(args.output) / aid
    try:
        local_path.mkdir(parents=True)
    except FileExistsError:
        raise SystemExit(f"Output path {local_path} already exists. Choose a different --output or --aid.")

    # Print AID to stdout for easy shell capture. All other logs are to stderr.
    print(aid)
    log(f'Importing from {args.input} for build {rh.build()}')

    path_suffix = Path(f"build={rh.build()}") / f"anno_name={args.anno_name}" / f"version={args.version}" / f"aid={aid}"

    pq_path = local_path / db.path.anno.data / path_suffix
    pq_path.mkdir(parents=True)

    with get_opener(args.input) as in_vcf:
        headers = vdb.anno_to_parquet(
            args.format,
            rh.build(),
            in_vcf,
            pq_path,
            nthreads=args.threads,
            tmpdir=args.tmp,
            verbose=args.debug
        )

    log("Conversion complete. Verifying parquet.")

    # These can get pretty big, and they count fairly quickly, so limit to max 8
    with multiprocessing.Pool(min(8, args.threads)) as p:
        pq_count = sum(p.map(get_count, [f for f in pq_path.glob('*')]))

    log(f"{pq_count} annotations converted to parquet")

    headers_path = local_path / db.path.anno.meta / path_suffix / 'headers.json'
    headers_path.parent.mkdir(parents=True)

    with open(f"{headers_path}", "w") as out_json:
        print(json.dumps(
            {
                "build": rh.build(),
                "refname": rh.common_name(),
                "refhash": rh.digest(),
                "header": '\n'.join(headers),
                "anno_name": args.anno_name,
                "aid": aid,
                "description": args.description,
                "imported_on": timestamp(),
                "variant_count": pq_count,
                "filename": Path(args.input).name
            }
        ).decode(), file=out_json)

    if pq_count == 0:
        shutil.rmtree(local_path)
        raise SystemExit(f"No variants found in {args.input}. Nothing to do.")

    log(f"Uploading to S3...")

    db.sync_to_s3(local_path / db.path.anno.root, db.path.anno.root)
    db.upload_to_s3(args.input, db.path.anno.files / path_suffix / 'original.vcf.gz')
    db.add_anno_partitions(rh.build(), args.anno_name, args.version, aid)

    # sanity check: remote pq matches local pq
    log(f'Validating upload...')

    if db.get_anno_variant_count(aid) != pq_count:
        # TODO: cleanup local fs, s3, partition
        raise SystemExit('Annotation count in VDB does not match local, aborting.')

    log(f"Annotation imported, aid={aid}")

    if args.keep_data:
        log(f'Import data retained in {local_path}/')
    else:
        shutil.rmtree(local_path)

def main(clargs):
    ''' Top level parser '''
    usage = f'''anno [COMMAND] [options]

Import and export variant annotation data.

Run any command with --help for additional information.

    import    {inspect.getdoc(CMDS['import'])}
    export    {inspect.getdoc(CMDS['export'])}

    list      {inspect.getdoc(CMDS['list'])}
    delete    {inspect.getdoc(CMDS['delete'])}

'''
    parser = argparse.ArgumentParser(prog="anno", usage=usage,
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
CMD = 'biograph vdb anno'

# module global CMDs
CMDS = {
    'import':   cmd_import,
    'export':   cmd_export,
    'list':     cmd_list,
    'delete':   cmd_delete,
}

if __name__ == '__main__':
    main(sys.argv[1:])
