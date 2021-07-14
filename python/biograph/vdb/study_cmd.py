#!/usr/bin/env python3
'''
Manage VDB studies
'''
import argparse
import inspect
import multiprocessing
import subprocess
import sys
import tempfile

from pathlib import Path

import orjson as json

import biograph.vdb.athena as athena
from biograph.tools.log import setup_logging, log
from biograph.tools.refhash import refhash


class AthenaTableReader(multiprocessing.Process):
    '''
    Parallel TSV reader. Read gz chunks directly from S3 and write lines to outq.
    '''
    def __init__(self, inq, out_file, sample_names):
        super().__init__()
        self.inq = inq
        self.out_file = out_file
        self.sample_names = sample_names
        self.db = athena.connect(allow_db_create=False)

    def merge_samples(self, sample_json, strict=False):
        '''
        Generate merged format and sample fields.
        Returns the format string, a list of sample columns in sorted order,
        and the number of samples with data.
        '''
        samples = {}
        sample_data = json.loads(sample_json)

        for sample in sample_data:
            if sample_data[sample] is None:
                continue
            samples[sample] = sample_data[sample]

        # square-off VCFs may have no sample data at all
        if not samples:
            return ("GT", ".")

        unique_fields = {field for sample in samples for field in samples[sample]}

        format_fields = sorted(unique_fields)

        # GT is always first
        format_fields.remove('GT')
        format_fields.insert(0, 'GT')

        sample_column = []
        for sample in self.sample_names:
            if sample not in samples:
                if strict:
                    samples[sample] = {}
                else:
                    sample_column.append('.')
                    continue
            for field in format_fields:
                if field not in samples[sample]:
                    samples[sample][field] = '.'
            sample_column.append(':'.join([samples[sample][field] for field in format_fields]))

        return (':'.join(format_fields), sample_column)

    def run(self):
        with open(self.out_file, "w") as f:
            while True:
                data = self.inq.get()
                if data is None:
                    break

                (chrom, prefix) = data
                with self.db.download_gz_fh(prefix) as in_vcf:
                    for line in in_vcf:
                        (pos, varid, ref, alt, qual, filt, info, sample_column) = line.decode().rstrip().split('\t')
                        fmt, samples = self.merge_samples(sample_column)
                        print(
                            chrom,
                            pos,
                            varid or '.',
                            ref,
                            alt,
                            qual,
                            filt,
                            info or 'NS=0',
                            fmt,
                            '\t'.join(samples),
                            sep='\t',
                            file=f
                        )

def write_vcf(in_path, out_fh, tmp=tempfile.gettempdir(), **kwargs):
    '''
    Sort headerless VCF files from in_path and append to out_file using GNU sort
    '''
    args = [
        kwargs.get("gnusort", "/usr/bin/sort"),
        "-k1,1V" if kwargs.get("chrom_sort", False) else "-k1,1d",
        "-k2,2n",
        "-T", tmp
    ] + [str(f) for f in Path(in_path).glob("*")]

    psort = subprocess.Popen(
        args,
        stdout=out_fh
    )
    psort.wait()

def add_common_arguments(parser):
    ''' common arguments '''
    parser.add_argument("study_name", help="Name of the study")

def cmd_create(clargs):
    ''' Create a new study '''
    parser = argparse.ArgumentParser(prog=f"{CMD} create", description=inspect.getdoc(cmd_create),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_arguments(parser)

    args = parser.parse_args(clargs)
    db = athena.connect()
    db.create_study(args.study_name)
    print(f"Study '{args.study_name}' created")

def cmd_meta(clargs):
    ''' Describe a study '''
    raise SystemExit('Not implemented yet.')

def cmd_add(clargs):
    ''' Add variants to a study '''

    description = f"""{inspect.getdoc(cmd_add)}

Specify a VCF id or sample name to include all of its variants.
Wildcard matching * is applied to match multiple sample names.

To copy variants from the most recent checkpoint of an existing study,
use --from and specify one or more sample names with optional wildcards.
Use --checkpoint to select an older checkpoint in the study.

To remove VCFs from a study, use the 'filter' or 'revert' study commands.

All variants in a study must be called against the same reference.

Examples:

 # Add a specific VCF id
 $ biograph vdb study add my_study 0d1da4fa-778d-4d1d-9700-45f56acba576

 # Sample name
 $ biograph vdb study add my_study HG002

 # Wildcard match. Wrap in '' to avoid accidental shell glob matching.
 $ biograph vdb study add my_study 'HG00*' 'NA*3'

 # Copy all variants from an existing study at the most recent checkpoint
 $ biograph vdb study add my_study --from another_study '*'

 # Copy sample HG003 from an existing study at a specific checkpoint
 $ biograph vdb study add my_study --from another_study --checkpoint 3 'HG003'

"""
    parser = argparse.ArgumentParser(prog=f"{CMD} add", description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_arguments(parser)

    parser.add_argument("sample", nargs="+", help="VCF Sample name or aid to add")
    parser.add_argument("--from", dest="src_study", help="Look for samples in this study")
    parser.add_argument("--checkpoint", type=int, help="When using --from, copy variants form this checkpoint (default: most recent)")

    args = parser.parse_args(clargs)

    db = athena.connect()

    if args.src_study:
        db.copy_from_study(args.src_study, args.checkpoint, args.study_name, args.sample)
    else:
        if not args.sample:
            raise SystemExit('You must specify at least one sample, aid, or --from')

        if args.checkpoint:
            raise SystemExit('You must specify --from when using --checkpoint.')

        db.add_to_study(args.study_name, args.sample)

def cmd_show(clargs):
    ''' Show details about a study '''
    parser = argparse.ArgumentParser(prog=f"{CMD} show", description=inspect.getdoc(cmd_show),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_arguments(parser)
    args = parser.parse_args(clargs)

    db = athena.connect()
    db.assert_study_exists(args.study_name)

    meta = db.scalar(
        db.query(
            f"""
                SELECT CAST(MAP_AGG(key, value) AS JSON) AS meta
                FROM {db.table.study.meta}
                WHERE study_name = %(study_name)s
            ;
            """,
            params={"study_name": args.study_name},
        )
    )
    print(f"{'study_name':>16}:", args.study_name)
    print(f"{'created_on':>16}:", meta.get('created_on', '')[:19])

    for k in sorted(meta):
        if k == "created_on" or k.startswith("checkpoint"):
            continue
        print(f"{k:>16}:", meta[k])

    checkpoint = db.get_current_study_checkpoint(args.study_name)
    if not checkpoint:
        print("\nNo variants have been added to this study.")
        return

    print("\ncheckpoints:")
    for k in sorted(meta):
        if k.startswith("checkpoint"):
            print(f"{k[11:]:>4}:", meta[k])

    print(f"\n{'sample_name':<17}variant_count")

    for (sample, count) in db.query(f"""
        SELECT sample_name, COUNT(*)
        FROM {db.table.study.data}
        WHERE study_name = %(study_name)s
        AND checkpoint = %(checkpoint)d
        GROUP BY sample_name
        ORDER BY sample_name ASC
        ;
        """, params={"study_name": args.study_name, "checkpoint": checkpoint}):
        print(f"{sample:<17}{count}")

def cmd_export(clargs): # pylint: disable=too-many-statements
    ''' Export a study to a VCF file '''
    parser = argparse.ArgumentParser(prog=f"{CMD} export", description=inspect.getdoc(cmd_export),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", default="/dev/stdout", help="Write output VCF to this file (default: STDOUT)")
    parser.add_argument("-f", "--force", action="store_true", help="Overwrite local output directory without confirmation")
    parser.add_argument("-a", "--anno", default=None, help="Annotate the output with this annotation")
    parser.add_argument("-r", "--remerge", action="store_true", help="Force a merge prior to export, required when changing --fields (default: use pre-merged data if possible)")
    parser.add_argument("-t", "--tmp", type=str, default=tempfile.gettempdir(), help="Temporary directory (%(default)s)")
    parser.add_argument("-c", "--chromosomal", action="store_true", help="Use natural order (1,2,3,10,22,X) instead of alphabetic order (1,10,2,22,3,X)")
    parser.add_argument("--fields", help="List of FORMAT fields to export, separated by : (default: all fields)")
    parser.add_argument("--checkpoint", type=int, help="Export the study from this checkpoint (default: latest)")
    parser.add_argument("--square-off", help="Create a 'square-off' VCF with this single sample column")
    parser.add_argument("--no-header", action="store_true", help="Do not write a VCF header")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (%(default)s)")
    parser.add_argument("--sort", default="/usr/bin/sort", type=str, help=argparse.SUPPRESS)

    add_common_arguments(parser)
    args = parser.parse_args(clargs)

    db = athena.connect()
    db.assert_study_exists(args.study_name)

    out_file = Path(args.output)

    if str(out_file) != "/dev/stdout" and out_file.exists() and not args.force:
        raise SystemExit(f"Output path {out_file} already exists, refusing to overwrite.")

    checkpoint = args.checkpoint or db.get_current_study_checkpoint(args.study_name)
    sample_names = db.get_study_sample_names(args.study_name, checkpoint)

    if args.square_off:
        if args.square_off in sample_names:
            sample_names = [args.square_off]
        else:
            raise SystemExit(f"sample '{args.square_off}' is not present in {args.study_name} at checkpoint {checkpoint}.")

    try:
        (header_path, variants_path) = \
            db.merge_study(
                args.study_name,
                force_merge=args.remerge,
                anno_name=args.anno,
                square_off=args.square_off,
                checkpoint=checkpoint,
                format_fields=args.fields.split(':') if args.fields else None
            )
    except KeyboardInterrupt:
        raise SystemExit('\nAborted.')

    inq = multiprocessing.Queue()

    log("Downloading VDB data")
    with tempfile.TemporaryDirectory(prefix=f"{args.tmp}/") as tmpdir:
        chunk_dir = Path(tmpdir) / "chunks"
        chunk_dir.mkdir()
        out_vcf = open(out_file, "wb")

        if not args.no_header:
            db.download_fileobj(header_path, out_vcf)
            out_vcf.write(b'\t'.join([s.encode() for s in sample_names]))
            out_vcf.write(b'\n')
            out_vcf.flush()

        reader_threads = []
        for fn in range(max(1, args.threads)):
            reader = AthenaTableReader(inq, f"{chunk_dir}/{fn}", sample_names)
            reader.start()
            reader_threads.append(reader)

        rh = refhash(lookup=db.get_metadata_from_study(args.study_name, 'refname'))
        for gz in db.ls(variants_path, '.gz'):
            # db_name/study_name/merged/_export/study_name=the_study/chrom=1/junk_uuid.gz
            # Chroms are stored internally in ebi style, so convert to native
            chrom = rh.to_native(Path(gz).parts[-2].split('=')[1], rh.build(), rh.style())
            inq.put((chrom, gz))

        for rt in reader_threads:
            inq.put(None)

        for rt in reader_threads:
            rt.join()

        log("Exporting VCF")
        write_vcf(chunk_dir, out_vcf, tmp=tmpdir, chrom_sort=args.chromosomal, gnusort=args.sort)
        out_vcf.close()

def cmd_filter(clargs):
    ''' Filter variants in a study '''
    description = f"""{inspect.getdoc(cmd_filter)}

Filter variants in a study using bcftools filter syntax. A new study
checkpoint will be created.

Use --include to include variants that match the filter.

Use --exclude to exclude variants that match the filter.

Examples:

 # PASS only
 $ biograph vdb study filter my_study --exclude "FILTER != 'PASS'"

 # High quality hets on chr16
 $ biograph vdb study filter my_study --include "chrom = '16' AND GT = 0/1 AND qual > 50"

 # Per-variant missingness
 $ biograph vdb study filter my_study --include "F_MISS > 0.2"

 # Per-sample missingness
 $ biograph vdb study filter my_study --exclude "SAMPLE_MISS > 0.1"

"""
    parser = argparse.ArgumentParser(prog=f"{CMD} filter", description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--include", help="Include only variants that match these criteria")
    group.add_argument("-e", "--exclude", help="Exclude all variants that match these criteria")

    add_common_arguments(parser)
    args = parser.parse_args(clargs)

    db = athena.connect()
    db.assert_study_exists(args.study_name)

    db.filter_study(
        study_name=args.study_name,
        the_filter=args.include or args.exclude,
        exclude=args.include is None
    )

def cmd_list(clargs):
    ''' List all available studies '''
    parser = argparse.ArgumentParser(prog=f"{CMD} list", description=inspect.getdoc(cmd_list),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.parse_args(clargs)

    db = athena.connect()
    print(f"{'study_name':<21} {'created_on':<21}")

    for (study_name, meta) in db.query(f"SELECT study_name, CAST(MAP_AGG(key, value) AS JSON) FROM {db.table.study.meta} GROUP BY study_name ORDER BY study_name ASC;"):
        print(f"{study_name:<21} {meta.get('created_on', '')[:19]:<21}")

def cmd_freeze(clargs):
    ''' Prevent changes to a study '''
    parser = argparse.ArgumentParser(prog=f"{CMD} freeze", description=inspect.getdoc(cmd_freeze),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_arguments(parser)

    args = parser.parse_args(clargs)
    db = athena.connect()
    db.study_freeze(args.study_name)
    print(f"Study '{args.study_name}' frozen")

def cmd_unfreeze(clargs):
    ''' Allow changes to a study '''
    parser = argparse.ArgumentParser(prog=f"{CMD} unfreeze", description=inspect.getdoc(cmd_unfreeze),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_arguments(parser)

    args = parser.parse_args(clargs)
    db = athena.connect()
    db.study_unfreeze(args.study_name)
    print(f"Study '{args.study_name}' unfrozen")

def cmd_delete(clargs):
    ''' Delete a study '''
    parser = argparse.ArgumentParser(prog=f"{CMD} delete", description=inspect.getdoc(cmd_delete),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_arguments(parser)
    args = parser.parse_args(clargs)
    db = athena.connect()
    db.delete_study(args.study_name)

    print(f"Study '{args.study_name}' deleted")

def cmd_revert(clargs):
    ''' Revert to a previous checkpoint '''
    parser = argparse.ArgumentParser(prog=f"{CMD} revert", description=inspect.getdoc(cmd_revert),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--checkpoint", type=int, help="Revert to this checkpoint (default: roll back one)")
    add_common_arguments(parser)
    args = parser.parse_args(clargs)
    db = athena.connect()

    current_checkpoint = db.get_current_study_checkpoint(args.study_name)

    if current_checkpoint == 0:
        raise SystemExit(f"No checkpoints yet in study {args.study_name}")
    if args.checkpoint:
        if args.checkpoint < 0:
            raise SystemExit(f"Invalid checkpoint {args.checkpoint}")
        if current_checkpoint < args.checkpoint:
            raise SystemExit(f"No checkpoint {args.checkpoint} in {args.study_name} (max {current_checkpoint})")
        if current_checkpoint == args.checkpoint:
            raise SystemExit(f"Study {args.study_name} already at checkpoint {current_checkpoint}, nothing to do.")
        target_checkpoint = args.checkpoint
    else:
        target_checkpoint = current_checkpoint - 1

    for chkpt in range(current_checkpoint, target_checkpoint, -1):
        db.delete_study(args.study_name, chkpt)

    print(f"Study '{args.study_name}' reverted to checkpoint {target_checkpoint}")

def main(clargs):
    ''' Top level parser '''
    usage = f'''study [COMMAND] [options]

Manage studies in the Spiral Variant DataBase (VDB).

Run any command with --help for additional information.

    create    {inspect.getdoc(CMDS['create'])}

    list      {inspect.getdoc(CMDS['list'])}
    show      {inspect.getdoc(CMDS['show'])}

    add       {inspect.getdoc(CMDS['add'])}
    filter    {inspect.getdoc(CMDS['filter'])}

    export    {inspect.getdoc(CMDS['export'])}

    freeze    {inspect.getdoc(CMDS['freeze'])}
    unfreeze  {inspect.getdoc(CMDS['unfreeze'])}

    revert    {inspect.getdoc(CMDS['revert'])}
    delete    {inspect.getdoc(CMDS['delete'])}

'''
    parser = argparse.ArgumentParser(prog="study", usage=usage,
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
CMD = 'biograph vdb study'

# module global CMDs
CMDS = {
    'create':   cmd_create,
    'add':      cmd_add,
    'filter':   cmd_filter,
    'list':     cmd_list,
    'show':     cmd_show,
    'delete':   cmd_delete,
    'freeze':   cmd_freeze,
    'unfreeze': cmd_unfreeze,
    'export':   cmd_export,
    'revert':   cmd_revert,
}

if __name__ == '__main__':
    try:
        main(sys.argv[1:])
    except KeyboardInterrupt:
        raise SystemExit('\nAborted.')
