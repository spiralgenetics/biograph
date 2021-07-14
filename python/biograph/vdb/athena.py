'''
Wrapper interface for the VDB Athena backend.

All vdb operations use the following environment variable overrides:

      VDB_DB: your VDB database name (default: vdb)
  VDB_BUCKET: your VDB S3 endpoint (default: s3://spiral-vdb)

You may also optionally set AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, and
AWS_DEFAULT_REGION. This is not required if you are using vdb from an instance
deployed with a role with sufficient privileges.
'''
import gzip
import os
import re
import uuid

from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
from time import sleep, time
from types import SimpleNamespace

import boto3
import orjson as json
import pandas as pd

from pyathena.connection import Connection
from pyathena.pandas.async_cursor import AsyncPandasCursor

from biograph.utils import timestamp, typed, plural, confirm, chunked
from biograph.vdb import create_table_sql
from biograph.vdb.cache import fetch_from_cache, get_table_mtime, update_table_mtime, clear_table_mtime
from biograph.vdb.filter import parser, ParseException

from biograph.tools.refhash import refhash
from biograph.tools.log import debug, log, error

class connect: # pylint: disable=too-many-lines
    '''
    Wrapper for the Athena VDB backend
    '''
    def __init__(
            self,
            database=None,
            bucket=None,
            aws_region=None,
            aws_access_key_id=os.environ.get('AWS_ACCESS_KEY_ID', None),
            aws_secret_access_key=os.environ.get('AWS_SECRET_ACCESS_KEY', None),
            allow_db_create=True
        ):
        ''' Set up the vdb connection '''

        self.aws_region = aws_region or os.environ.get('AWS_DEFAULT_REGION', 'us-west-2')
        os.environ['AWS_DEFAULT_REGION'] = self.aws_region

        self.s3 = boto3.resource('s3')
        self.athena = boto3.client('athena')

        self.database = self.get_database_name(database or os.environ.get('VDB_DB', 'main'))

        self.bucket = bucket or os.environ.get('VDB_BUCKET', 'vdb-demo')
        self.bucket = self.bucket.rstrip('/')

        # root should never equal 'meta' or 'data' as it confuses AWS
        self.path = SimpleNamespace(
            # path names for VCF data, stored under self.bucket/self.path.vcf.root/
            vcf=SimpleNamespace(
                root=Path(f'{self.database}/vcf'),              # top level local path
                meta=Path(f'{self.database}/vcf/headers'),      # VCF headers and metadata path
                data=Path(f'{self.database}/vcf/variants'),     # VCF variants path
                files=Path(f'{self.database}/vcf/files'),       # raw VCF files
            ),
            # path names for study data, stored under self.bucket/self.path.study.root/
            study=SimpleNamespace(
                root=Path(f'{self.database}/study'),            # top level local path
                meta=Path(f'{self.database}/study/meta'),       # study metadata
                data=Path(f'{self.database}/study/variants'),   # VCF variants path
                merged=Path(f'{self.database}/study/merged'),   # optional merged path
                frozen=Path(f'{self.database}/study/_frozen'),  # frozen flag, _ files are ignored
                header=Path('_header'),                         # merged header file (relative to current study)
                export=Path('_export'),                         # export prefix (relative to current study)
            ),
            # path names for annotations and metadata, stored under self.bucket/self.path.annotation.root/
            anno=SimpleNamespace(
                root=Path(f'{self.database}/annotations'),             # top level local path
                meta=Path(f'{self.database}/annotations/anno_meta'),   # annotation metadata path
                data=Path(f'{self.database}/annotations/anno'),        # actual annotations path
                files=Path(f'{self.database}/annotations/files'),      # raw VCF files
            ),
            # cache paths, stored under self.bucket/self.path.results.root/
            results=SimpleNamespace(
                root=Path(f'{self.database}/results'),                  # Results root
                stage=Path(f'{self.database}/results/stage'),           # Athena query results stage
                cache=Path(f'{self.database}/results/cache'),           # VDB query cache
                mtime=Path(f'{self.database}/results/mtime'),           # VDB partition modified times
            ),
            ready=Path(f'{self.database}/_ready')                       # VDB is ready to go flag
        )

        # Athena table names
        self.table = SimpleNamespace(
            vcf=SimpleNamespace(
                meta='headers',             # one global headers table
                data='variants',            # one global variants table
            ),
            anno=SimpleNamespace(
                meta='anno_meta',           # one global annotations metadata table
                data='anno',                # one global annotations data table
            ),
            study=SimpleNamespace(
                meta='study_meta',          # study metadata
                data='study_variants',      # study variants
                merged='study_merged',      # optional merged table
            ),
        )

        self.cursor = Connection(
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name=self.aws_region,
            schema_name=self.database,
            s3_staging_dir=f"s3://{self.bucket}/{self.path.results.stage}/",
            cursor_class=AsyncPandasCursor
        ).cursor(max_workers=10)

        # Create Athena tables if needed.
        if allow_db_create:
            self.create_tables()

    # __init__() complete

    # input validation methods
    @staticmethod
    def validate_aid(aid):
        '''
        If aid is a valid UUID, return it lowercased.
        If aid is invalid, raise SystemExit.
        '''
        try:
            return str(uuid.UUID(aid)).lower()
        except (RuntimeError, ValueError):
            raise SystemExit(f"Malformed aid '{aid}'. Must be of the form: {uuid.uuid4()}")

    @staticmethod
    def validate_study_name(study_name):
        ''' Ensure a well-formed study name '''
        if len(study_name) > 64:
            raise SystemExit('Study names must be <= 64 characters')

        if not re.match(r'^[a-zA-Z0-9_]+$', study_name):
            raise SystemExit(f"Study names may only consist of alphanumerics or _: '{study_name}'")

        if study_name in ('meta', 'data'):
            raise SystemExit(f"'{study_name}' is not a valid study name.")

        return study_name

    @staticmethod
    def validate_sample_name(sample_name):
        ''' Ensure a well-formed sample name '''
        if len(sample_name) > 64:
            raise SystemExit('Sample names must be <= 64 characters')

        if not re.match(r'^[a-zA-Z0-9_-]+$', sample_name):
            raise SystemExit(f"Sample names may only consist of alphanumerics, -, or _: '{sample_name}'")

        return sample_name

    @staticmethod
    def get_database_name(db_name):
        ''' Ensure a well-formed database name. If valid, db_name is returned with vdb_ prepended to it. '''

        # Don't accidentally stack vdb_s
        if db_name.lower().startswith('vdb_'):
            db_name = db_name[4:]

        # The actual max db name length is 64 characters, but leave room for vdb_
        if len(db_name) > 60:
            raise SystemExit('Database names must be <= 60 characters')

        if not re.match(r'^[a-z0-9_-]+$', db_name):
            raise SystemExit(f"Database names may only consist of lowercase letters, numbers, - or _: '{db_name}'")

        if db_name.endswith(('-', '_')):
            raise SystemExit(f"Database name must not end with - or _: '{db_name}'")

        return f"vdb_{db_name}"

    @staticmethod
    def validate_wildcard(wildcard):
        ''' Ensure a reasonable wildcard '''
        if len(wildcard) > 64:
            raise SystemExit('Wildcards must be <= 64 characters')

        if not re.match(r'^[a-zA-Z0-9*_-]+$', wildcard):
            raise SystemExit(f"Wildcards may only consist of alphanumerics, -, _, or *: '{wildcard}'")

        return wildcard

    @staticmethod
    def validate_format_field(format_field):
        ''' Ensure a well-formed format field '''
        if len(format_field) > 64:
            raise SystemExit('FORMAT fields must be <= 64 characters')

        if not re.match(r'^[a-zA-Z0-9]+$', format_field):
            raise SystemExit(f"FORMAT fields may only consist of alphanumerics: '{format_field}'")

        return format_field

    def quoted_sample_list(self, items):
        ''' Return a SQL friendly string quoting every sample name in a sequence '''
        return ','.join([f"'{self.validate_sample_name(i)}'" for i in items])

    def quoted_format_list(self, items):
        ''' Return a SQL friendly string quoting every format field in a sequence '''
        return ','.join([f"'{self.validate_format_field(i)}'" for i in items])

    def quoted_aid_list(self, items):
        ''' Return a SQL friendly string quoting all aids in a sequence '''
        return ','.join([f"'{self.validate_aid(i)}'" for i in items])

    @staticmethod
    def scalar(result):
        '''
        For db results that boil down to a single scalar, just return the single value
        If no results, return None
        '''
        if result:
            return result[0][0]
        return None

    def query(self, query, params=None, cache=True, block=True):
        ''' Execute query_with_id and strip the query_id '''
        return self.query_with_id(query, params, cache, block)[1]

    def query_with_id(self, query, params=None, cache=True, block=True):
        '''
        Execute a query and return a tuple of (query_id, all_rows).
        If cache == True, attempt to fetch from the cache, falling back a direct query
        If block == False, immediately return (query_id, future)
        '''
        if cache:
            return fetch_from_cache(self, query, params=params, block=block)

        try:
            query_id, future = self.cursor.execute(query, params)
            if block:
                result = future.result()
                if result.state != 'SUCCEEDED':
                    raise SystemExit(result.state_change_reason)
                return query_id, result.fetchall()
            # async error handling is up to the caller
            return query_id, future
        except Exception as e:
            raise SystemExit(f"Query failed:\n{query}\n{e}")

    @staticmethod
    def collect_futures(futures):
        '''
        Collect all results and pending query IDs.
        Removes completed jobs from futures and requeues throttled requests.
        '''
        results = []
        pending = []
        for i, f in enumerate(futures):
            query_id, future = f
            if future.done():
                futures.pop(i)
                result = future.result()
                if result.state != 'SUCCEEDED':
                    if any(err in result.state_change_reason for err in ['ThrottlingException', 'SlowDown']):
                        sleep(1)
                        pending.append(query_id)
                    else:
                        raise RuntimeError(result.state_change_reason)
                results.append(result.fetchall())
                sleep(0.1)
            else:
                pending.append(query_id)

        return (results, pending)

    def athena_query_status(self, pending):
        '''
        Return the status of a list of Athena jobs as a dict of { 'STATE': [ list of query ids ] }
        '''
        ret = {'QUEUED': [], 'RUNNING': [], 'SUCCEEDED': [], 'FAILED': [], 'CANCELLED': []}
        status = self.athena.batch_get_query_execution(QueryExecutionIds=pending)
        if status and 'QueryExecutions' in status:
            for job in status['QueryExecutions']:
                ret[job['Status']['State']].append(job['QueryExecutionId'])

        return ret

    def parallel_query(self, query, iterate_param, const_params=None, cache=True, parallel=10, complain_time=30):
        '''
        Run up to parallel queries at a time.

        iterate_param is a dict with one key pointing to a list of param values to substitute.
        const_params is optional and are constant for every query.

        No attempt is made to do combinatorial expansion of iterate_param; it should contain
        only one key (eg. 'chrom') pointing to a list of values.

        Raises SystemExit on any query failure. Returns all results and blocks until done.

        If complain_time seconds pass and queries are still queued on Athena, log it.
        '''
        merged_params = []
        for k in iterate_param:
            for v in iterate_param[k]:
                if const_params:
                    merged_params.append({k: v, **const_params})
                else:
                    merged_params.append({k: v})

        results = []
        futures = []
        last_log = datetime.now()

        for params in merged_params:
            try:
                futures.append(self.query_with_id(query, params=params, cache=cache, block=False))
            except Exception as e:
                raise RuntimeError(f"Query failed:\n{query}\n{e}")

            # Don't let too many futures accumulate
            while len(futures) >= parallel:
                new_results, pending = self.collect_futures(futures)
                for r in new_results:
                    results.append(r)
                if pending:
                    if (datetime.now() - last_log).total_seconds() > complain_time:
                        last_log = datetime.now()
                        queued = len(self.athena_query_status(pending)['QUEUED'])
                        if queued:
                            log(f"{queued} job{plural(queued)} queued / {len(pending) - queued} running / {len(results)} completed / {len(merged_params) - len(results)} to go")
                    sleep(1)

        while futures:
            new_results, pending = self.collect_futures(futures)
            for r in new_results:
                results.append(r)
            if pending:
                if (datetime.now() - last_log).total_seconds() > complain_time:
                    last_log = datetime.now()
                    queued = len(self.athena_query_status(pending)['QUEUED'])
                    if queued:
                        log(f"{queued} job{plural(queued)} queued / {len(pending) - queued} running / {len(results)} completed / {len(merged_params) - len(results)} to go")
                sleep(1)

        return results

    def parallel_queries(self, queries, cache=True, parallel=10, complain_time=30):
        '''
        Run up to parallel queries at a time. No parameters are allowed.

        This is useful for running several precomposed queries in parallel.

        Raises SystemExit on any query failure. Returns all results and blocks until done.

        If complain_time seconds pass and queries are still queued on Athena, log it.
        '''
        results = []
        futures = []
        last_log = datetime.now()

        for query in queries:
            try:
                futures.append(self.query_with_id(query, cache=cache, block=False))
            except Exception as e:
                raise RuntimeError(f"Query failed:\n{query}\n{e}")

            # Don't let too many futures accumulate
            while len(futures) >= parallel:
                new_results, pending = self.collect_futures(futures)
                for r in new_results:
                    results.append(r)
                if pending:
                    if (datetime.now() - last_log).total_seconds() > complain_time:
                        last_log = datetime.now()
                        queued = len(self.athena_query_status(pending)['QUEUED'])
                        if queued:
                            log(f"{queued} job{plural(queued)} queued / {len(pending) - queued} running / {len(results)} completed / {len(queries) - len(results)} to go")
                    sleep(1)

        while futures:
            new_results, pending = self.collect_futures(futures)
            for r in new_results:
                results.append(r)
            if pending:
                if (datetime.now() - last_log).total_seconds() > complain_time:
                    last_log = datetime.now()
                    queued = len(self.athena_query_status(pending)['QUEUED'])
                    if queued:
                        log(f"{queued} job{plural(queued)} queued / {len(pending) - queued} running / {len(results)} completed / {len(queries) - len(results)} to go")
                sleep(1)

        return results

    def parallel_query_template(self, query_template, iterate_param, const_params=None, cache=True, parallel=10, complain_time=30):
        '''
        Run up to parallel queries at a time using a template.

        iterate_param is a dict with a single k pointing to a list of values
        to substitute and inject into the query_template using string
        substitution. This is necessary for IN (list...) constructs where
        simple parameter substition can't be used.

        const_params is optional and are constant for every query.

        No attempt is made to do combinatorial expansion of iterate_param; it should contain
        only one key (eg. 'chrom') pointing to a list of values.

        Raises SystemExit on any query failure. Returns all results and blocks until done.

        If complain_time seconds pass and queries are still queued on Athena, log it.
        '''
        queries = []
        for k in iterate_param:
            for v in iterate_param[k]:
                if isinstance(v, str):
                    queries.append(query_template.replace(f"%({k})s", v))
                elif isinstance(v, float):
                    queries.append(query_template.replace(f"%({k})f", v))
                elif isinstance(v, int):
                    queries.append(query_template.replace(f"%({k})d", v))
                else:
                    raise SystemExit(f"parallel_query_template: Unknown type: {v}")

        results = []
        futures = []
        last_log = datetime.now()
        for query in queries:
            try:
                futures.append(self.query_with_id(query, params=const_params, cache=cache, block=False))
            except Exception as e:
                raise RuntimeError(f"Query failed:\n{query}\n{e}")

            # Don't let too many futures accumulate
            while len(futures) >= parallel:
                new_results, pending = self.collect_futures(futures)
                for r in new_results:
                    results.append(r)
                if pending:
                    if (datetime.now() - last_log).total_seconds() > complain_time:
                        last_log = datetime.now()
                        queued = len(self.athena_query_status(pending)['QUEUED'])
                        if queued:
                            log(f"{queued} job{plural(queued)} queued / {len(pending) - queued} running / {len(results)} completed / {len(queries) - len(results)} to go")
                    sleep(1)

        while futures:
            new_results, pending = self.collect_futures(futures)
            for r in new_results:
                results.append(r)
            if pending:
                if (datetime.now() - last_log).total_seconds() > complain_time:
                    last_log = datetime.now()
                    queued = len(self.athena_query_status(pending)['QUEUED'])
                    if queued:
                        log(f"{queued} job{plural(queued)} queued / {len(pending) - queued} running / {len(results)} completed / {len(queries) - len(results)} to go")
                sleep(1)

        return results

    def query_pandas(self, query, params=None, cache=True, block=True):
        ''' Execute query_pandas_with_id and strip the query_id '''
        return self.query_pandas_with_id(query, params, cache, block)[1]

    def query_pandas_with_id(self, query, params=None, cache=True, block=True):
        '''
        Execute a query and return the entire result as a pandas dataframe.
        Note: this bypasses the database fetch and uses the CSV result directly,
        which is faster for larger query results that still fit in memory, but
        introduces the overhead of converting to pandas.

        Fetch from the cache if cache == True.
        '''
        if cache:
            return fetch_from_cache(self, query, params=params, output='pandas', block=block)

        try:
            query_id, future = self.cursor.execute(query, params)
            if not block:
                return query_id, future
            return query_id, future.result().as_pandas()
        except Exception as e:
            raise SystemExit(f"Query failed:\n{query}\n{e}")

    def query_fetch_csv(self, query, out_csv, params=None, cache=True, block=True):
        ''' Execute query_fetch_csv and strip the query_id '''
        return self.query_fetch_csv_with_id(query, out_csv, params, cache, block)[1]

    def query_fetch_csv_with_id(self, query, out_csv, params=None, cache=True, block=True):
        '''
        Execute a query and save the results to csv.
        Fetch from the cache if cache == True.
        '''
        if cache:
            debug(f"cache == {cache}, fetch_from_cache({query}, {params}, {out_csv}, {block})")
            return fetch_from_cache(self, query, params=params, output=out_csv, block=block)

        try:
            debug(f"cursor.execute({query}, {params}")
            query_id, future = self.cursor.execute(query, params)
            if not block:
                return query_id, future
            future.result()
            self.fetch_result_to_csv(query_id, out_csv)

        except Exception as e:
            raise SystemExit(f"Query failed:\n{query}\n{e}")

        return query_id, out_csv

    def fetch_result_to_csv(self, query_id, out_csv):
        ''' Fetch Athena CSV results to a file '''
        debug(f"{self.path.results.stage}/{query_id}.csv to {out_csv}")
        return query_id, self.download_file(f"{self.path.results.stage}/{query_id}.csv", out_csv)

    def upload_to_s3(self, local_file, dest_path):
        ''' Upload a single file to S3 '''
        self.s3.meta.client.upload_file(
            Filename=str(local_file),
            Bucket=self.bucket,
            Key=str(dest_path)
        )

    @staticmethod
    def upload_to_s3_parallel(args):
        ''' Helper to allow parallel upload via Pool.map() '''
        boto3.resource('s3').meta.client.upload_file(
            Filename=args['local_file'],
            Bucket=args['bucket'],
            Key=args['dest_path']
        )

    @staticmethod
    def s3_cp_parallel(args):
        ''' Helper to allow parallel copy via Pool.map() '''
        boto3.resource('s3').meta.client.copy(
            {'Bucket': args['bucket'], 'Key': args['src_path']},
            args['bucket'],
            args['dest_path']
        )

    def sync_to_s3(self, local_path, dest_path=None, parallel=True):
        ''' Recursively copy the contents of a local path to S3 '''
        if not dest_path:
            dest_path = self.path.vcf.root

        skip = len(str(local_path)) + 1

        queue = []
        for (path, _, files) in os.walk(local_path):
            for f in files:
                queue.append({
                    "local_file": str(Path(path) / f),
                    "dest_path": str(Path(dest_path) / path[skip:] / f),
                    "bucket": self.bucket
                })

        if parallel:
            with Pool() as p:
                p.map(connect.upload_to_s3_parallel, queue)
        else:
            for t in queue:
                self.upload_to_s3(t["local_file"], t["dest_path"])

    def download_s3_path(self, s3_path, out_path):
        '''
        Download an arbitrary s3_path from any bucket. Returns the local
        filename saved under out_path.
        '''
        parts = Path(s3_path).parts
        out_file = str(Path(out_path) / parts[-1])
        self.s3.meta.client.download_file(
            Bucket=parts[1],
            Key='/'.join(parts[2:]),
            Filename=out_file
        )
        return out_file

    def download_file(self, prefix, out, full_path=False):
        '''
        Download an arbitrary prefix from the current s3 bucket to out.
        If out is a directory, save to the original filename in that location.
        If full_path is True, save the whole prefix locally, creating directories as needed.
        '''
        prefix = Path(prefix)
        out_path = Path(out)
        if out_path.is_dir():
            if full_path:
                out = out_path / prefix
                out.parent.mkdir(parents=True, exist_ok=True)
            else:
                out = out_path / prefix.parts[-1]

        debug(f"{prefix} -> {out}")

        self.s3.meta.client.download_file(
            Bucket=self.bucket,
            Key=str(prefix),
            Filename=str(out)
        )
        return out

    def download_fileobj(self, prefix, out):
        ''' Download an arbitrary prefix from s3 to a filehandle '''
        self.s3.meta.client.download_fileobj(
            Bucket=self.bucket,
            Key=prefix,
            Fileobj=out
        )
        return out

    def download_gz_fh(self, prefix):
        ''' Download an arbitrary prefix from s3 and return an open filehandle, gzip on the fly '''
        return gzip.GzipFile(fileobj=self.s3.Object(self.bucket, prefix).get()["Body"])

    def download_aid(self, aid, out, dest_path, full_path=False):
        ''' Download an aid from s3 to local '''
        obj = self.ls(str(dest_path), f"aid={aid}")
        if obj:
            for f in obj:
                self.download_file(f, str(out), full_path)
        else:
            raise SystemExit(f"Could not find aid {aid}")

    def ls(self, prefix, filt=None):
        ''' List everything at a prefix with an optional filter '''
        objects = []
        for obj in self.s3.Bucket(self.bucket).objects.filter(Prefix=str(prefix)):
            if filt is None or str(filt) in obj.key:
                objects.append(obj.key)

        return objects

    def s3_path_exists(self, prefix):
        ''' Returns True if any object with the given prefix exists, otherwise False. '''
        return bool(list(self.s3.Bucket(self.bucket).objects.filter(Prefix=prefix).limit(1)))

    def s3_rm_recursive(self, prefix, filt=None):
        ''' Remove all objects matching the prefix '''
        prefix = str(prefix)
        if not prefix or len(prefix) < 3:
            raise SystemExit(f'Refusing s3_rm_recursive() without a proper prefix ({prefix})')

        if filt:
            count = 0
            for obj in self.s3.Bucket(self.bucket).objects.filter(Prefix=str(prefix)):
                if filt in obj.key:
                    obj.delete()
                    count += 1

            return count

        # This method is much faster when the exact prefix is known
        resp = self.s3.Bucket(self.bucket).objects.filter(Prefix=str(prefix)).delete()
        if resp:
            return len(resp[0]['Deleted'])
        return 0

    def s3_cp_recursive(self, src_prefix, dest_prefix):
        ''' Recursively copy all objects from src_prefix to dest_prefix '''
        count = 0
        for obj in self.ls(src_prefix):
            count += 1
            self.s3.meta.client.copy(
                {'Bucket': self.bucket, 'Key': obj},
                self.bucket,
                f"{dest_prefix.rstrip('/')}/{obj[len(src_prefix):]}"
            )
        return count

    def s3_cp_aid(self, aid, dest_prefix):
        '''
        Copy an aid to a new prefix. Similar to s3_cp_recursive() but drops the build=.../ partition.

            vdb_rob/vcf/variants/sample_name=VDB004/build=GRCh37/aid=xxx/yyy.parquet
        ->  vdb_rob/vcf/variants/sample_name=VDB004/aid=xxx/yyy.parquet

        '''
        count = 0
        queue = []
        for obj in self.ls(self.path.vcf.data, f"/aid={aid}"):
            count += 1
            dest_obj = Path(obj)
            src_len = len(self.path.vcf.data.parts)
            queue.append(
                {
                    'bucket': self.bucket,
                    'src_path': obj,
                    'dest_path': f"{dest_prefix.rstrip('/')}/{'/'.join(dest_obj.parts[src_len:src_len+1] + dest_obj.parts[src_len+2:])}"
                }
            )

        with Pool() as p:
            p.map(connect.s3_cp_parallel, queue)

        return count

    def delete_aid(self, aids, aid_type="vcf"):
        ''' Delete an aid and drop relevant database partitions '''

        if aid_type == "vcf":
            s3_path = self.path.vcf.root
            tables = (self.table.vcf.meta, self.table.vcf.data)
        elif aid_type == "anno":
            s3_path = self.path.anno.root
            tables = (self.table.anno.meta, self.table.anno.data)
        else:
            raise SystemExit(f"Unknown aid_type: {aid_type}")

        if isinstance(aids, str):
            aids = [aids]

        for aid in aids:
            self.validate_aid(aid)

        log(f"Deleting {aid_type} data")
        for aid in aids:
            if not self.s3_rm_recursive(s3_path, f"/aid={aid}"):
                if len(aids) == 1:
                    raise SystemExit(f"No aid found in {aid_type}: {aid}")
                # If bulk deleting, just complain
                log(f"No such aid: {aid}")
                continue

        debug(f"Dropping partitions")
        for table in tables:
            update_table_mtime(self, table)
            self.parallel_query(
                f"ALTER TABLE `{table}` DROP IF EXISTS PARTITION (aid=%(aid)s);",
                iterate_param={'aid': aids}
            )

    @staticmethod
    def get_study_path(study_name, base):
        ''' Return the correct s3 path for a given study and base path '''
        return str(base / Path(f'study_name={study_name}'))

    def create_tables(self):
        ''' Create all necessary VDB tables '''
        db_exists = self.database in [d[0] for d in self.query("SHOW DATABASES;")]

        if db_exists:
            if self.ls(self.path.ready):
                return

        if not confirm(f"VDB database {self.database} does not exist in AWS region {self.aws_region}. Create it? n/Y: ", default=True):
            raise SystemExit("Aborted.")

        log(f"Initializing new VDB '{self.database}' at s3://{self.bucket}/{self.database}/")

        self.query(f"CREATE DATABASE IF NOT EXISTS `{self.database}`;")

        for table, s3path in (
                (self.table.vcf.data, self.path.vcf.data),
                (self.table.vcf.meta, self.path.vcf.meta),
                (self.table.study.data, self.path.study.data),
                (self.table.study.meta, self.path.study.meta),
                (self.table.study.merged, self.path.study.merged),
                (self.table.anno.data, self.path.anno.data),
                (self.table.anno.meta, self.path.anno.meta)
            ):

            update_table_mtime(self, table)
            self.query(
                create_table_sql(table),
                params={
                    "location": f"s3://{self.bucket}/{s3path}/"
                }
            )

        self.s3.Object(self.bucket, str(self.path.ready)).put(
            Body=json.dumps({"tables_created_on": timestamp()})
        )

    def get_annotation_query(self, study_name, anno):
        '''
        Return a query suffix for annotations, or None if no annotation is requested.
        TODO: expand this to include anno version and optional aid
        '''
        if not anno:
            return None

        aid = self.query(
            f"""
                SELECT am.aid
                FROM {self.table.anno.meta} AS am, {self.table.study.meta} AS sm
                WHERE
                    am.anno_name = %(anno)s
                    AND sm.key = 'build'
                    AND am.build = sm.value
                    AND study_name = %(study_name)s
                ;
            """,
            params={"study_name": study_name, "anno": anno}
        )
        if not aid:
            raise SystemExit(f"There is no annotation named {anno} with a matching reference build.")
        if len(aid) > 1:
            raise SystemExit(f"There are multiple matching annotations for {anno}. Please specify a version or aid.")

        if anno == "Ensembl":
            return f"a.aid = '{self.scalar(aid)}' AND a.feature = 'gene'"

        return f"a.aid = '{self.scalar(aid)}'"

    def get_current_study_checkpoint(self, study_name):
        ''' Get the most recent checkpoint for this study '''
        self.assert_study_exists(study_name)

        checkpoint = self.scalar(self.query(
            f"SELECT max(checkpoint) FROM {self.table.study.data} WHERE study_name = %(study_name)s ;",
            params={"study_name": study_name}
            ))
        if pd.isna(checkpoint):
            return 0
        return checkpoint

    def get_study_chroms(self, study_name, checkpoint):
        ''' Fetch all chroms in a study at a given checkpoint '''
        return [chrom[0] for chrom in self.query(
            f"""
                SELECT DISTINCT(chrom)
                FROM {self.table.study.data}
                WHERE study_name = %(study_name)s
                AND checkpoint = %(checkpoint)d
                ORDER BY chrom
                ;
            """,
            params={"study_name": study_name, "checkpoint": checkpoint}
        )]

    @staticmethod
    def get_merge_partition(study_name, checkpoint):
        ''' Return the partition to be used for merged variants '''
        return f"study_name={study_name}/checkpoint={checkpoint}"

    def merge_study(self, study_name, force_merge=False, anno_name=None, checkpoint=None, square_off=None, format_fields=None): # pylint: disable=too-many-statements
        ''' Merge all samples in this study at the given checkpoint '''
        self.assert_study_is_unfrozen(study_name)

        if checkpoint is None:
            checkpoint = self.get_current_study_checkpoint(study_name)
        elif checkpoint > self.get_current_study_checkpoint(study_name):
            raise SystemExit(f"Requested checkpoint {checkpoint} does not exist in study {study_name}.")

        count = self.scalar(
            self.query(
                f"SELECT count(*) FROM {self.table.study.data} WHERE study_name = %(study_name)s AND checkpoint = %(checkpoint)d ;",
                params={"study_name": study_name, "checkpoint": checkpoint}
            )
        )
        if not count:
            raise SystemExit(f"No variants found in study {study_name}")

        if format_fields:
            # unique fields only
            format_fields = set(format_fields)
            # GT is mandatory
            format_fields.add('GT')

        partition = self.get_merge_partition(study_name, checkpoint)

        if square_off:
            if format_fields:
                samples = f"""
                    transform_values(
                        map(array['{square_off}'], array[element_at(m.samples, '{square_off}')]),
                        (k, v) -> map_filter(v, (k, v) -> k in ({self.quoted_format_list(format_fields)}) )
                    )
                    """
            else:
                samples = f"map(array['{square_off}'], array[element_at(m.samples, '{square_off}')])"
            export_location = f"{self.path.study.merged}/{self.path.study.export}/{partition}/_sample_name={square_off}/"
        else:
            if format_fields:
                samples = f"transform_values(m.samples, (k, v) -> map_filter(v, (k, v) -> k in ({self.quoted_format_list(format_fields)})))"
            else:
                samples = "m.samples"
            export_location = f"{self.path.study.merged}/{self.path.study.export}/{partition}/_sample_name=ALL_SAMPLES/"

        # annotations are exported separately
        annotation_query = self.get_annotation_query(study_name, anno_name)
        if annotation_query:
            export_location = f"{export_location}_anno={anno_name}/"

        # Only need to merge once for each checkpoint, since they can't be changed.
        if force_merge or not self.s3_path_exists(f"{self.path.study.merged}/{partition}"):
            self.s3_rm_recursive(f"{self.path.study.merged}/{partition}/")
            self.s3_rm_recursive(export_location, partition)

            log(f"Merging variants for checkpoint {checkpoint}")
            self.merge_study_variants(study_name, checkpoint)

        header_path = f"{self.path.study.merged}/{self.path.study.export}/{partition}/{self.path.study.header}"
        if force_merge or not self.s3_path_exists(header_path):
            log("Merging headers")
            self.merge_study_headers(study_name, checkpoint, header_path)

        # Reuse the existing export if possible
        if self.s3_path_exists(export_location):
            log("No updates since the previous export, reusing existing merge.")
            return (header_path, export_location)

        chroms = self.get_study_chroms(study_name, checkpoint)
        debug(chroms)

        # Create a temporary TSV table
        table_id = f"{study_name}_merge_{str(hex(int(time()*10000000))[8:])}"
        self.query(f"DROP TABLE IF EXISTS {table_id};")
        self.query(
            f"""
                CREATE EXTERNAL TABLE {table_id} (
                    `pos` bigint,
                    `varid` string,
                    `ref` string,
                    `alt` string,
                    `qual` float,
                    `filt` string,
                    `info` string,
                    `samples` string
                )
                    PARTITIONED BY (`chrom` STRING)
                    ROW FORMAT DELIMITED
                    FIELDS TERMINATED BY '\\t'
                    ESCAPED BY '\\\\'
                    LINES TERMINATED BY '\\n'
                    NULL DEFINED AS ''
                    LOCATION %(location)s
                ;
            """,
            params={
                "location": f"s3://{self.bucket}/{export_location}",
            }
        )

        if annotation_query:
            log(f"Annotating variants with {anno_name}")

            # partition on chrom, samples are JSON
            self.parallel_query(
                f"""
                    INSERT INTO "{table_id}"
                        (chrom, pos, varid, ref, alt, qual, filt, info, samples)
                    WITH annospan AS
                        (
                            SELECT anno.*, spanpos
                            FROM anno
                            CROSS JOIN unnest(sequence(pos/100, varend/100)) span(spanpos)
                            WHERE pos <= 999999999
                                    AND varend >= 0
                                    AND chrom = %(chrom)s
                        )
                    SELECT
                         m.chrom AS chrom
                        ,m.pos + 1 AS pos
                        ,array_join(array_agg(a.varid),';') AS varid
                        ,m.ref AS ref
                        ,m.alt AS alt
                        ,arbitrary(m.qual) AS qual
                        ,arbitrary(array_join(m.filters, ';')) AS filt
                        ,arbitrary(array_join(zip_with(map_keys(m.infos), map_values(m.infos), (v1, v2) -> concat(v1, '=', v2)), ';')) AS info
                        ,arbitrary(json_format(cast({samples} AS JSON))) AS samples
                    FROM
                        {self.table.study.merged} as m
                    LEFT JOIN annospan AS a
                        ON {annotation_query}
                        AND m.chrom = a.chrom
                        AND m.pos >= a.pos
                        AND m.varend <= a.varend
                        AND m.pos / 100 = a.spanpos
                    WHERE
                        m.chrom = %(chrom)s
                        AND study_name = %(study_name)s
                        AND checkpoint = %(checkpoint)d
                        AND m.pos BETWEEN 0 AND 999999999
                    GROUP BY
                        m.chrom, m.pos, m.ref, m.alt
                    ;
                """,
                iterate_param={"chrom": chroms},
                const_params={"study_name": study_name, "checkpoint": checkpoint},
            )
        else:
            # partition on chrom, samples are JSON
            log(f"Writing variants")
            self.parallel_query(
                f"""
                    INSERT INTO "{table_id}"
                      (chrom, pos, ref, alt, qual, filt, info, samples)
                    SELECT
                      m.chrom,
                      m.pos + 1,
                      m.ref,
                      m.alt,
                      m.qual,
                      array_join(m.filters, ';'),
                      array_join(zip_with(map_keys(m.infos), map_values(m.infos), (v1, v2) -> concat(v1, '=', v2)), ';'),
                      json_format(CAST({samples} AS json))
                    FROM
                      {self.table.study.merged} m
                    WHERE
                      chrom = %(chrom)s
                      AND study_name = %(study_name)s
                      AND checkpoint = %(checkpoint)d
                    ;
                """,
                iterate_param={"chrom": chroms},
                const_params={"study_name": study_name, "checkpoint": checkpoint}
            )

        # The table only exists to generate the TSV, so drop it.
        self.query(f"DROP TABLE IF EXISTS {table_id};")

        return (header_path, export_location)

    def merge_study_headers(self, study_name, checkpoint, header_path):
        '''
        Merge VCF headers in the given study and save to header_path.
        Returns the merged header up to (but not including) the sample column names.
        '''
        study_headers = self.query(
            f"""
                SELECT refhash, header FROM {self.table.vcf.meta} WHERE aid IN (
                    SELECT DISTINCT aid
                    FROM {self.table.study.data}
                    WHERE study_name = %(study_name)s
                    AND checkpoint = %(checkpoint)d
                )
            ;
            """,
            params={
                "study_name": study_name, "checkpoint": checkpoint
            }
        )

        contigs = []
        for line in study_headers[0][1].splitlines():
            if line.startswith('##contig='):
                contigs.append(line)

        headers = []
        the_date = datetime.today()
        headers.append(f'##fileformat=VCFv4.1')
        headers.append(f'##fileDate={the_date.year}{the_date.month:02}{the_date.day:02}')
        headers.append(f'''##source="Spiral Genetics VDB",description="biograph vdb study export {study_name} --checkpoint {checkpoint}"''')

        study_meta = self.get_metadata_from_study(study_name)
        for chkpt in [c for c in sorted(study_meta) if c.startswith('checkpoint ')]:
            _, rev = chkpt.split()
            if int(rev) > checkpoint:
                continue
            headers.append(f'''##checkpoint="{chkpt}: {study_meta[chkpt]}"''')

        # contigs must appear in the original order
        for contig in contigs:
            headers.append(contig)

        # Merge all INFO, FORMAT, FILTER, and ALT lines
        old_headers = set()
        for header in study_headers:
            for line in header[1].splitlines():
                if line.startswith(('##source=', '##contig=', '##fileDate=', '##fileformat=', '#CHROM')):
                    continue
                old_headers.add(line)

        # Add in computed fields
        old_headers.add('##INFO=<ID=N_MISS,Number=1,Type=Integer,Description="Number of samples missing this variant">')
        old_headers.add('##INFO=<ID=F_MISS,Number=1,Type=Float,Description="Fraction of samples missing this variant">')

        for line in sorted(old_headers):
            headers.append(line)

        headers.append(
            f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        )

        header = '\n'.join(headers)
        self.s3.Object(self.bucket, header_path).put(
            Body=header
        )
        return header

    def merge_study_variants(self, study_name, checkpoint):
        ''' Merge variants in a study. This only needs to be done once per checkpoint. '''
        partition = self.get_merge_partition(study_name, checkpoint)
        self.s3_rm_recursive(f"{self.table.study.merged}/{partition}")
        update_table_mtime(self, self.table.study.merged, partition=partition)

        self.query(f"ALTER TABLE {self.table.study.merged} DROP IF EXISTS PARTITION (study_name='{study_name}', checkpoint={checkpoint});")
        study_variants_mtime = get_table_mtime(self, self.table.study.data, auto_update=True, partition=partition)

        # INFO['N_MISS'] == number of individuals missing this variant
        # INFO['F_MISS'] == N_MISS / samples_in_study
        samples_in_study = self.scalar(self.query(
            f"""
                SELECT COUNT(DISTINCT sample_name)
                FROM {self.table.study.data}
                WHERE study_name = %(study_name)s
                AND checkpoint = %(checkpoint)d
                ;
            """,
            params={"study_name": study_name, "checkpoint": checkpoint}
        ))

        chroms = self.get_study_chroms(study_name, checkpoint)
        update_table_mtime(self, self.table.study.merged, partition=partition, ts=study_variants_mtime)
        self.query(f"""ALTER TABLE {self.table.study.merged} ADD IF NOT EXISTS PARTITION (study_name='{study_name}', checkpoint={checkpoint});""")

        #   QUAL is the maximum value for any sample
        #   FILT is a sorted array of distinct filter entries from all samples
        #   INFO is an aggregate of all unique info map entries from all samples
        #     with the proper NS count and N_MISS + F_MISS added
        #   FORMAT does not exist; it is later derived from SAMPLES
        #   SAMPLES is a map of sample names to all sample fields
        #
        # Counting NS as count(DISTINCT(sample_name)) instead of count(sample)
        # is necessary since multiple calls at the same site in the same
        # individual would be counted multiple times, but subsequently collapsed
        # into a single sample entry.
        the_query = f"""
            INSERT INTO "{self.table.study.merged}"
                (spans, reflen, chrom, pos, varend, varid, ref, alt, qual, filters, infos, samples, study_name, checkpoint)
            SELECT
                arbitrary(spans),
                arbitrary(reflen),
                arbitrary(chrom),
                arbitrary(pos),
                arbitrary(varend),
                arbitrary(varid),
                arbitrary(ref),
                arbitrary(alt),
                max(qual),
                array_sort(array_distinct(array_agg(filt))),
                map_concat(
                    map_union(map_filter(info, (k, v) -> k != 'NS')),
                    map(ARRAY['NS'], ARRAY[CAST(count(DISTINCT(sample_name)) AS VARCHAR)]),
                    map(
                        ARRAY['N_MISS', 'F_MISS'],
                        ARRAY[
                            CAST({samples_in_study} - count(DISTINCT(sample_name)) AS VARCHAR),
                            CAST(round(({samples_in_study} - count(DISTINCT(sample_name))) / {samples_in_study}.0, 5) AS VARCHAR)
                        ]
                    )
                ),
                map_agg(sample_name, sample),
                %(study_name)s,
                %(checkpoint)d
            FROM
                {self.table.study.data}
            WHERE
                chrom = %(chrom)s
                AND study_name = %(study_name)s
                AND checkpoint = %(checkpoint)d
            GROUP BY
                chrom, pos, ref, alt
        """
        # Run the chroms in parallel
        self.parallel_query(
            the_query,
            iterate_param={"chrom": chroms},
            const_params={"study_name": study_name, "checkpoint": checkpoint},
        )
        update_table_mtime(self, self.table.study.merged, partition=partition, ts=study_variants_mtime)

    def study_freeze(self, study_name):
        ''' Freeze a study '''
        self.assert_study_exists(study_name)
        self.s3.Object(self.bucket, f"{self.get_study_path(study_name, self.path.study.meta)}/{self.path.study.frozen}").put(
            Body=json.dumps({"frozen_on": timestamp()})
        )

    def study_unfreeze(self, study_name):
        ''' Unfreeze a study '''
        self.assert_study_exists(study_name)
        self.s3.Object(self.bucket, f"{self.get_study_path(study_name, self.path.study.meta)}/{self.path.study.frozen}").delete()

    def study_is_frozen(self, study_name):
        ''' Check if a study is frozen '''
        self.assert_study_exists(study_name)
        return self.s3_path_exists(f"{self.get_study_path(study_name, self.path.study.meta)}/{self.path.study.frozen}")

    def assert_study_is_unfrozen(self, study_name):
        ''' Raise if the current study is frozen '''
        if self.study_is_frozen(study_name):
            raise SystemExit(f"The study '{study_name}' is frozen and cannot be altered.")

    def study_exists(self, study_name):
        ''' Returns True if study exists, otherwise False. '''
        return self.s3_path_exists(f"{self.get_study_path(self.validate_study_name(study_name), self.path.study.meta)}/")

    def assert_study_exists(self, study_name, must_exist=True):
        ''' Check if a study exists. If must_exist is True or False, raise if the expectation is not met. '''
        exists = self.study_exists(study_name)

        if must_exist == exists:
            return exists

        if must_exist:
            raise SystemExit(f"No such study '{study_name}'.")

        raise SystemExit(f"Study '{study_name}' already exists.")

    def create_study(self, study_name):
        ''' Returns True if study exists, otherwise False. '''
        self.assert_study_exists(study_name, must_exist=False)
        update_table_mtime(self, self.table.study.meta)
        self.add_metadata_to_study(study_name, 'created_on', timestamp(), 'timestamp')

    def add_metadata_to_study(self, study_name, key, value, dtype='str'):
        ''' Create a new tiny tsv for this metadata entry '''
        partition = f"study_name={study_name}"
        update_table_mtime(self, self.table.study.meta, partition=partition)
        self.s3.Object(self.bucket, f"{self.get_study_path(study_name, self.path.study.meta)}/{key}").put(
            Body=f'{key}\t{value}\t{dtype}\n'
        )
        self.query(f"ALTER TABLE `{self.table.study.meta}` ADD IF NOT EXISTS PARTITION (study_name='{study_name}');")

    def remove_metadata_from_study(self, study_name, key):
        ''' Delete a metadata entry '''
        partition = f"study_name={study_name}"
        update_table_mtime(self, self.table.study.meta, partition=partition)
        self.s3.Object(self.bucket, f"{self.get_study_path(study_name, self.path.study.meta)}/{key}").delete()

    def get_metadata_from_study(self, study_name, key=None):
        '''
        Get all metadata from a study as a dict.
        If key is specified, return only the value for that key.
        '''
        if key:
            debug(f"query {study_name} for {key}")
            row = self.query(
                f"""
                    SELECT value, dtype
                    FROM {self.table.study.meta}
                    WHERE study_name = %(study_name)s
                    AND key = %(key)s
                    ;
                """,
                params={"study_name": study_name, "key": key},
            )
            if not row:
                debug(f"{key} not found")
                return {}

            value, dtype = row[0]
            debug(f"{typed(dtype, value)}")
            return typed(dtype, value)

        ret = {}
        debug(f"query {study_name}")
        for row in self.query(
                f"""
                    SELECT key, value, dtype
                    FROM {self.table.study.meta}
                    WHERE study_name = %(study_name)s
                    ;
                """,
                params={"study_name": study_name},
            ):
            key, value, dtype = row
            ret[key] = typed(dtype, value)

        debug(ret)
        return ret

    def get_study_sample_names(self, study_name, checkpoint=None):
        '''
        Return a sorted list of all sample names from the latest checkpoint for this study
        '''
        if checkpoint is None:
            checkpoint = self.get_current_study_checkpoint(study_name)

        return [
            sample_name[0] for sample_name in self.query(
                f"""
                    SELECT DISTINCT sample_name
                    FROM {self.table.study.data}
                    WHERE study_name = %(study_name)s
                    AND checkpoint = %(checkpoint)d
                    ORDER BY sample_name ASC
                    ;
                """,
                params={"study_name": study_name, "checkpoint": checkpoint}
            )
        ]

    def delete_study(self, study_name, checkpoint=None):
        '''
        Remove all objects under the given study prefix and drop all associated tables.

        If checkpoint is given, only delete the given checkpoint.

        Returns the number of objects deleted from S3.
        '''
        self.assert_study_is_unfrozen(study_name)

        if checkpoint:
            partition = f"study_name={study_name}/checkpoint={checkpoint}"
            if not (self.get_metadata_from_study(study_name, f"checkpoint {checkpoint}") or self.s3_path_exists(f"{self.path.study.data}/{partition}")):
                raise SystemExit(f"No such checkpoint {checkpoint} for study {study_name}")
            log(f"Removing S3 data for '{study_name}' at checkpoint {checkpoint}")
        else:
            partition = f"study_name={study_name}"
            log(f"Removing S3 data for '{study_name}'")

        # NOTE: trailing / is important here so we don't delete studies with the same name prefix
        count = self.s3_rm_recursive(self.path.study.root, f"/{partition}/")

        # drop the tables
        if checkpoint:
            log(f"Dropping SQL partitions for '{study_name}' at checkpoint {checkpoint}")
            for table in (self.table.study.meta, self.table.study.data, self.table.study.merged):
                update_table_mtime(self, table)
                update_table_mtime(self, table, partition=f"study_name={study_name}")
                clear_table_mtime(self, table, partition=partition)
            for table in (self.table.study.data, self.table.study.merged):
                self.query(
                    f"ALTER TABLE `{table}` DROP IF EXISTS PARTITION (study_name=%(study_name)s, checkpoint=%(checkpoint)d);",
                    params={"study_name": study_name, "checkpoint": checkpoint}
                )
            self.remove_metadata_from_study(study_name, f"checkpoint {checkpoint}")
        else:
            log(f"Dropping SQL partitions for '{study_name}'")
            for table in (self.table.study.meta, self.table.study.data, self.table.study.merged):
                update_table_mtime(self, table)
                update_table_mtime(self, table, partition=partition)
                self.query(f"ALTER TABLE `{table}` DROP IF EXISTS PARTITION (study_name='{study_name}');")

        return count

    def find_all_matching_vcf_aids(self, items):
        '''
        Turn a list of potentially redundant aids, sample names, and wildcard
        globs into a sorted list of (sample_name, [aids]) tuples for matching VCFs.
        '''
        debug(items)
        if not isinstance(items, list):
            items = [items]

        aids = set()
        samples = set()
        wildcards = set()
        for item in items:
            try:
                aids.add(self.validate_aid(item))
            except (SystemExit, ValueError):
                if '*' in item:
                    wildcards.add(self.validate_wildcard(item))
                else:
                    samples.add(self.validate_sample_name(item))

        if samples:
            for aid in self.query(f"SELECT aid FROM {self.table.vcf.meta} WHERE sample_name IN ({self.quoted_sample_list(samples)}) ;"):
                aids.add(aid[0])

        if wildcards:
            reply = self.parallel_query(
                f"SELECT aid FROM {self.table.vcf.meta} WHERE sample_name LIKE %(wildcard)s ;",
                iterate_param={"wildcard": [wc.replace('*', '%') for wc in wildcards]}
            )
            for r in reply:
                for aid in r:
                    if not aid:
                        continue
                    aids.add(aid[0])

        debug("samples:", samples)
        debug("wildcards:", wildcards)
        debug("aids:", aids)

        if not aids:
            return None

        return self.query(
            f"""
                SELECT sample_name, cast(array_agg(DISTINCT(aid)) AS JSON)
                FROM {self.table.vcf.meta}
                WHERE aid IN ({self.quoted_aid_list(aids)})
                GROUP BY sample_name
                ORDER BY sample_name ASC
                ;
            """)

    def find_all_matching_study_aids(self, study_name, checkpoint, items):
        '''
        Turn a list of potentially redundant aids, sample names, and wildcard
        globs into a sorted list of (sample_name, [aids]) tuples for matching
        study variants.
        '''
        debug(items)
        if not isinstance(items, list):
            items = [items]

        aids = set()
        samples = set()
        wildcards = set()
        for item in items:
            try:
                aids.add(self.validate_aid(item))
            except (SystemExit, ValueError):
                if '*' in item:
                    wildcards.add(self.validate_wildcard(item))
                else:
                    samples.add(self.validate_sample_name(item))

        if samples:
            for aid in self.query(
                    f"""
                        SELECT DISTINCT(aid)
                        FROM {self.table.study.data}
                        WHERE study_name = %(study_name)s
                        AND checkpoint = %(checkpoint)d
                        AND sample_name IN ({self.quoted_sample_list(samples)})
                        ;
                    """,
                    params={"study_name": study_name, "checkpoint": checkpoint}):
                aids.add(aid[0])

        if wildcards:
            reply = self.parallel_query(
                f"""
                    SELECT DISTINCT(aid)
                    FROM {self.table.study.data}
                    WHERE study_name = %(study_name)s
                    AND checkpoint = %(checkpoint)d
                    AND sample_name LIKE %(wildcard)s ;
                """,
                iterate_param={"wildcard": [wc.replace('*', '%') for wc in wildcards]},
                const_params={"study_name": study_name, "checkpoint": checkpoint}
            )
            for r in reply:
                for aid in r:
                    if not aid:
                        continue
                    aids.add(aid[0])

        debug("samples:", samples)
        debug("wildcards:", wildcards)
        debug("aids:", aids)

        if not aids:
            return None

        return self.query(
            f"""
                SELECT sample_name, cast(array_agg(DISTINCT(aid)) AS JSON)
                FROM {self.table.study.data}
                WHERE study_name = %(study_name)s
                AND aid IN ({self.quoted_aid_list(aids)})
                GROUP BY sample_name
                ORDER BY sample_name ASC
                ;
            """,
            params={"study_name": study_name})

    def check_matching_refnames(self, aids):
        '''
        Ensure that all aids use the same refname.
        Returns the refname, or raises if more than one is found.
        '''
        reply = self.query(
            f"SELECT aid, refname FROM {self.table.vcf.meta} WHERE aid IN ({self.quoted_aid_list(aids)});",
        )

        # [(f6f8fee8-c488-4154-9286-988c124498e2, grch38), da5279e9-9ccc-4a40-bf3b-40fad12b98c1, hs37d5)]
        if len({r[1] for r in reply}) > 1:
            log(f"All variants in a study must be called against the same reference:")
            for ref in reply:
                log(f"  {ref[0]} -> {ref[1]}")
            raise SystemExit("\nAborted.")

        return reply[0][1]

    def add_to_study(self, study_name, items):
        '''
        Add VCF variants to a study
        '''
        self.assert_study_is_unfrozen(study_name)

        samples = self.find_all_matching_vcf_aids(items)
        if not samples:
            raise SystemExit(f"No matching VCFs found.")
        debug(samples)

        # [('HG002', ['6bc81054-d0e6-4b2a-92d1-82dd2950a33d']),
        #  ('HG003',
        #   ['cabceac9-932f-42aa-8043-7a93dacc13bf',
        #    '9649df68-5b62-4993-8e67-89de88694606']),
        #  ('HG004', ['59890231-b506-4098-acb0-c51bf270536d'])]
        #
        # -> {
        #      '6bc81054-d0e6-4b2a-92d1-82dd2950a33d': 'HG002',
        #      'cabceac9-932f-42aa-8043-7a93dacc13bf': 'HG003',
        #      '9649df68-5b62-4993-8e67-89de88694606': 'HG003',
        #      '59890231-b506-4098-acb0-c51bf270536d': 'HG004'
        #    }
        aids = {aid:sample for sample, aids in samples for aid in aids}
        debug(aids)

        log(f"Matching VCFs:")
        for sample in samples:
            log(f"  {sample[0]}: {', '.join(sample[1])}")
        log("")

        vcf_ref = self.check_matching_refnames(aids.keys())
        vcf_build = refhash(lookup=vcf_ref).build()

        study_ref = self.get_metadata_from_study(study_name, 'refname')

        if not study_ref:
            debug(f"{study_name} has no ref, set to {vcf_ref}")
            self.add_metadata_to_study(study_name, 'refname', vcf_ref)
            self.add_metadata_to_study(study_name, 'build', vcf_build)
        elif study_ref != vcf_ref:
            raise SystemExit(f"Study {study_name} uses reference {study_ref}, but the specified VCFs use {vcf_ref}.")

        checkpoint = self.get_current_study_checkpoint(study_name)

        in_study = {a[0] for a in self.query(
            f"""
                SELECT DISTINCT aid
                FROM {self.table.study.data}
                WHERE study_name = %(study_name)s
                AND checkpoint = {checkpoint}
            ;
            """,
            params={"study_name": study_name}
        )}

        if set(aids.keys()).intersection(in_study):
            log(f"The following VCFs are already in this study at checkpoint {checkpoint} and will be skipped:")
            for aid in set(aids.keys()).intersection(in_study):
                log(f"  {aids[aid]}: {aid}")
                aids.pop(aid)
            log("")

        if not aids:
            raise SystemExit("Nothing left to import.")

        count = self.scalar(
            self.query(
                f"SELECT count(*) FROM {self.table.vcf.data} WHERE aid IN ({self.quoted_aid_list(aids.keys())});"
            )
        )
        if not count:
            raise SystemExit(f"No variants found.")

        log(f"Adding {count:,} variants from {len(aids)} VCF{plural(len(aids))} to study {study_name}")
        new_checkpoint = checkpoint + 1
        update_table_mtime(self, self.table.study.data, partition=f"study_name={study_name}")
        # carry previous rev forward
        # TODO: make this a parallel query
        if checkpoint:
            self.query(
                f"""
                    INSERT INTO "{self.table.study.data}"
                        (spans, reflen, chrom, pos, varend, varid, ref, alt, qual, filt, info, sample, study_name, checkpoint, sample_name, aid)
                    SELECT
                        spans,
                        reflen,
                        chrom,
                        pos,
                        varend,
                        varid,
                        ref,
                        alt,
                        qual,
                        filt,
                        info,
                        sample,
                        %(study_name)s,
                        %(new_checkpoint)d,
                        sample_name,
                        aid
                    FROM {self.table.study.data}
                    WHERE study_name = %(study_name)s
                    AND checkpoint = %(old_checkpoint)d
                    ;
                """,
                params={"study_name": study_name, "old_checkpoint": checkpoint, "new_checkpoint": new_checkpoint}
            )

        # Create partitions pointing to the variants. This saves the time and cost of making a copy.
        queries = []
        parts = []
        for aid, sample_name in aids.items():
            parts.append(
                f"""(study_name='{study_name}', checkpoint={new_checkpoint}, sample_name='{sample_name}', aid='{aid}') LOCATION 's3://{self.bucket}/{self.path.vcf.data}/sample_name={sample_name}/build={vcf_build}/aid={aid}/'""")
            # Keep the query to a reasonable size
            if len(parts) > 100:
                queries.append(
                    f"""ALTER TABLE {self.table.study.data} ADD IF NOT EXISTS PARTITION {' PARTITION '.join(parts)};""")
                parts = []

        if parts:
            queries.append(
                f"""ALTER TABLE {self.table.study.data} ADD IF NOT EXISTS PARTITION {' PARTITION '.join(parts)};""")
            parts = []

        # Add partitions in parallel
        self.parallel_queries(queries)
        # Since we no longer copy data, create a placeholder for the path on s3
        self.s3.Object(self.bucket, f"{self.get_study_path(study_name, self.path.study.data)}/checkpoint={new_checkpoint}/_study_add").put(
            Body=json.dumps(aids)
        )
        # Add metadata
        self.add_metadata_to_study(study_name, f"checkpoint {new_checkpoint}", f"added {'; '.join([f'{v}: {k}' for k, v in aids.items()])}")

    def copy_from_study(self, src_study, src_checkpoint, dest_study, items):
        '''
        Copy variants from one study to another
        '''
        self.assert_study_is_unfrozen(dest_study)
        self.assert_study_exists(src_study)

        if src_checkpoint:
            if src_checkpoint > self.get_current_study_checkpoint(src_study):
                raise SystemExit(f"No such checkpoint {src_checkpoint} for study {src_study}")
        else:
            src_checkpoint = self.get_current_study_checkpoint(src_study)

        dest_ref = self.get_metadata_from_study(dest_study, 'refname')
        src_ref = self.get_metadata_from_study(src_study, 'refname')

        if not all([src_checkpoint, src_ref]):
            raise SystemExit(f"Study {src_study} has no variants.")

        if not dest_ref:
            debug(f"{dest_study} has no ref, set to {src_ref}")
            self.add_metadata_to_study(dest_study, 'refname', src_ref)
            self.add_metadata_to_study(dest_study, 'build', refhash(lookup=src_ref).build())

        if dest_ref and src_ref != dest_ref:
            raise SystemExit(f"Studies use different references ({src_study}:{src_ref} vs. {dest_study}:{dest_ref})")

        samples = self.find_all_matching_study_aids(src_study, src_checkpoint, items)
        if not samples:
            raise SystemExit(f"No matching VCFs found.")
        debug(samples)

        # [('HG002', ['6bc81054-d0e6-4b2a-92d1-82dd2950a33d']),
        #  ('HG003',
        #   ['cabceac9-932f-42aa-8043-7a93dacc13bf',
        #    '9649df68-5b62-4993-8e67-89de88694606']),
        #  ('HG004', ['59890231-b506-4098-acb0-c51bf270536d'])]
        #
        # -> {
        #      '6bc81054-d0e6-4b2a-92d1-82dd2950a33d': 'HG002',
        #      'cabceac9-932f-42aa-8043-7a93dacc13bf': 'HG003',
        #      '9649df68-5b62-4993-8e67-89de88694606': 'HG003',
        #      '59890231-b506-4098-acb0-c51bf270536d': 'HG004'
        #    }
        aids = {aid:sample for sample, aids in samples for aid in aids}
        debug(aids)

        log(f"Matching variants from {src_study}:{src_checkpoint}")
        for sample in samples:
            log(f"  {sample[0]}: {', '.join(sample[1])}")
        log("")

        dest_checkpoint = self.get_current_study_checkpoint(dest_study)
        if dest_checkpoint:
            in_study = {a[0] for a in self.query(
                f"""
                    SELECT DISTINCT aid
                    FROM {self.table.study.data}
                    WHERE study_name = %(study_name)s
                    AND checkpoint = {dest_checkpoint}
                ;
                """,
                params={"study_name": dest_study}
            )}

            if set(aids.keys()).intersection(in_study):
                log(f"The following VCFs are already in this study at checkpoint {dest_checkpoint} and will be skipped:")
                for aid in set(aids.keys()).intersection(in_study):
                    log(f"  {aids[aid]}: {aid}")
                    aids.pop(aid)
                log("")

            if not aids:
                raise SystemExit("Nothing left to import.")

        count = self.scalar(
            self.query(
                f"""
                    SELECT count(*)
                    FROM {self.table.study.data}
                    WHERE study_name = %(study_name)s
                    AND checkpoint = %(checkpoint)d
                    AND aid IN ({self.quoted_aid_list(aids.keys())})
                    ;
                """,
                params={"study_name": src_study, "checkpoint": src_checkpoint}
            )
        )
        if not count:
            raise SystemExit(f"No variants found.")

        log(f"Adding {count:,} variants from {len(aids)} VCF{plural(len(aids))} to study {dest_study}")
        new_checkpoint = dest_checkpoint + 1
        update_table_mtime(self, self.table.study.data, partition=f"study_name={dest_study}")
        # carry previous rev forward
        # TODO: make this a parallel query
        self.query(
            f"""
                INSERT INTO "{self.table.study.data}"
                    (spans, reflen, chrom, pos, varend, varid, ref, alt, qual, filt, info, sample, study_name, checkpoint, sample_name, aid)
                SELECT
                    spans,
                    reflen,
                    chrom,
                    pos,
                    varend,
                    varid,
                    ref,
                    alt,
                    qual,
                    filt,
                    info,
                    sample,
                    %(dest_study)s,
                    %(new_checkpoint)d,
                    sample_name,
                    aid
                FROM {self.table.study.data}
                WHERE study_name = %(dest_study)s
                AND checkpoint = %(old_checkpoint)d
                ;
            """,
            params={"dest_study": dest_study, "old_checkpoint": dest_checkpoint, "new_checkpoint": new_checkpoint}
        )
        # insert new variants
        self.parallel_query(
            f"""
                INSERT INTO "{self.table.study.data}"
                    (spans, reflen, chrom, pos, varend, varid, ref, alt, qual, filt, info, sample, study_name, checkpoint, sample_name, aid)
                SELECT
                    spans,
                    reflen,
                    chrom,
                    pos,
                    varend,
                    varid,
                    ref,
                    alt,
                    qual,
                    filt,
                    info,
                    sample,
                    %(dest_study)s,
                    %(new_checkpoint)d,
                    sample_name,
                    aid
                FROM {self.table.study.data}
                WHERE study_name = %(src_study)s
                AND checkpoint = %(src_checkpoint)d
                AND aid = %(aid)s
                ;
            """,
            iterate_param={"aid": aids.keys()},
            const_params={"src_study": src_study, "src_checkpoint": src_checkpoint, "dest_study": dest_study, "new_checkpoint": new_checkpoint}
        )

        self.add_metadata_to_study(dest_study, f"checkpoint {new_checkpoint}", f"added {'; '.join([f'{v}: {k}' for k, v in aids.items()])} from study {src_study} checkpoint {src_checkpoint}")

    def sample_missingness(self, study_name, the_filter, checkpoint):
        '''
        Return a filter clause for sample missingness
        '''
        reply = self.query(
            f"""
                WITH uv AS (
                        SELECT count(*) AS total FROM {self.table.study.merged}
                        WHERE study_name = %(study_name)s
                        AND checkpoint = %(checkpoint)d
                )
                SELECT sample_name from uv, (
                SELECT sample_name, COUNT(*) as ct
                        FROM {self.table.study.data}
                        WHERE study_name = %(study_name)s
                        AND checkpoint = %(checkpoint)d
                        GROUP BY sample_name
                        ORDER BY sample_name ASC
                  )
                WHERE 1 - (ct / cast(total AS double)) {the_filter[len('S_MISS'):]}
                ;
            """,
            params={"study_name": study_name, "checkpoint": checkpoint}
        )
        if not reply:
            raise SystemExit("That filter would eliminate all samples, aborting.")

        return f"sample_name IN ({self.quoted_sample_list([s[0] for s in reply])})"

    def filter_study(self, study_name, the_filter, exclude=False): # pylint: disable=too-many-statements
        '''
        Create a new study checkpoint after applying a filter.
        If exclude is True, exclude variants that match the_filter.
        If exclude is False, include variants that match the_filter.
        '''
        self.assert_study_is_unfrozen(study_name)

        current_checkpoint = self.get_current_study_checkpoint(study_name)
        new_checkpoint = current_checkpoint + 1
        missing = False

        try:
            # Missingness filters are more restricted. They can only be applied
            # one at a time, with a simple comparison like F_MISS > 0.2
            if '_MISS' in the_filter.upper():
                missing = True
                # validate the filter first
                filter_clause = parser(the_filter, parser_type='missingness')

                # missingness requires merge
                self.merge_study(study_name, checkpoint=current_checkpoint)

                if 'S_MISS' in filter_clause:
                    prefix = ""
                    postfix = ""
                    filter_term = f'{"NOT" if exclude else ""} ( {self.sample_missingness(study_name, filter_clause, current_checkpoint)} )'

                else:
                    # NOTE: the inversed exclude is intentional here, since we're in an EXCEPT clause.
                    prefix = f"""
                        WITH uv AS (
                            SELECT chrom, pos, ref, alt FROM {self.table.study.data}
                            WHERE study_name = '{study_name}'
                            AND checkpoint = {current_checkpoint}
                            EXCEPT
                            SELECT chrom, pos, ref, alt FROM {self.table.study.merged}
                            WHERE study_name = '{study_name}'
                            AND checkpoint = {current_checkpoint}
                            AND {"" if exclude else "NOT"} ( {filter_clause} )
                        )
                    """
                    postfix = f"""
                        RIGHT JOIN uv ON
                            uv.chrom = sv.chrom
                            AND uv.pos = sv.pos
                            AND uv.ref = sv.ref
                            AND uv.alt = sv.alt
                    """
                    # All other conditions should be True for the INSERT INTO
                    filter_term = '1=1'
            else:
                prefix = ""
                postfix = ""
                filter_term = f'{"NOT" if exclude else ""} ( {parser(the_filter)} )'

            log("Applying filter")
            update_table_mtime(self, self.table.study.data, partition=f"study_name={study_name}")

            # Create partitions
            potential_aids = self.query(
                f"""
                    SELECT aid, sample_name
                    FROM {self.table.study.data} sv
                    WHERE study_name = %(study_name)s
                    AND checkpoint = {current_checkpoint}
                    GROUP BY aid, sample_name
                    ;
                """,
                params={"study_name": study_name}
            )
            queries = []
            parts = []
            for (aid, sample_name) in potential_aids:
                parts.append(
                    f"""(study_name='{study_name}', checkpoint={new_checkpoint}, sample_name='{sample_name}', aid='{aid}')""")
                # Keep the query to a reasonable size
                if len(parts) > 100:
                    queries.append(
                        f"""ALTER TABLE {self.table.study.data} ADD IF NOT EXISTS PARTITION {' PARTITION '.join(parts)};""")
                    parts = []

            if parts:
                queries.append(
                    f"""ALTER TABLE {self.table.study.data} ADD IF NOT EXISTS PARTITION {' PARTITION '.join(parts)};""")
                parts = []

            # Add partitions in parallel
            self.parallel_queries(queries)

            # Athena has a 100 writer limit, but leave some headroom and chunk INSERTs by aid
            queries = []
            for chunk in chunked([a[0] for a in potential_aids], 95):
                queries.append(f"""
                    INSERT INTO "{self.table.study.data}"
                        (spans, reflen, chrom, pos, varend, varid, ref, alt, qual, filt, info, sample, study_name, checkpoint, sample_name, aid)
                    {prefix}
                    SELECT
                        spans,
                        reflen,
                        sv.chrom,
                        sv.pos,
                        varend,
                        varid,
                        sv.ref,
                        sv.alt,
                        qual,
                        filt,
                        info,
                        sample,
                        '{study_name}',
                        {new_checkpoint},
                        sample_name,
                        aid
                    FROM {self.table.study.data} sv
                    {postfix}
                    WHERE study_name = '{study_name}'
                    AND checkpoint = {current_checkpoint}
                    AND aid IN ({self.quoted_aid_list(chunk)})
                    AND {filter_term}
                    ;
                """)

            self.parallel_queries(queries)
            self.add_metadata_to_study(study_name, f"checkpoint {new_checkpoint}", f"{'exclude' if exclude else 'include'} {'missingness ' if missing else ''}{the_filter}")

        except ParseException as err:
            error(f"Could not parse{'missingness' if missing else ''} filter:\n")
            error(err.syntax)
            raise SystemExit(err.msg)

        update_table_mtime(self, self.table.study.data, partition=f"study_name={study_name}")

        reply = self.query(
            f"""
                SELECT checkpoint, count(*), array_sort(array_agg(DISTINCT(sample_name)))
                FROM {self.table.study.data}
                WHERE study_name = %(study_name)s
                and checkpoint >= {current_checkpoint}
                GROUP BY checkpoint
                ORDER BY checkpoint
                ;
            """,
            params={"study_name": study_name}
        )

        if len(reply) == 1:
            log(f"This filter removed all variants from the study. Rolling back to previous checkpoint.")
            self.delete_study(study_name, new_checkpoint)
        elif reply[0][1] == reply[1][1]:
            log(f"Study {study_name} variants: no change ({reply[0][1]})")
        else:
            log(f"Study {study_name}:")
            log(f"  variants: {reply[0][1]} -> {reply[1][1]}")
            if reply[0][2] != reply[1][2]:
                log(f"   samples: {reply[0][2]} -> {reply[1][2]}")

    def add_vcf_partitions(self, sample_name, build, aid):
        ''' Add new partitions to the vcf tables '''
        for table in (self.table.vcf.meta, self.table.vcf.data):
            update_table_mtime(self, table)
            self.query(f"""ALTER TABLE {table} ADD IF NOT EXISTS PARTITION (sample_name='{sample_name}', build='{build}', aid='{aid}');""")

    def add_anno_partitions(self, build, anno_name, version, aid):
        ''' Add new partitions to the anno tables '''
        for table in (self.table.anno.meta, self.table.anno.data):
            update_table_mtime(self, table)
            self.query(f"""ALTER TABLE {table} ADD IF NOT EXISTS PARTITION (build='{build}', anno_name='{anno_name}', version='{version}', aid='{aid}');""")

    def get_anno_variant_count(self, aid):
        ''' Return the count of variants for the annotation with the given aid '''
        return self.scalar(self.query(
            f"""SELECT count(*) FROM {self.table.anno.data} WHERE aid = %(aid)s""",
            params={"aid": aid}
        ))

    def get_vcf_variant_count(self, aid):
        ''' Return the count of variants with the given aid '''
        return self.scalar(self.query(
            f"""SELECT count(*) FROM {self.table.vcf.data} WHERE aid = %(aid)s""",
            params={"aid": aid}
        ))
