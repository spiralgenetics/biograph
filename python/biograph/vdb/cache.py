'''
VDB query caching. All functions expect an athena.py:connect() db object as the first parameter.
'''
import hashlib
import re

from biograph.utils import timestamp
from biograph.tools.log import debug

def get_cache_entry(db, query, params=None):
    '''
    The cache_entry is an s3 Object with the location where this query
    would point if it is cacheable.

    If the query is not cacheable, return (None, None)
    If the current cache entry is valid, return (query_id, cache_entry)
    Otherwise return (None, cache_entry)
    '''
    if query.lstrip()[:6].upper() != 'SELECT':
        debug(f"non-SELECT is uncacheable")
        return (None, None)

    sha = hashlib.sha256()

    # The cache is only valid if the tables it queries have not been updated since
    # the last time this query was run.

    # The mtime is updated each time a partition is added or dropped.
    for table_class in [vars(tc).values() for tc in vars(db.table).values()]:
        for table_name in table_class:
            if re.search(r"\b" + table_name + r"\b", query):
                if params and table_name in vars(db.table.study).values() and 'study_name' in params:
                    mtime = str(get_table_mtime(db, table_name, auto_update=True, partition=f"study_name={params['study_name']}"))
                    debug(f"table {table_name} study_name {params['study_name']} {mtime}")
                    sha.update(mtime.encode())
                else:
                    mtime = str(get_table_mtime(db, table_name, auto_update=True))
                    debug(f"table {table_name} {mtime}")
                    sha.update(mtime.encode())

    # The query must match verbatim, sans leading/trailing whitespace
    sha.update(query.strip().encode())

    # The params (if any) must match verbatim in sorted order
    if params:
        for param in sorted(params):
            p = f"|{param}={params[param]}|"
            debug(f"param {p}")
            sha.update(p.encode())

    debug(f"sha: {sha.hexdigest()}")
    cache_entry = db.s3.Object(db.bucket, f"{db.path.results.cache}/{sha.hexdigest()}")
    try:
        # max s3 path length is 1024, so next() is safer than get()
        query_id = next(cache_entry.get()['Body']).decode()
        if db.s3_path_exists(f"{db.path.results.stage}/{query_id}.csv"):
            return query_id, cache_entry

    except db.s3.meta.client.exceptions.NoSuchKey:
        # miss (cache_entry does not exist)
        return (None, cache_entry)

    # miss (cache_entry points to a missing target) so delete it
    cache_entry.delete()
    return (None, cache_entry)

def save_to_cache(db, query_id, query, params=None):
    ''' Save the current query_id and result to the cache if possible '''
    _, cache_entry = get_cache_entry(db, query, params)

    if cache_entry:
        debug(f"updating cache with {query_id}")
        cache_entry.put(Body=query_id)

def fetch_from_cache(db, query, params=None, output='list', block=True):
    '''
    Retrieve results from the cache if available. If not, send the query to Athena.

    If output = list, return a simple list of results.
    If output = pandas, return a dataframe
    If output is anything else, treat it as a filename and write CSV results to it.
    '''
    debug(f"checking cache for {query.lstrip()}")

    query_id, cache_entry = get_cache_entry(db, query, params)
    # no cache_entry means not cacheable
    if cache_entry is None:
        return db.query_with_id(query, params, cache=False, block=block)

    if query_id:
        debug(f"cache hit! {query_id}")
        if output == 'pandas':
            # Return the entire contents as a dataframe
            return query_id, db.cursor._collect_result_set(query_id=query_id).as_pandas() # pylint: disable=protected-access
        if output == 'list':
            # convert to a simple list
            return query_id, db.cursor._collect_result_set(query_id=query_id).fetchall() # pylint: disable=protected-access
        # anything else is a filename
        return db.fetch_result_to_csv(query_id, output)

    if output == 'pandas':
        query_id, result = db.query_pandas_with_id(query, params=params, cache=False, block=block)
    elif output == 'list':
        query_id, result = db.query_with_id(query, params=params, cache=False, block=block)
    else:
        query_id, result = db.query_fetch_csv_with_id(query, output, params=params, cache=False, block=block)

    debug(f"no hit, save {query_id} to cache")

    # save the result to the cache if possible
    if block:
        save_to_cache(db, query_id, query, params)

    return query_id, result

def get_table_mtime(db, table, auto_update=False, partition=None):
    '''
    Return the last modified time for the table. This isn't managed by AWS so we do it ourselves.
    If no timestamp exists, and auto_update is true, set the mtime to the current time.
    '''
    if partition:
        mtime = db.s3.Object(db.bucket, f"{db.path.results.mtime}/{table}/{partition}")
    else:
        mtime = db.s3.Object(db.bucket, f"{db.path.results.mtime}/{table}")

    try:
        ret = next(mtime.get()['Body'])
        debug(f"mtime for {table} / {partition}: {ret}")
        return ret.decode()
    except db.s3.meta.client.exceptions.NoSuchKey:
        if not auto_update:
            return None
        debug(f'no mtime found for {table}/{partition}, setting to now')
        return update_table_mtime(db, table, partition)

def update_table_mtime(db, table, partition=None, ts=None):
    '''
    Update the last modified time for the table (or partition) to now.
    Returns the timestamp set on the table.
    '''
    if ts is None:
        ts = timestamp()

    if partition:
        partition.lstrip('/')
        db.s3.Object(db.bucket, f"{db.path.results.mtime}/{table}/{partition}").put(Body=str(ts))
    else:
        db.s3.Object(db.bucket, f"{db.path.results.mtime}/{table}").put(Body=str(ts))

    return ts

def clear_table_mtime(db, table, partition=None):
    '''
    Clear the last modified time for the table (or partition).
    '''
    if partition:
        partition.lstrip('/')
        db.s3.Object(db.bucket, f"{db.path.results.mtime}/{table}/{partition}").delete()
    else:
        db.s3.Object(db.bucket, f"{db.path.results.mtime}/{table}").delete()
