#!/usr/bin/env python3
'''
Query the VDB
'''
import argparse
import logging
import sys
import pprint

import pandas as pd
import orjson as json
from biograph.internal import vdb

def parse_args(clargs):
    ''' biograph vdb query args '''
    parser = argparse.ArgumentParser(
        description='Query the Spiral Variant DataBase (VDB)'
    )

    if not clargs:
        clargs.append('--help')

    parser.add_argument('sample', type=str, nargs='?', help='Match this sample or annotation name')
    parser.add_argument('-a', '--aid', type=str, help='Match this analysis id')
    parser.add_argument('-g', '--group', help='Use this VDB group', default=None)
    parser.add_argument('--all', action='store_true', help='List every available analysis')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show full details for each analysis')
    parser.add_argument('--annotation', action='store_true', help='Show a list of all available annotations')
    parser.add_argument('--lookup', type=str, help='Look up annotations with this VARID')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    return parser.parse_args(clargs)

def list_entries(entries, verbose, annotation):
    ''' Pretty print query results '''

    if verbose:
        for entry in entries:
            if annotation:
                aid, sample, refname, build, imported_on, description, version, refhash, header = entry
            else:
                aid, sample, refname, build, imported_on, description, refhash, header = entry

            print(f"""
name: {sample}
aid: {aid}
build: {build}
refname: {refname}
refhash: {refhash}
imported_on: {imported_on}"""
                 )

            if not pd.isna(description):
                print(f'description: {description}')
            if annotation:
                print(f'version: {version}')

            print()

            for line in header.split('\n'):
                if line.lower().startswith(
                        (
                            '##filter',
                            '##filedate',
                            '##reference',
                            '##info',
                            '##format',
                            '##alt',
                            '##contig',
                            '##refhash',
                            '##fileformat',
                            '#chrom',
                            '##sequence-region',
                            '##gff-version',
                        )
                    ):
                    continue
                if line.lower().startswith('##source='):
                    subs = line.split(',')
                    print(subs[0][2:])
                    for sub in subs[1:]:
                        print(' ', sub)
                else:
                    print(line[2:])

            print('-=' * 10)
        return

    # not verbose
    if annotation:
        print(f"{'name':<16}  {'version':<12}  {'build':<7}  {'analysis_id':<36}  {'imported_on':<24}  description")
    else:
        print(f"{'name':<16}  {'refname':<10}  {'build':<7}  {'analysis_id':<36}  {'imported_on':<24}  description")

    for entry in entries:
        if annotation:
            aid, sample, _, build, imported_on, description, version = entry
            print(f"{sample:<16}  {version:<12}  {build:<7}  {aid:<36}  {imported_on.ctime():<24}  {'' if pd.isna(description) else description}")
        else:
            aid, sample, refname, build, imported_on, description = entry
            print(f"{sample:<16}  {refname:10}  {build:<7}  {aid:<36}  {imported_on.ctime():<24}  {'' if pd.isna(description) else description}")

def lookup(db, varid_lookup):
    ''' look up an annotation '''
    urls = {
        'ClinVar': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/{varid}/',
        'ensGene': 'https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={varid}',
        'Ensembl': 'https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={varid}',
        'knownGene': 'https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={varid}',
        'GCF_000001405.39': 'https://www.ncbi.nlm.nih.gov/search/all/?term={varid}',
        'ncbiRefSeq': 'https://www.ncbi.nlm.nih.gov/search/all/?term={varid}',
        'GeneCards': 'https://www.genecards.org/cgi-bin/carddisp.pl?{varid}',
        'UniProt': 'https://www.uniprot.org/uniprot/{varid}',
        'OMIM': 'https://omim.org/entry/{varid}',
    }
    pp = pprint.PrettyPrinter(indent=4)
    query = f'''
        SELECT name, chrom, pos, varid, build, t_info
        FROM {db.path.annotation.data_table}
        WHERE varid = '{varid_lookup}'
    '''
    for (name, chrom, pos, varid, build, t_info) in db.query(query):
        if name in urls:
            url = urls[name].format(varid=varid)
        else:
            url = ''
        print(f"{name}:{varid} at {chrom}:{pos} on {build} {url}")
        pp.pprint(json.loads(t_info))
        print('')

def main(clargs):
    ''' the main event '''
    args = parse_args(clargs)

    logLevel = logging.DEBUG if args.debug else logging.WARNING
    logging.basicConfig(stream=sys.stderr, level=logLevel, format='%(message)s')

    db = vdb.connect(group=args.group)

    if args.lookup:
        lookup(db, args.lookup)
        exit(0)

    if args.annotation:
        table = db.path.annotation.meta_table
        sample_field = 'name'
        ready = db.get_crawler_state(db.path.annotation.crawler) == 'READY'
    else:
        table = db.path.vcf.meta_table
        sample_field = 'sample'
        ready = db.get_crawler_state(db.path.vcf.crawler) == 'READY'

    if not ready:
        logging.warning(f"NOTE: The crawler is currently running. Some data may not yet be indexed.")

    if not db.query(f"SHOW TABLES LIKE '{table}'"):
        raise SystemExit(f"{db.group} is empty. Run 'biograph vdb import --group {db.group}' to import VCF data,\nor 'biograph vdb group --crawl {db.group}' to update the index.")

    if not args.annotation:
        logging.warning(f"vdb group '{db.group}' {'(frozen)' if db.is_frozen() else ''}")

    query_filters = []
    if not args.all:
        if args.aid:
            query_filters.append(f'''
                AND aid = '{args.aid}'
            ''')

        if args.sample:
            if args.annotation:
                query_filters.append(f"AND name = '{args.sample}'")
            else:
                query_filters.append(f"AND sample = '{args.sample}'")

    query_filter = '\n'.join(query_filters)
    fields = [
        'aid',
        sample_field,
        'refname',
        'build',
        'imported_on',
        'description'
    ]

    if args.annotation:
        fields.append('version')

    if args.verbose:
        fields.extend(['refhash', 'header'])

    query = f'''
    SELECT
        {','.join(fields)}
    FROM
        {table}
    WHERE
        1=1
    {query_filter}
    ORDER BY {sample_field} ASC, imported_on DESC
    ;
    '''

    logging.debug(query)

    list_entries(db.query(query, cache=False), args.verbose, args.annotation)

if __name__ == '__main__':
    main(sys.argv[1:])
