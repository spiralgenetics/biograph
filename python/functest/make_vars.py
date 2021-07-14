#!/usr/bin/env python3
"""
make_vars.py: Create a variant with make_vars() and Art Illumina
"""
import subprocess
import os
import json

# anchored configuration variables
opts = {
    'ref_src':            '/data/e_coli_k12.ASM584v1.fasta',
    'reference':        'e_coli_k12_ASM584v1',
    'chromosome':        'F',
    'variant_list':        None,
    'het':                False,

    'art_path':            '/opt/art/Linux64',
    'art_do_alignment':    True,
    'art_read_length':    '100',
    'art_coverage':        '50',
    'art_mean_size':    '500',
    'art_stddev_size':    '50',
    'art_seed':            '13',
    'art_qs':            '0',

    'leeway':            '400',
    'snp_count':        '0',
    'ins_count':        '0',
    'ins_size':            '1000',
    'del_count':        '0',
    'del_size':            '1000',

    'spiral_cfg':        '/opt/spiral/etc/unittest.json',
    'out':                '.',
}

# cmdline help for all options. A richer struct for both would be neater here,
# but messier on lookup.
cmdhelp = {
    'ref_src':                'path to a fasta which is the source of a reference (usually under /data/)',
    'reference':            'reference name to use (must exist under /reference/)',
    'chromosome':            'limit variations to a specific chromosome',
    'variant_list':            'path to the variant_list SVS file created by make_vars (required for --premade_reads)',
    'het':                    'produce heterozygous variants',

    'art_path':            'full path to art_illumina installation',
    'art_do_alignment':    'output SAM when generating reads (can be large)',
    'art_read_length':    'read length to be simulated',
    'art_coverage':        'read coverage to be simulated',
    'art_mean_size':    'the mean size of DNA fragments for simulation',
    'art_stddev_size':    'the standard deviation of DNA fragment size for simulation',
    'art_seed':            'random seed for simulation',
    'art_qs':            'art simulation quality shift (negative for more error)',

    'leeway':            'leeway parameter for variant creation',
    'snp_count':        'number of SNPs',
    'ins_count':        'number of insertions',
    'ins_size':            'size of each insertion',
    'del_count':        'number of deletions',
    'del_size':            'size of each deletion',
    'out':               'output directory',

    'spiral_cfg':        'path to spiral config.json'
}

def generate_reads(reads_file, sequence, coverage):
    """
        Run art_illumina to generate simulated reads
    """
    art_path = os.path.expanduser(opts['art_path'])

    do_alignment = '-sam' if opts['art_do_alignment'] else '-na'

    cmd = [
        '%s/art_illumina' % art_path,
        '-1', '%s/%s' % (art_path, 'Illumina_profiles/EmpMiSeq250R1.txt'),
        '-2', '%s/%s' % (art_path, 'Illumina_profiles/EmpMiSeq250R2.txt'),
        '-o', '%s/%s' % (opts['out'], reads_file),
        '-i', sequence,
        '-l', opts['art_read_length'],
        '-f', coverage,
        '-m', opts['art_mean_size'],
        '-s', opts['art_stddev_size'],
        '-rs', opts['art_seed'],
        '-qs', opts['art_qs'],
        '-qs2', opts['art_qs'],
        '-p',
        do_alignment
    ]

    subprocess.call(cmd)

if __name__ == '__main__':

    import argparse

    # Build out a dynamic argument list based on opts
    PARSER = argparse.ArgumentParser(description='Anchored Assembly functional tests')
    for opt in sorted(opts):
        if opts[opt] in [True, False]:
            action = 'store_true'
        else:
            action = 'store'

        PARSER.add_argument('--%s' % opt, default=opts[opt], action=action, help='%s (default: %s)' % (cmdhelp[opt], opts[opt]))

    (ARGS, ARGV) = PARSER.parse_known_args()

    for arg in vars(ARGS):
        newarg = getattr(ARGS, arg)
        if opts[arg] != newarg:
            print 'Using %s: %s' % (arg, newarg)
        opts[arg] = newarg


    # Import libspiral C++ library.
    # pylint: disable=import-error
    import libspiral

    variant_sequence = '%s/variant.fasta' % opts['out']
    variant_list = '%s/variant_list' % opts['out']

    with open(opts['spiral_cfg']) as cfg_file:
        cfg = json.load(cfg_file)

    # libspiral.load_config(opts['spiral_cfg'])

    # Build a local ref if needed
    if not os.path.exists(os.path.join(cfg['reference_path'], opts['reference'])):
        print 'Local reference needed for make_vars. Importing ref_src %s' % opts['ref_src']
        libspiral.build_ref(opts['reference'], opts['ref_src'])

    for i in [
        'leeway',
        'snp_count',
        'ins_count', 'ins_size',
        'del_count', 'del_size',
    ]:
        opts[i] = int(opts[i])

    # Assinging an empty string to a dict value always results in None.
    # But we need to pass a literal empty string (not None) to make_vars
    # if no chromosome is selected, so use a real variable.
    chromosome = opts['chromosome']
    if chromosome == None:
        chromosome = ''

    with open(variant_list, 'w') as list_file:

        make_vars = libspiral.make_vars(opts['reference'], opts['leeway'], int(opts['art_read_length']), list_file, chromosome)

        # snp is just a name, but doesn't seem to show up in vars_list...?
        for i in range(opts['snp_count']):
            make_vars.snp('snp_%d' % i)

        # insertions and deletions are name + size

        for i in range(opts['ins_count']):
            make_vars.random_insert('ins_%d:%d' % (i, opts['ins_size']), opts['ins_size'])

        for i in range(opts['del_count']):
            make_vars.random_delete('del_%d:%d' % (i, opts['del_size']), opts['del_size'])

        # save the variant sequence
        with open(variant_sequence, 'w') as seq_file:
            make_vars.print_sequence(seq_file)


    variant_sequence = '%s/variant.fasta' % opts['out']

    # We already have the variants. Now generate a variant sequence with no changes.
    # The reference file can't be used directly, since we may have specified a single chromosome.
    if opts['het']:
        print 'Generating non-variant heterozygous complements'

        no_variant_sequence = '%s/no_variants.fasta' % opts['out']

        # No changes here, so no need to save the variant list
        with open(os.devnull, 'w') as list_file:
            no_make_vars = libspiral.make_vars(opts['reference'], opts['leeway'], int(opts['art_read_length']), list_file, chromosome)

            with open(no_variant_sequence, 'w') as seq_file:
                no_make_vars.print_sequence(seq_file)

        generate_reads('%s/reads' % opts['out'], variant_sequence, opts['art_coverage'] / 2)

        print 'Generating additional reads (with no variants) for heterozygous variant complements.'

        # This must match the definition in test_00_make_vars
        no_variant_sequence = '%s/no_variants.fasta' % opts['out']
        normal_reads = '%s/no_variant_reads' % opts['out']
        generate_reads(normal_reads, no_variant_sequence, opts['art_coverage'] / 2)

    else:
        generate_reads('%s/reads' % opts['out'], variant_sequence, opts['art_coverage'])

