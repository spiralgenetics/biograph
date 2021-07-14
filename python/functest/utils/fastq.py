"""
    fastq.py

    Utility functions for manipulating fastq files.
"""
from __future__ import print_function

import re
from python.functest.utils.fileops import opener


def count_fastq_records(file_name):
    """
        Count all of the fastq records in a file. Uncompress the file if needed.

        Returns the number of records found. Raises an exception if an invalid
        fastq file is detected.

        This may eventually be replaced by the UMich validator. See SG-608.
    """
    records = 0

    with opener(file_name, 'r') as f:
        pos = 1
        for num, line in enumerate(f):
            if pos == 1:
                # Line 1 must start with an @
                if not re.search(r'^@', line):
                    raise Exception('%s is not a valid fastq file, invalid identifier line %d' % (file_name, num + 1))
                pos = 2

            elif pos == 2:
                # Line 2 is the sequence
                if not re.search(r'^[ACGTN]', line):
                    raise Exception('%s is not a valid fastq file, invalid sequence line %d' % (file_name, num + 1))
                pos = 3

            elif pos == 3:
                # Line 3 starts with a +
                if not re.search(r'^\+', line):
                    raise Exception('%s is not a valid fastq file, invalid plus line %d' % (file_name, num + 1))
                pos = 4

            elif pos == 4:
                # There are a variety of encodings for the quality line. Let's just count it.
                records += 1
                pos = 1

            else:
                raise Exception('Invalid state: %d' % pos)

    return records


def fetch_fastq_record(file_handle):
    """
        Fetch four lines from the given filehandle. Return only lines 1, 2,
        and 4 (line 3 is assumed to be worthless, since Spiral discards
        everything after the +)

        Returns the munged record, or None if EOF is reached.
    """
    line_one = file_handle.readline()

    if not line_one:
        return None

    # line 2
    line_two = file_handle.readline()

    # discard line 3
    file_handle.readline()

    # line 4
    line_four = file_handle.readline()

    # quick and dirty record check.
    if not line_one.startswith('@') or not line_two or not line_four:
        raise Exception('Not a valid fastq record')

    # Return the whole thing
    return '%s%s%s' % (line_one, line_two, line_four)


def compare_paired_fastq_files(file_one, file_two, paired_file):
    """
        Ensure that paired_file consists of interleaved records from file_one and file_two.

        If any files are compressed, uncompress them on the fly.

        Returns True if the files match, False otherwise.
    """
    fh_1 = opener(file_one)
    fh_2 = opener(file_two)
    paired = opener(paired_file)

    record = 1

    while True:
        paired_record = fetch_fastq_record(paired)
        # EOF?
        if not paired_record:
            # If there are records left in either of the other files, we lose
            if fetch_fastq_record(fh_1) or fetch_fastq_record(fh_2):
                print('Record count mismatch')
                return False
            return True

        if record % 2:
            # odd records read file_one
            compare_record = fetch_fastq_record(fh_1)
        else:
            # even records read file_two
            compare_record = fetch_fastq_record(fh_2)

        # The records should match
        if paired_record != compare_record:
            print('Record #%d mismatch' % record)
            print('Expected: "%s"' % compare_record)
            print('Actual  : "%s"' % paired_record)
            return False

        # Increment record count and repeat
        record += 1


def compare_fastq_files(file_one, file_two):
    """
        Ensure that the two fastq files are identical.

        If any files are compressed, uncompress them on the fly.

        Returns True if the files match, False otherwise.
    """
    fh_1 = opener(file_one)
    fh_2 = opener(file_two)

    while True:
        file_one_record = fetch_fastq_record(fh_1)
        # EOF?
        if not file_one_record:
            # If there are records left in either of the other files, we lose
            return not fetch_fastq_record(fh_2)

        file_two_record = fetch_fastq_record(fh_2)

        # The records should match
        if file_one_record != file_two_record:
            return False
