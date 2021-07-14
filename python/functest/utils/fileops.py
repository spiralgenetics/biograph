"""
    fileops.py

    File operation utility functions go here.
"""
import gettext
import subprocess
import hashlib
import gzip
import re
import os

from python.functest.utils.defaults import (
    FTEST_OUTPUT
)

gettext.NullTranslations().install()

def opener(file_name, options='r', test=None):
    """
        Open gzip or other file. Returns a file handle to the opened file.
        Files are opened read-only by default.

        'options' is not required, and can be any valid open or gzip.open flag.

        If file_name does not contain a '/', assume it is in test.data_dir or FTEST_OUTPUT
    """
    if not "/" in file_name:
        if test:
            file_name = test.data_dir + "/" + file_name
        else:
            file_name = FTEST_OUTPUT + "/" + file_name

    if file_name.endswith('.gz'):
        return gzip.open(file_name, options)
    return open(file_name, options)

def sha1_file(file_name, test=None):
    """
        Return the printable sha1 digest of a given file.

        If file_name doesn't contain a / then assume it's in test.data_dir or FTEST_OUTPUT
    """
    if not "/" in file_name:
        if test:
            file_name = test.data_dir + "/" + file_name
        else:
            file_name = FTEST_OUTPUT + "/" + file_name

    sha1 = hashlib.sha1()
    with open(file_name, "rb") as f:
        while True:
            data = f.read(160)
            if not data:
                break
            sha1.update(data)
    return sha1.hexdigest()

def match_all_uncommented_lines(file_name, pattern, negative=False, test=None):
    """
        re.search all uncommented lines for pattern.

        If negative is False (default), pattern must match every uncommented line.
        If negative is True, pattern must not match ANY uncommented line.

        Return True if pattern match succeeds, False otherwise.
    """
    if test:
        file_name = test.data_dir + "/" + file_name
    else:
        file_name = FTEST_OUTPUT + "/" + file_name

    with open(file_name, "r") as f:
        for line in f:
            # comment? do nothing.
            if re.search('^#', line):
                pass
            elif negative:
                # negative match and we matched? fail.
                if re.search(pattern, line):
                    # print "matched pattern %s when negative, failing." % pattern
                    return False
            # positive match and we didn't match? fail.
            elif not re.search(pattern, line):
                # print "matched pattern %s when positive, failing." % pattern
                return False

    return True

def samtools_flagstat(file_name=None, test=None):
    """
        run 'samtools flagstat' on a file and return the results.
    """
    if not file_name:
        return False

    if not "/" in file_name:
        if test:
            file_name = test.data_dir + "/" + file_name
        else:
            file_name = FTEST_OUTPUT + "/" + file_name

    return subprocess.check_output(['/usr/bin/samtools', 'flagstat', file_name])

def cat_file(file_name=None, test=None):
    """
        Return the contents of file_name.
    """
    if not file_name:
        return False

    if not "/" in file_name:
        if test:
            file_name = test.data_dir + "/" + file_name
        else:
            file_name = FTEST_OUTPUT + "/" + file_name

    with open(file_name, 'r') as f:
        return f.read()

def which(program):
    """
        Like unix which(): find the exe in your PATH, if any
    """
    def is_exe(fpath):
        """ returns true if the file exists and is executable """
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath = os.path.split(program)[0]
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

if __name__ == '__main__':
    print("File operation utility library.")
