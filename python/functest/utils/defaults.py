"""
    defaults.py

    Global defaults
"""
from os.path import (
    dirname,
    abspath
)
from os import makedirs, getenv

import time
import sys


def rchop(in_string, ending):
    """ Strip a matching end pattern from a string. """
    if in_string.endswith(ending):
        return in_string[:-len(ending)]
    return in_string

# Tests should be run in this order. Anything not in this list is run
# alphabetically once this list completes.
TEST_ORDER = [
    'basic',
    'spec',
    'seqset',
    'sam'
]

# Timeout for long-running commands.
# A half hour should be plenty for any single operation
TIMEOUT = '1800'

# Compute the canonical path to functest, python, golden_dir, etc.
FUNCTEST_PATH = rchop(dirname(abspath(__file__)), '/utils')
PY_PATH = rchop(FUNCTEST_PATH, '/functest')
SPIRAL_ROOT = "/src"
if getenv("TEST_UNDECLARED_OUTPUTS_DIR"):
    # running under bazel
    FTEST_OUT_BASE = getenv("TEST_UNDECLARED_OUTPUTS_DIR")
    GOLDEN_DIR = "golden/ftest"
elif getenv("TEST_TMPDIR"):
    # running under bazel
    FTEST_OUT_BASE = getenv("TEST_TMPDIR")
    GOLDEN_DIR = "golden/ftest"
else:
    FTEST_OUT_BASE = "/out"
    GOLDEN_DIR = '%s/golden/ftest' % SPIRAL_ROOT
FTEST_OUTPUT = '%s/%s' % (FTEST_OUT_BASE, time.strftime('%Y-%m-%d.%H_%M_%S'))

makedirs(FTEST_OUTPUT)
with open('%s/ftest_cmdline' % FTEST_OUTPUT, 'w') as cmdline:
    cmdline.write(' '.join(sys.argv) + '\n')
