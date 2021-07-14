#!/usr/bin/env python3

"""Wrapper around pylint to use a different config file depending on pylint version"""

from __future__ import print_function

pylintrc = "tools/pylintrc"

import tempfile

temp_pylintrc = tempfile.mktemp()

import os
from shutil import copyfile

copyfile(pylintrc, temp_pylintrc)
os.environ["PYLINTRC"] = temp_pylintrc
# Pylint likes to store temporary files in $HOME
#os.environ["HOME"] = "/"

import pylint

if pylint.__pkginfo__.version != "2.3.1":
    print("Non-hermetic pylint being used?  Should use the one from tool_requirements.txt.")
    import sys
    sys.exit(1)

pylint.run_pylint()
