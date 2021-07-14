"""
    __init__.py

    Path fix-ups
"""
import os
import sys

def rchop(in_string, ending):
    """ Strip a matching end pattern from a string. """
    if in_string.endswith(ending):
        return in_string[:-len(ending)]
    return in_string

# Include the full path to python/ and python/functest/
FUNCTEST_PATH = rchop(os.path.dirname(os.path.abspath(__file__)), '/utils')
