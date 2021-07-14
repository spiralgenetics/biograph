"""

The ``biograph`` library provides methods for querying files in the Spiral
BioGraph format. BioGraphs and references can be loaded and rapidly searched
for aribtrary k-mers. Reads can also be assembled against a reference and
checked for variants.

For a quick start guide and full reference for the SDK, see:

    https://www.spiralgenetics.com/user-documentation

Usage:
    >>> from biograph import BioGraph, Reference
    >>> my_biograph = BioGraph('datasets/lambdaToyData/benchmark/family_lambda.bg')
    >>> my_ref = Reference('/reference/hs37d5/')
    >>> results = my_biograph.seqset.find('ACTG')
"""

import importlib
import sys

# Load the appropriate extension module for this version of python
# pylint:disable=import-error
_capi = importlib.import_module(f"biograph._capi_{sys.version_info.major}{sys.version_info.minor}")
sys.modules["biograph._capi"] = _capi

# pylint:disable=wrong-import-position
from biograph._capi import (
    BioGraph,
    CacheStrategy,
    Metadata,
    Readmap,
    ReadmapRead,
    ReadmapPairStats,
    Reference,
    ReferenceContext,
    ReferenceRange,
    Seqset,
    SeqsetEntry,
    Sequence,

    build_revision,
    version,

    #set_spiral_logging_target,
    #log_build_stamp,
)

from biograph.utils import (
    Assembly,
    find_breakpoint_variants,
    find_region_variants,
    genotyper,
    visualize,
)

__all__ = [
    'BioGraph',
    'CacheStrategy',
    'Metadata',
    'Readmap',
    'ReadmapRead',
    'ReadmapPairStats',
    'Reference',
    'ReferenceContext',
    'ReferenceRange',
    'Seqset',
    'SeqsetEntry',
    'Sequence',
    'build_revision',
    'version',
]

def _remove_unsightly_capi(target_dict, new_module):
    """Removes the unsightly ._capi component from classes exposed through the C API"""
    for v in target_dict.values():
        if hasattr(v, "__module__") and v.__module__ == "biograph._capi":
            try:
                v.__module__ = new_module
            except AttributeError:
                # Some objects (like boost function objects) don't allow modifying
                # __module__.
                pass

_remove_unsightly_capi(globals(), "biograph")
