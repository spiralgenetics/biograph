"""Internal SDK calls"""

# TODO(nils): Make sure this stuff doesn't get packaged when building
# a release SDK.

# pylint: disable=import-error
from biograph._capi import (
    Anchor,
    Variant,
    VarEdge,
    VarNode,
    VarCoverage,
    VarGraph,
    assemble,
    find_anchors,
)

from biograph.internal.GenomeGraph import (
    GNode,
    GenomeGraph
)

# TODO(nils): Fix aligner to not require swalign, which isn't compatible with python3.
#from biograph.internal.Aligner import (
#    aligner,
#    aln_to_vcf
#)
#
#from biograph.internal.KmerAsm import (
#    KmerAsm,
#    make_anchor,
#    kmer_assemble
#)

from biograph import _remove_unsightly_capi
_remove_unsightly_capi(globals(), "biograph.internal.internal")
