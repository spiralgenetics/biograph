"""Experimental SDK calls relating to variants"""

# pylint: disable=import-error
from biograph._capi import (
    EdgeCoverage,
    Assembly,
    apply_edges,
    generate_pair_edge_cov,
    generate_read_cov,
    generate_pair_cov,
    add_ref_assemblies,
    trim_ref,
    dedup_cov_reads,
    ReadCoverage,
    ReadCoverageRead,
    PairEdgeCovGenerator,
    ParallelDiscover,
    RefMap,
    ReadIdSet,
    BigReadIdSet,
    join_phases,
    JoinPhasesGenerator,
    propagate_subassembly_coverage,
    split_phases,
    SplitPhasesGenerator,
    resolve_phase_conflicts,
    ResolvePhaseConflictsGenerator,
    PhaseSet,
    limit_alleles,
    LimitAllelesGenerator,
    filter_dup_align,
    FilterDupAlignGenerator,
    align_reads,
    AlignReadsGenerator,
    AlignedRead,
    discover_branch,
    discover_push_to_pair,
    update_rc_seqset_entries,
    make_ref_assemblies,
    graph_trim_ref,
    place_pair_cov,
    PlacePairCovGenerator,
)

# Remove unsightly _capi from published class names.
from biograph import _remove_unsightly_capi
_remove_unsightly_capi(globals(), "biograph.variants")
