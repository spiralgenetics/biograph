"""
Resolves conflicts between variants with the same phase id that overlap
"""

import logging
import biograph.variants as bgexvar

class Resolver:
    "Counts the phase conflicts per region"
    def __init__(self):
        self.num_conflicts = 0

    @classmethod
    def is_preferred(cls, a, b):
        "Returns true if a phase conflict between a and b should prefer a"
        # Prefer things with shorter size changes.
        abs_svlen_diff = abs(cls.get_svlen(a)) - abs(cls.get_svlen(b))
        if abs_svlen_diff:
            return abs_svlen_diff < 0
        # Prefer things with longer sequences of either reference or variant
        max_svlen_diff = max(len(a.seq), cls.get_reflen(a)) - max(len(b.seq), cls.get_reflen(b))
        if max_svlen_diff:
            return max_svlen_diff > 0
        # Prefer thigns that encompass more of reference
        reflen_diff = cls.get_reflen(a) - cls.get_reflen(b)
        if reflen_diff:
            return reflen_diff > 0
        svlen_diff = cls.get_svlen(a) - cls.get_svlen(b)
        if svlen_diff:
            # Prefer deletes (smaller svlen) to inserts
            return svlen_diff < 0
        phasecount_diff = len(a.phase_ids) - len(b.phase_ids)
        if phasecount_diff:
            # Prefer things that are part of many phases
            return phasecount_diff > 0
        return False

    @staticmethod
    def get_reflen(a):
        "Returns the number of bases in reference spanned by the give assembly"
        return a.right_offset - a.left_offset

    @classmethod
    def get_svlen(cls, a):
        "Returns the change in the number of bases with respect to reference"
        return len(a.seq) - cls.get_reflen(a)

    def on_conflict(self, a, b, ids):
        "Resolves a phase id conflict between two assemblies"
        self.num_conflicts += 1

        # Verbosely log the first conflict of each region.
        if self.num_conflicts == 1:
            a_phase_ids = str(a.phase_ids)
            b_phase_ids = str(b.phase_ids)

        if self.is_preferred(a, b):
            b.phase_ids = b.phase_ids.difference(ids)
            if self.num_conflicts == 1:
                logging.warning(f"Phase conflict detected between {a} and {b}, with a phases = {a_phase_ids} b phases = {b_phase_ids} -> {b.phase_ids} conflicts = {ids}; disabling further phase warnings in this region")
        else:
            a.phase_ids = a.phase_ids.difference(ids)
            if self.num_conflicts == 1:
                logging.warning(f"Phase conflict detected between {a} and {b}, with a phases = {a_phase_ids} -> {a.phase_ids} b phases = {b_phase_ids} conflicts = {ids}; disabling further phase warnings in this region")


class PhaseConflictResolver:
    """
    Resolves conflicts between variants with the same phase id that overlap.
    """

    @staticmethod
    def get_header():
        "Returns VCF header lines added"
        return []
    @staticmethod
    def get_format_tags():
        "Returns VCF format tags added"
        return []
    @staticmethod
    def get_blank_format():
        "Returns blank format values for format tags added"
        return {}

    @staticmethod
    def parse(_, entries):
        "Resolves phases conflicts between the given entries"
        resolver = Resolver()
        for a in bgexvar.resolve_phase_conflicts(resolver.on_conflict, entries):
            yield a
        if resolver.num_conflicts:
            logging.warning(f"{resolver.num_conflicts} phase conflicts resolved")
