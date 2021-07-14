"""
Helper object for tracking specific VCF entry information while processing associated assemblies
"""
import logging

import tabix

import biograph.variants as bgexvar


class VcfEntryInfo:

    """
    Holds tracking information for a specific VCF entry while processing all the associated assemblies
    call static method parse_vcf over a fetched region
    """
    __slots__ = ["vcf_entry", "total_alleles", "alleles_remaining", "alleles", "new_fmt", "sample_field"]
    def __init__(self, vcf_entry, sample_field, total_alleles):
        self.vcf_entry = vcf_entry
        self.sample_field = sample_field
        self.total_alleles = total_alleles
        self.alleles_remaining = total_alleles
        self.alleles = dict()
        self.new_fmt = {}  # new format dict that's populated

    def collect_allele(self, allele):
        """
        Marks that this allele has been processed
        """
        self.alleles[allele.orig_alt] = allele
        self.alleles_remaining = self.alleles_remaining - 1

    @staticmethod
    def parse_region(vcf_file, chrom, start, end, sample_index, passonly=None, phasing=False):
        """
        iterates each of the variants in chrom/start/end region,
        returns VCFEntryInfo generator
        if passonly is None: return all entries
        if passonly is 1: return all PASS entries
        if passonly is -1 return all non-PASS entries
        """
        try:
            fetched = tabix.open(vcf_file).query(chrom, start, end)
        except tabix.TabixError as e:
            logging.warning(f"No variants found in region {chrom}:{start}-{end} ({e})")
            return []

        return VcfEntryInfo.parse_entries(fetched, start, end, sample_index, passonly, phasing)

    @staticmethod
    def parse_entries(items, start, end, sample_index, passonly=None, phasing=False):
        """
        Parse list of vcf entries. Must have vcf columns 0, 4
        if passonly is None: return all entries
        if passonly is True: return all PASS entries
        if passonly is False return all non-PASS entries
        """
        for asm_id, vcf_entry in enumerate(items):
            # tabix returns any variant that overlaps a region.
            # Only report variants that start within the region.
            if int(vcf_entry[1]) <= start or int(vcf_entry[1]) > end:
                continue
            if passonly is not None and passonly != ("PASS" in vcf_entry[6]):
                continue

            if "N" in vcf_entry[3] or "N" in vcf_entry[4]:
                logging.debug("Entry Skipped for 'N' in REF/ALT @ %s:%s", vcf_entry[0], vcf_entry[1])
                continue

            alt_alleles = vcf_entry[4].split(',')
            left = int(vcf_entry[1]) - 1
            right = left + len(vcf_entry[3])
            if sample_index is not None:
                sample_field = vcf_entry[9 + sample_index]
            else:
                sample_field = None
            vcf_entry_info = VcfEntryInfo(vcf_entry[:9], sample_field, len(alt_alleles))

            phase_ids = None
            if phasing:
                formats = vcf_entry[8].split(':')
                phasing_format_offset = formats.index("PI")
                phase_ids = bgexvar.PhaseSet.from_format_fields(phasing_format_offset, vcf_entry[9:])

                if len(alt_alleles) > 1:
                    logging.warning("Phase ids not supported for multiallelic VCF entries")

            for alt in alt_alleles:
                allele = bgexvar.Assembly(left, right, alt.replace("N", ""), asm_id)
                allele.vcf_entry_info = vcf_entry_info
                allele.orig_alt = alt

                if phase_ids is not None:
                    allele.phase_ids = phase_ids
                yield allele
