"""
Genotype events in a VCF

VCF must have one entry per chrom:start-end. Use `bcftools norm -m+any`
By default, check the BioGraph file for the sample to annotate.
If BG sample name is not in the vcf's sample fields, `--sample` must be specified.

Outputs all VCF records and the single SAMPLE column  (i.e. other samples are removed)
Note - for grch37d5, it seems like 16 is a pretty sweet spot for cores
    - Will run in about 40 minutes and will use all the cores most of the time
    - Give it 48 cores and it'll run for like 30 minutes, but will not be using all cores the whole time.
"""

import itertools

import biograph

import biograph.variants as bgexvar


class GTAnno:

    """
    Given a vcf_entry and a set of var.Assembly from the pair_coverage generator
    figure out the GT and GQ.
    I'll update PCMP writer to do the work of actually setting the fields
    """

    @staticmethod
    def build_probs(allele_covs):
        """
        Build the probabilities for all the gentype possibilitis

        returns a dictionary with key being the VCFs PL ordering index
        and value being a list of [Genotype, GenotypeLiklihood]
        Where GenotypeLiklihood is the phred-scaled probability that this is the correct genotype
        """
        probs = {}
        for a1_id, a2_id in itertools.combinations_with_replacement(range(len(allele_covs)), 2):
            if a1_id > a2_id:
                a1_id, a2_id = a2_id, a1_id
            gt = "{0}/{1}".format(a1_id, a2_id)
            my_gq = 0
            for myall in allele_covs:
                if myall[0] == a1_id or myall[0] == a2_id:
                    if a1_id == a2_id:  # Hom, prob of 2 copy
                        my_gq += myall[2][2]
                    else:  # Het, prob of 1 copy
                        my_gq += myall[2][1]
                else:
                    my_gq += myall[2][0]
            pl = 10 * (-my_gq/10)
            f_idx = (a2_id * (a2_id + 1) / 2) + a1_id

            probs[f_idx] = [gt, my_gq, min(99, int(round(pl)))]
        return probs

    @staticmethod
    def get_always_alt_coverages(entry):
        """
        # redundant reads cannot support reference alleles
        //... maybe ...//
        10      380     9.41%   9.93% << so a few more recalled
        01      213     5.27%   5.82% << at tiny more cost
        """
        alts = entry.vcf_entry_info.vcf_entry[4].split(',')
        ref_reads = bgexvar.ReadIdSet()
        redundant_reads = bgexvar.ReadIdSet()
        seen_reads = bgexvar.ReadIdSet()
        var_reads_list = list()
        for alt in alts:
            ec = entry.vcf_entry_info.alleles[alt].edge_coverage
            ref_reads = ref_reads + ec.reference_start
            ref_reads = ref_reads + ec.reference_end
            var_reads = ec.variant_start + ec.variant_end
            var_reads_list.append(len(var_reads))

            redundant_reads = redundant_reads + (seen_reads & var_reads)
            seen_reads = seen_reads + var_reads

        redundant_reads = redundant_reads + (seen_reads & ref_reads)

        # consolidate
        allele_covs = []
        tot_coverage = 0
        # reference coverage
        ref_unique = len(ref_reads - redundant_reads)
        allele_covs.append([0, ref_unique])
        tot_coverage += ref_unique
        for altid, alt in enumerate(alts):
            var_reads = var_reads_list[altid]
            allele_covs.append([altid + 1, var_reads])
            tot_coverage += var_reads
        return tot_coverage, allele_covs, len(redundant_reads)

    def calc_gt(self, entry):
        """
        Calculates the genotype information of this set of alleles in a vcf_entry
        returns the list of AD and a list of lists of [Genotype, RawGenotypeProb, GenotypeLiklihood]
        """
        # tot_coverage, allele_covs = self.get_uniq_allele_coverages() # Intermediate - 89.1% accuracy
        tot_coverage, allele_covs, redundant_cov = self.get_always_alt_coverages(entry)  # Intermediate2 - 89.7% accuracy
        # Over calls INS of 50-100bp as Hom alt

        if tot_coverage == 0:
            return allele_covs, None, '.'

        # add list of probabilities of 0, 1, 2 of the alleles at this position
        for allele in allele_covs:
            allele.append(biograph.genotyper(tot_coverage, allele[1]))

        probs = self.build_probs(allele_covs)

        return allele_covs, probs, redundant_cov

    def build_vcf(self, entry):
        """
        Edits the allele.vcf_entry_info.new_fmt in-place
        """
        sample_dict = GTAnno.get_blank_format()
        gt_cov, gt_data, redun_cov = self.calc_gt(entry)
        # always populated
        sample_dict["AD"] = ",".join([str(x[1]) for x in gt_cov])
        sample_dict["DP"] = str(sum([x[1] for x in gt_cov]))
        sample_dict["RC"] = str(redun_cov)

        # sometimes we have no coverage, just get rid of it
        if gt_data is None:
            entry.vcf_entry_info.new_fmt.update(sample_dict)
            return
        # This is annotation work. shouldn't be in the writer
        mgt = None
        max_gq = None
        max_pl = None
        for dat in gt_data.values():
            if max_gq is None or dat[1] > max_gq:
                mgt = dat[0]
                max_gq = dat[1]
                max_pl = dat[2]

        sample_dict["GT"] = mgt
        sample_dict["PL"] = ",".join([str(gt_data[x][2]) for x in gt_data])
        sample_dict["GQ"] = max_pl

        entry.vcf_entry_info.new_fmt.update(sample_dict)

    def parse(self, rmap, entries): #pylint: disable=unused-argument
        """
        Annotate a set of bgexvar.Assembly entries by editing the VcfEntryInfo
        yields the new bgexvar.Assembly objects
        Assumes entries have alread been put through bgexvar.generate_read_cov or generate_pair_cov
        """
        # entries = bgexvar.generate_pair_cov(rmap, entries, self.min_insert, self.max_insert)
        entries = bgexvar.generate_pair_edge_cov(entries)

        for asm in entries:
            if asm.matches_reference:
                continue
            asm.vcf_entry_info.alleles[asm.orig_alt] = asm
            asm.vcf_entry_info.alleles_remaining -= 1
            if asm.vcf_entry_info.alleles_remaining == 0:
                self.build_vcf(asm)
                yield asm
                # Remove circular reference to self (through .vcf_entry_info) to avoid memory leaks
                asm.vcf_entry_info.alleles = None

    @staticmethod
    def get_header():
        """
        Returns header information for GT annotations
        """
        return ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample Depth">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description=' +
                '"Allelic depths for the ref and alt alleles in the order listed">',
                '##FORMAT=<ID=PL,Number=G,Type=Integer,Description=' +
                '"Phred-scaled likelihoods of the genotypes at a locus">',
                '##FORMAT=<ID=RC,Number=1,Type=Integer,Description=' +
                '"Number of reads supporting redundant alleles at a locus">']

    @staticmethod
    def get_format_tags():
        """
        Pull up the tags populated by GTAnno
        returns list of tag strings
        """
        return GTAnno.get_blank_format().keys()

    @staticmethod
    def get_blank_format():
        """
        Returns the default values
        """

        return {"GT": "./.",
                "GQ": ".",
                "PL": ".",
                "DP": ".",
                "AD": ".",
                "RC": "."}
