"""
Annotate VCFs based on counting how many alignment sections there are for each read
"""

import itertools

import biograph

import biograph.variants as bgexvar


class ACAnno:
    def build_vcf(self, entry):
        """
        Edits the allele.vcf_entry_info.new_fmt in-place
        """
        sample_dict = ACAnno.get_blank_format()
        
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
        Assumes entries have alread been put through bgexvar.generate_read_cov and generate_pair_cov
        """
        # entries = bgexvar.generate_pair_cov(rmap, entries, self.min_insert, self.max_insert)
        entries = bgexvar.align_count(entries)

        for asm in entries:
            if not asm.matches_reference:
                sample_dict = self.get_blank_format()
                sample_dict["AC_LR"] = asm.align_count.local_read_lens
                sample_dict["AC_AB"] = asm.align_count.local_aligned_bases
                sample_dict["AC_TA"] = asm.align_count.tot_aligned_bases
                asm.vcf_entry_info.new_fmt.update(sample_dict)
            yield asm

    @staticmethod
    def get_header():
        """
        Returns header information for GT annotations
        """
        return ['##FORMAT=<ID=AC_LR,Number=1,Type=Integer,Description="Sum of lengths of distinct reads in the variant">',
                '##FORMAT=<ID=AC_AB,Number=1,Type=Integer,Description="Sum of total aligned bases in this variant">',
                '##FORMAT=<ID=AC_TA,Number=1,Type=Integer,Description="Sum of total aligned bases in other strands for reads in this variant">',
        ]

    @staticmethod
    def get_format_tags():
        """
        Pull up the tags populated by GTAnno
        returns list of tag strings
        """
        return ACAnno.get_blank_format().keys()

    @staticmethod
    def get_blank_format():
        """
        Returns the default values
        """

        return {"AC_LR": ".",
                "AC_AB": ".",
                "AC_TA": "."}
