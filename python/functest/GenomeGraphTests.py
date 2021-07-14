"""
Functional tests for GenomeGraph
"""
import logging
import unittest
from io import StringIO
import vcf
from biograph import Sequence
from biograph.internal import GenomeGraph

class FakeReference:

    """
    Behaves like a biograph.Reference, but doesn't actually need to be created so that we can do testing on
    reference sequences created on the fly
    """

    def __init__(self, chrom, sequence):
        self.ref_dict = {chrom: sequence}

    def make_range(self, chrom, start, end, use_exact_loci=True):
        """
        Fakes the make_range object
        """
        class FakeRefRange:
            # pylint: disable=too-few-public-methods
            """
            #'chromosome', 'end', 'scaffold', 'sequence', 'size', 'start'
            """

            def __init__(self, chrom, start, end, sequence, exact=False):
                self.chromosome = chrom
                self.start = start
                self.end = end
                self.size = end - start
                self.scaffold = chrom
                self.sequence = Sequence(sequence)
                self.exact = exact

        return FakeRefRange(chrom, start, end, self.ref_dict[chrom][start:end], use_exact_loci)

    def find_ranges(self, chrom, start, end, use_exact_loci=True):
        """
        Mockup of find_ranges for the fake ref
        """
        return [self.make_range(chrom, start, end, use_exact_loci)]

    @property
    def scaffold_lens(self):
        """
        Returns a dict of the references' lengths
        """
        ret = {}
        for key in self.ref_dict:
            ret[key] = len(self.ref_dict[key])
        return ret


class GenomeGraphTestCases(unittest.TestCase):

    """
    Tests building and using a GenomeGraph
    """

    def setUp(self):
        self.chrom = "f"
        self.ref_seq = "ATCAAGCACTA"
        self.fake_reference = FakeReference(self.chrom, self.ref_seq)
        # Header for doing vcf work - this should probably be somewhere else
        self.vcf = StringIO()
        self.vcf.write("##fileformat=VCFv4.1\n")
        self.vcf.write("##fileDate=20161123\n")
        self.vcf.write("##source=Spiral Genetics v2.1\n")
        self.vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        self.vcf.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples\">\n")
        self.vcf.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural Variant Type\">\n")
        self.vcf.write("##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakends\">\n")
        self.vcf.write(
            "##INFO=<ID=AID,Number=.,Type=Integer,Description=\"Assembly IDs used in constructing this variant\">\n")
        self.vcf.write(
            "##INFO=<ID=AMBCOUNT,Number=1,Type=Integer,Description=\"Count of alternate locations for this end of an ambiguous breakend\">\n")
        self.vcf.write(
            "##INFO=<ID=AMBMATES,Number=1,Type=Integer,Description=\"Count of possible mate locations of an ambiguous breakend\">\n")
        self.vcf.write(
            "##INFO=<ID=ENTROPYALT,Number=A,Type=Float,Description=\"Shannon entropy of alt allele if longer than 100 bp\">\n")
        self.vcf.write(
            "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
        self.vcf.write(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
        self.vcf.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
        self.vcf.write(
            "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
        self.vcf.write(
            "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n")
        self.vcf.write(
            "##INFO=<ID=TRANSPOSE,Number=1,Type=String,Description=\"Transposon FASTA sequence ID that this breakpoint anchor matches\">\n")
        self.vcf.write(
            "##INFO=<ID=SAS,Number=1,Type=Float,Description=\"Simple alignment score. Likelihood a breakend is not structural but rather aligns to reference simply\">\n")
        self.vcf.write("##INFO=<ID=FW,Number=A,Type=Float,Description=\"Percent of forward reads\">\n")
        self.vcf.write("##INFO=<ID=BQ,Number=A,Type=Integer,Description=\"Average base quality at this position\">\n")
        self.vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        self.vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample Depth\">\n")
        self.vcf.write(
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
        self.vcf.write("##FORMAT=<ID=ED,Number=1,Type=Integer,Description=\"Edit distance\">\n")
        self.vcf.write("##FORMAT=<ID=OV,Number=A,Type=Integer,Description=\"Minimum read overlap in assembly\">\n")
        self.vcf.write(
            "##FILTER=<ID=homologous_breakends,Description=\"The edit_distance between sides of breakpoints was below the minimum allowed threshold\">\n")
        self.vcf.write(
            "##FILTER=<ID=too_many_alleles,Description=\"The set of possible alleles was too large/supported to be called\">\n")
        self.vcf.write("##FILTER=<ID=dust_mask,Description=\"At least 45 bp were considered masked out by DUST\">\n")
        self.vcf.write(
            "##FILTER=<ID=missing_assembly,Description=\"No hits were reported by BLAST for this assembly\">\n")
        self.vcf.write(
            "##FILTER=<ID=non_structural_alignment,Description=\"BLAST alignment indicates probable SNP\">\n")
        self.vcf.write("##FILTER=<ID=missing_anchor,Description=\"One or more anchors were not found by BLAST\">\n")
        self.vcf.write(
            "##FILTER=<ID=no_unique_anchor,Description=\"BLAST could not uniquely identify at least one of the anchors for this assembly\">\n")
        self.vcf.write(
            "##FILTER=<ID=ambiguous_anchor,Description=\"A BLAST query for this anchor reported multiple ambiguous hits\">\n")
        self.vcf.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
        self.vcf.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
        self.vcf.write("##contig=<ID=1,length=249250621>\n")
        self.vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

    def test_cut(self):
        """
        ATCAA-G-CACTA
        """
        test_spl = GenomeGraph(self.fake_reference)
        test_spl.split_ref_pos(self.chrom, 5)
        src = test_spl.get_node(self.chrom, 0)
        snk = test_spl.get_node(self.chrom, len(self.ref_seq))
        paths = test_spl.all_paths(src, snk)
        cnt = 0
        for p in paths:
            cnt += 1
            self.assertEqual(test_spl.get_path_seq(p), self.ref_seq)
        self.assertTrue(cnt, 1)
        self.assertEqual(test_spl.graph.number_of_nodes(), 4)
        self.assertEqual(len([x for x in test_spl.graph.nodes() if x.strand]), 2)

        test_spl.split_ref_pos(self.chrom, 6)
        src = test_spl.index[0]
        snk = test_spl.index[-1]
        paths = test_spl.all_paths(src, snk)
        cnt = 0
        for p in paths:
            cnt += 1
            self.assertEqual(test_spl.get_path_seq(p), self.ref_seq)
        self.assertEqual(test_spl.graph.number_of_nodes(), 6)
        self.assertEqual(len([x for x in test_spl.graph.nodes() if x.strand]), 3)
        self.assertTrue(cnt, 2)

    def test_snp(self):
        r"""
        ATCAA-G-CACTA
             \T/
        """
        # snp
        alt_seq = self.ref_seq[:5] + "T" + self.ref_seq[6:]
        test_snp = GenomeGraph(self.fake_reference)
        test_snp.add_var(self.chrom, 5, 6, 'T')

        src = test_snp.get_node(self.chrom, 0)
        snk = test_snp.get_node(self.chrom, len(self.ref_seq))
        paths = test_snp.all_paths(src, snk)
        cnt = 0
        for p in paths:
            cnt += 1
            self.assertTrue(test_snp.get_path_seq(p) in [self.ref_seq, alt_seq])
        self.assertEqual(cnt, 2)

    def test_del(self):
        r"""
        ATCAA-G-CACTA
             \_/
        """
        alt_seq = self.ref_seq[:5] + self.ref_seq[6:]
        test_del = GenomeGraph(self.fake_reference)
        test_del.add_var(self.chrom, 5, 6, "")
        src = test_del.get_node(self.chrom, 0)
        snk = test_del.get_node(self.chrom, len(self.ref_seq))
        cnt = 0
        for p in test_del.all_paths(src, snk):
            cnt += 1
            self.assertTrue(test_del.get_path_seq(p) in [self.ref_seq, alt_seq])
        self.assertEqual(cnt, 2)

    def test_del2(self):
        r"""
        ATCAA-GC-ACTA
             \__/
        """
        alt_seq = self.ref_seq[:5] + self.ref_seq[7:]
        test_del = GenomeGraph(self.fake_reference)
        test_del.add_var(self.chrom, 5, 7, "")
        src = test_del.get_node(self.chrom, 0)
        snk = test_del.get_node(self.chrom, len(self.ref_seq))
        cnt = 0
        for p in test_del.all_paths(src, snk):
            cnt += 1
            self.assertTrue(test_del.get_path_seq(p) in [self.ref_seq, alt_seq])
        self.assertEqual(cnt, 2)

    def test_repl(self):
        r"""
        ATCAA-GC-ACTA
             \T_/
        """
        allele = "T"
        alt_seq = self.ref_seq[:5] + allele + self.ref_seq[7:]
        test_del = GenomeGraph(self.fake_reference)
        test_del.add_var(self.chrom, 5, 7, allele)
        src = test_del.get_node(self.chrom, 0)
        snk = test_del.get_node(self.chrom, len(self.ref_seq))
        cnt = 0
        for p in test_del.all_paths(src, snk):
            cnt += 1
            self.assertTrue(test_del.get_path_seq(p) in [self.ref_seq, alt_seq])
        self.assertEqual(cnt, 2)

    def test_ins(self):
        r"""
        ATCAAG---CACTA
             \CAT/
        """
        allele = "CAT"
        alt_seq = self.ref_seq[:6] + allele + self.ref_seq[6:]
        test_ins = GenomeGraph(self.fake_reference)
        test_ins.add_var(self.chrom, 6, 6, allele)
        src = test_ins.get_node(self.chrom, 0)
        snk = test_ins.get_node(self.chrom, len(self.ref_seq))
        paths = test_ins.all_paths(src, snk)
        cnt = 0
        for p in paths:
            cnt += 1
            self.assertTrue(test_ins.get_path_seq(p) in [self.ref_seq, alt_seq])
        self.assertEqual(cnt, 2)

    def test_dbl1(self):
        r"""
        Testing adding two events
        ATCAAGCACTA
            \T/\G/
        I expect 4 alleles returned
        (all paths)

        I'll also test getting 'in the middle' src/snk
            e.g. from ref_seq[6] and ref_seq[6:]
        """
        exp1 = [self.ref_seq, self.ref_seq[:5] + "T" + self.ref_seq[6:]]
        exp2 = exp1 + [exp1[-1][:9] + "G" + exp1[-1][10:], self.ref_seq[:9] + "G" + self.ref_seq[10:]]
        test_dbl = GenomeGraph(self.fake_reference)
        test_dbl.add_var(self.chrom, 5, 6, "T")
        src = test_dbl.get_node(self.chrom, 0)
        snk = test_dbl.get_node(self.chrom, len(self.ref_seq))
        for path in test_dbl.all_paths(src, snk):
            seq = test_dbl.get_path_seq(path)
            self.assertTrue(seq in exp1)
            exp1.remove(seq)
        self.assertEqual(len(exp1), 0)

        test_dbl.add_var(self.chrom, 9, 10, "G")
        src = test_dbl.get_node(self.chrom, 0)
        snk = test_dbl.get_node(self.chrom, len(self.ref_seq))
        for path in test_dbl.all_paths(src, snk):
            seq = test_dbl.get_path_seq(path)
            self.assertTrue(seq in exp2)
            exp2.remove(seq)
        self.assertEqual(len(exp2), 0)

    def test_dbl_ovl(self):
        r"""
        Testing adding two events -- But they OVERLAP!
             /T\
        ATCAA-G-CACTA
             \_/

        I expect 3 alleles returned
        (all paths)

        When I delete 5,6 (or replace), my ref_up_node
        is the G. it needs to be the A
            I want to say an anchor base helps this, but
            what if that anchor base...? nope.
        when this base is being altered, and it's the entire node(?)
            then I
            if start != abs_start for md and new I just made..

        I'll also test getting 'in the middle' src/snk
            e.g. from ref_seq[6] and ref_seq[6:]
        """
        exp1 = [self.ref_seq, self.ref_seq[:5] + "T" + self.ref_seq[6:],
                self.ref_seq[:5] + self.ref_seq[6:]]

        test = GenomeGraph(self.fake_reference)
        test.add_var(self.chrom, 5, 6, "T")
        test.add_var(self.chrom, 5, 6, "")
        src = test.get_node(self.chrom, 0)  # shortcut this method?
        snk = test.get_node(self.chrom, len(self.ref_seq))
        for path in test.all_paths(src, snk):
            seq = test.get_path_seq(path)
            self.assertTrue(seq in exp1)
            exp1.remove(seq)
        self.assertEqual(len(exp1), 0)

    def test_dbl_ovl2(self):
        r"""
        Testing adding two events -- Overlap contained
              /T\
        ATC-AA-G-CA-CTA
           \_______/

        I expect 3 alleles returned
        (all paths)
        I'm losing the 4th base.. also test_dbl errors
        """
        exp1 = [self.ref_seq, self.ref_seq[:5] + "T" + self.ref_seq[6:],
                self.ref_seq[:4] + self.ref_seq[8:]]

        test = GenomeGraph(self.fake_reference)
        test.add_var(self.chrom, 5, 6, "T")
        test.add_var(self.chrom, 4, 8, "")
        src = test.get_node(self.chrom, 0)  # shortcut this method? could just test.index[0] and -1
        snk = test.get_node(self.chrom, len(self.ref_seq))
        for path in test.all_paths(src, snk):
            seq = test.get_path_seq(path)
            self.assertTrue(seq in exp1)
            exp1.remove(seq)
        self.assertEqual(len(exp1), 0)

    def test_paths(self):
        """
        check that all/ref/alt paths returns correctly
        """
        testg = GenomeGraph(FakeReference('f', "ATC"))
        testg.add_var('f', 1, 2, "G")
        cnt = 0
        for _ in testg.all_paths(testg.get_node('f', 0), testg.get_node('f', 5)):
            cnt += 1
        self.assertEqual(cnt, 2)

        cnt = 0
        for _ in testg.ref_path(testg.get_node('f', 0), testg.get_node('f', 5)):
            cnt += 1
        self.assertEqual(cnt, 1)

        cnt = 0
        for _ in testg.alt_paths(testg.get_node('f', 0), testg.get_node('f', 5)):
            cnt += 1
        self.assertEqual(cnt, 1)

    def test_fail_badref(self):
        """
        Need to fail when we give a bad chromosome or a bad coordinate within the refs
        """


    def test_index(self):
        """
        Do we approprately add to and access the reference
        """
        ref = FakeReference("chr1", "ATA")
        ref.ref_dict["chr2"] = "TAT"
        test = GenomeGraph(ref)
        test.add_var('chr1', 1, 2, "G")
        test.add_var('chr2', 1, 2, "C")
        for i in test.index:
            self.assertEqual(i, test.index[test.get_node_idx(i)])

    def test_multi_ref(self):
        """
        Can we correctly edit two references?
        """
        ref = FakeReference("chr1", "ATA")
        ref.ref_dict["chr2"] = "TAT"
        test = GenomeGraph(ref)

        test.add_var('chr1', 1, 2, "G")
        test.add_var('chr2', 1, 2, "C")
        rseq1 = ["ATA", "AGA"]
        rseq2 = ["TAT", "TCT"]

        for i in test.all_paths(test.get_node('chr1', 0), test.get_node('chr1', 5)):
            p = test.get_path_seq(i)
            self.assertTrue(p in rseq1)
            rseq1.remove(p)

        for i in test.all_paths(test.get_node('chr2', 0), test.get_node('chr2', 5)):
            p = test.get_path_seq(i)
            self.assertTrue(p in rseq2)
            rseq2.remove(p)

    def test_abs_shifted(self):
        r"""
        Do we properly add variation when insert sub-sequences of graphs
        this is important for kmer-analysis
            /---------\
        ATCA-AGC-A-CTA-AT----CAA-G-CACTA
                \T/      \TC/   \T/
        0123 456 7 890 12    345 6 78901
        0            1                2

        all paths:
        ATCA-------ATTCCAAGCACTA
        ATCA-------ATTCCAATCACTA
        ATCA-------AT--CAAGCACTA
        ATCA-------AT--CAATCACTA
        ATCAAGCACTAATTCCAAGCACTA
        ATCAAGCACTAATTCCAATCACTA
        ATCAAGCACTAAT--CAAGCACTA < ref
        ATCAAGCACTAAT--CAATCACTA
        ATCAAGCTCTAATTCCAAGCACTA
        ATCAAGCTCTAATTCCAATCACTA
        ATCAAGCTCTAAT--CAAGCACTA
        ATCAAGCTCTAAT--CAATCACTA
               *          *
        """
        exp = [x.replace('-', '') for x in ["ATCA-------ATTCCAAGCACTA",
                                            "ATCA-------ATTCCAATCACTA",
                                            "ATCA-------AT--CAAGCACTA",
                                            "ATCA-------AT--CAATCACTA",
                                            "ATCAAGCACTAATTCCAAGCACTA",
                                            "ATCAAGCACTAATTCCAATCACTA",
                                            "ATCAAGCACTAAT--CAAGCACTA",
                                            "ATCAAGCACTAAT--CAATCACTA",
                                            "ATCAAGCTCTAATTCCAAGCACTA",
                                            "ATCAAGCTCTAATTCCAATCACTA",
                                            "ATCAAGCTCTAAT--CAAGCACTA",
                                            "ATCAAGCTCTAAT--CAATCACTA"]]
        ref_seq = "ATCAAGCACTAATCAAGCACTA"

        test = GenomeGraph(FakeReference('crazy', ref_seq))
        offset = 0  # we don't offset anymore
        test.add_var('crazy', 4 + offset, 11 + offset, "")
        test.add_var('crazy', 7 + offset, 8 + offset, "T")
        test.add_var('crazy', 13 + offset, 13 + offset, "TC")
        test.add_var('crazy', 16 + offset, 17 + offset, "T")

        src = test.get_node('crazy', 0)  # shortcut this method? could just test.index[0] and -1
        snk = test.get_node('crazy', len(ref_seq) + offset)
        for path in test.all_paths(src, snk):
            seq = test.get_path_seq(path)
            self.assertTrue(seq in exp, msg="unk_seq = {0}".format(seq))
            exp.remove(seq)
        self.assertEqual(len(exp), 0)

    def test_vcftesting(self):
        """
        Let's make sure this vcf testing thing I'm making works
        """
        self.vcf.write("1\t10\t.\tA\tG\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        cnt = 0
        for entry in my_vcf:
            cnt += 1
            logging.debug(str(entry))
        self.assertEqual(cnt, 1)

    def test_vcf_snp(self):
        """
        Tests a vcf entry's SNP
        ATCACTA - ref
        ATCGCTA - alt
        """
        expected = ["ATCACTA", "ATCGCTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write("1\t4\t.\tA\tG\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_clumped(self):
        """
        Tests a vcf entry with multiple SNPs in the alt - almost a MNP. repl is my term for it
        ATCATACTA - ref
        ATCGTGCTA - alt
        """
        expected = ["ATCATACTA", "ATCGTGCTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(
            "1\t4\t.\tATA\tGTG\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_insertion(self):
        """
        Tests a fully resolved insertion
        ATCATACTA - ref
        ATCAGGGTACTA - alt
        """
        expected = ["ATCATACTA", "ATCAGGGTACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(
            "1\t4\t.\tA\tAGGG\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_deletion(self):
        """
        Tests a fully resolved deletion
        ATCAGGGTACTA - ref
        ATCATACTA - alt
        """
        expected = ["ATCAGGGTACTA", "ATCATACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(
            "1\t4\t.\tAGGG\tA\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_del2(self):
        """
        The other representation of a deletion
        """
        expected = ["ATCAGGGTACTA", "ATCATACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(("1\t4\tsv_3212\tA\t<DEL>\t100\thomologous_breakends\t"
                        "NS=1;DP=12;SVTYPE=DEL;END=7;SVLEN=-1761;AID=3161131;"
                        "IMPRECISE;CIPOS=0,0;CIEND=0,0\tGT:DP:AD:ED:OV\t./.:12:3,9:17:86\n"))
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_ins_repl(self):
        """
        Tests a fully resolved insertion but with some of the ref sequence deleted
        """
        expected = ["ATCATACTA", "ATCAGGGACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(
            "1\t4\t.\tAT\tAGGG\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_del_repl(self):
        """
        Tests a fully resolved deletion, but it has some extra ref sequence (this is the
        same as ins_repl, but with max size between ALT and REF switched
        """
        expected = ["ATCAGGGTACTA", "ATCACTACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(
            "1\t4\t.\tAGGG\tAC\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_bnd(self):
        """
        BND test 1 of breakpoints where there is no strand switch
        #s   t[p[ piece extending to the right of p is joined after t
        #s   ]p]t piece extending to the left of p is joined before t
        """
        expected = ["ATCAGGGTACTA", "ATCATACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(("1\t4\t.\tA\tA[1:7[\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))
        self.vcf.write(("1\t7\t.\tA\t]1:4]A\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))

        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)

        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_bnd2(self):
        """
        BND test 2 of breakpoints where there is no strand switch
        #s   t[p[ piece extending to the right of p is joined after t
        #s   ]p]t piece extending to the left of p is joined before t
        We'll be putting in a repl
        """
        expected = ["ATCAGGGTACTA", "ATCACTACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(("1\t4\t.\tA\tAC[1:7[\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))
        self.vcf.write(("1\t7\t.\tA\t]1:4]AC\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))

        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)

        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_bnd3(self):
        r"""
        BND test 3 of breakpoints where there is strand switch
        Full inversion
        # s   t]p] reverse comp piece extending left of p is joined after t
        # s   [p[t reverse comp piece extending right of p is joined before t


        G-CAAAA-GC
         C   <-^
         \  C<-/
        """
        expected = ["GCAAAAGC", "GCTTTTGC"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(("1\t2\t.\tC\tC]1:7]\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))
        self.vcf.write(("1\t7\t.\tG\t[1:2[G\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))

        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        for entry in my_vcf:
            ggraph.add_vcf(entry)
        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 1)

    def test_vcf_bnd4(self):
        r"""
        BND test 3 of breakpoints where there is strand switch
        Full inversion
        # s   t]p] reverse comp piece extending left of p is joined after t
        # s   [p[t reverse comp piece extending right of p is joined before t


        G-CAAAA-GC
         C   <-^
         \  C<-/
        """
        expected = ["CGATGCAGGCTAGATCGAT", "CGATGCctagcctATCGAT".upper(),
                    #       *                          *
                    "CGATGCATGCTAGATCGAT", "CGATGCctagcatATCGAT".upper()]
        fake_ref = FakeReference('1', expected[0])
        ggraph = GenomeGraph(fake_ref)

        self.vcf.write(("1\t6\t.\tC\tC]1:14]\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))
        self.vcf.write("1\t8\t.\tG\tT\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        self.vcf.write(("1\t14\t.\tA\t[1:6[A\t100\tPASS\tNS=1;DP=906;SVTYPE=BND;AID=3043759;MATEID=bnd_6087519"
                        "\tGT:DP:AD:ED:OV\t./.:907:.,907:16:100\n"))
        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)

        for entry in my_vcf:
            ggraph.add_vcf(entry)

        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)

    def test_vcf_node_data(self):
        """
        add variants and then recover the edge-data (vcf_records)
        """
        # can I add vcf_record in edge_data and some variants (even same edge on some of them)
        # and recover the edge_data correctly - mainly worried about edge_data.append(each initital):w
        expected = ["ATCAGGGTACTA",
                    #   *     *
                    "ATCTGGGTAGTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))
        # back
        self.vcf.write(
            ("1\t4\t.\tA\tT\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n"
             "1\t10\t.\tC\tG\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n")
        )

        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        my_entries = []
        for entry in my_vcf:
            ggraph.add_vcf(entry)
            my_entries.append(entry)
        alt_node_cnt = 0
        for node in ggraph.alt_node_iter():
            if not node.strand:
                continue
            alt_node_cnt += 1
            self.assertTrue(node.var_data in my_entries)
            my_entries.remove(node.var_data)
        self.assertEqual(alt_node_cnt, 2)
        self.assertEqual(len(my_entries), 0)

    def test_vcf_edge_data_ovl(self):
        """
        add variants and then recover the edge-data (vcf_records)
        some of the records will have the same overlap
        """
        expected = ["ATCAGGGTACTA",
                    #   *
                    "ATCTGGGTACTA",
                    "ATCCGGGTACTA"]
        ggraph = GenomeGraph(FakeReference('1', expected[0]))

        self.vcf.write(
            "1\t4\t.\tA\tT,C\t100\tPASS\tNS=1;DP=32;AID=2993215;FW=0.37;BQ=65\tGT:DP:AD:OV\t0/1:32:7,25:87\n"
        )

        self.vcf.seek(0)
        my_vcf = vcf.Reader(self.vcf)
        my_entries = []
        for entry in my_vcf:
            ggraph.add_vcf(entry)
            my_entries.append(entry)

        alt_node_cnt = 0
        for j in ggraph.alt_edge_iter():
            alt_node_cnt += 1
            alt_edge_cnt = 0
            for _, snk, _ in j:
                if snk.is_alt:
                    self.assertTrue(snk.var_data in my_entries)
                    alt_edge_cnt += 1
            self.assertEqual(alt_edge_cnt, 2)

        self.assertEqual(alt_node_cnt, 1)

        paths = ggraph.all_paths(ggraph.index[0], ggraph.index[-1])
        for p in paths:
            seq = ggraph.get_path_seq(p)
            self.assertTrue(seq in expected)
            expected.remove(seq)
        self.assertEqual(len(expected), 0)
    # Test GNode sequence - ref and regular string.., with and without Ns
    # Test GenomeGraph.get_chr_end_node
    # Test index remove/add
    # Test get allele_seqs
    # Test get_node_kmer
    # Test get_edge_kmer

if __name__ == '__main__':
    unittest.main()
