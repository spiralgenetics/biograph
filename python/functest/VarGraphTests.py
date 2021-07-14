import os
import unittest

import vcf
import biograph as bgsdk


LAMBDIR="/home/english/lambda_testing/lambdaToyData"
class VarGraphTestCases(unittest.TestCase):
    """
    Test using the biograph.VarGraph
    """
    def setUp(self):
        self.ref = bgsdk.Reference(os.path.join(LAMBDIR, "benchmark/ref_lambda"))
        self.bg = bgsdk.BioGraph(os.path.join(LAMBDIR, "benchmark/father_lambda.bg"))
        self.variants = vcf.Reader(filename=os.path.join(LAMBDIR, "variants/family.vcf"))
    def print_nodes(self, nodes):
        ret = "Incorrect Genotype\nis_ref start end pair_upcov pair_dncov unpair_upcov unpair_dncov seq_cov\n"
        for nod in nodes:
            ntext = " ".join([str(x) for x in [nod.is_ref, nod.start, nod.end, 
                        nod.paired.upstream_edge, nod.paired.downstream_edge, 
                        nod.unpaired.upstream_edge, nod.unpaired.downstream_edge]])
            ret += ntext + "\n"
        return ret
    
    def test005_readmap_test(self):
        """
        I'm seeing something weird in the reads
        """
        b1 = bgsdk.BioGraph(os.path.join(LAMBDIR, "benchmark/father_lambda.bg"))
        self.bg.set_readmap("father")
        seq = "TATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATGTTTTTATTTAAAATACCCTCTGAAAAG"
        a = set([b1.read_to_entry(x).sequence for x in b1.read_search(seq)])
        b = set([self.bg.read_to_entry(x).sequence for x in self.bg.read_search(seq)])
        self.assertEqual(len(a - b), 0)

    def test010_fullvargraph_single(self):
        for sample in ["proband", "father", "mother"]:
            print "analyzing", sample
            self.bg.set_readmap(sample)
            #Only the one for lambda
            location = self.ref.supercontigs[0]
            varg = bgsdk.VarGraph(location, 10000)
            vlookup = {}
            tlookup = {}
            for entry in self.variants:
                #VarGraph over a section of the reference
                #I need to ensure that the variant doesn't span across contigs
                #They will, though
                #overload add_variant so it can take a sequence OR a string
                varg.add_variant(entry.start, entry.end, bgsdk.Sequence(str(entry.ALT[0])))
                vlookup[entry.start] = entry.genotype(sample)["GT"], entry.INFO["vtype"], str(entry.ALT[0]), entry.start, entry.end
    
            varg.trace(self.bg)
            cov = self.bg.seq_coverage(location.sequence)
            #would prefer to have an accessor to nodes[name]
            nodes = varg.get_nodes() #this is something
            for key in vlookup:
                if vlookup[key][1] == "INS":
                    continue
                #should only be a single ref, right?
                #Then multiple nodes.. 
                #I should organize them.. ugh that's hard
                for m_node in nodes[key]:
                    if m_node.is_ref:
                        ref_cov_u = m_node.unpaired.upstream_edge + m_node.unpaired.downstream_edge
                        ref_cov_p = m_node.paired.upstream_edge + m_node.paired.downstream_edge
                    else:
                        alt_cov_u = m_node.unpaired.upstream_edge + m_node.unpaired.downstream_edge
                        alt_cov_p = m_node.paired.upstream_edge + m_node.paired.downstream_edge
                #self.print_nodes(nodes[key])
                #print "variant",  tlookup[key]
                #print "expected", vlookup[key], "found", bgsdk.genotyper(ref_cov_u + alt_cov_u, alt_cov_u)[0]
                #print "expected", vlookup[key],"found", bgsdk.genotyper(ref_cov_p + alt_cov_p, alt_cov_p)[0]
                m_gt = bgsdk.genotyper(ref_cov_u + alt_cov_u, alt_cov_u)[0]
                self.assertEqual(vlookup[key][0], m_gt, self.print_nodes(nodes[key]) + ("%s != %s" % (str(vlookup[key]), m_gt)))
                m_gt = bgsdk.genotyper(ref_cov_p + alt_cov_p, alt_cov_p)[0]
                self.assertEqual(vlookup[key], m_gt, self.print_nodes(nodes[key]) + ("%s != %s" % (str(vlookup[key]), m_gt)))
            #have to reset the vcf
            self.variants = vcf.Reader(filename=os.path.join(LAMBDIR, "variants/family.vcf"))
            
        
        #subset trace
        #check the redundancy thing
    
    def test_020_partial_vargraph_single(self):
        """
        Test one individual tracing just a single variant
        """
        return
        sub_start = 16619-500
        sub_end = 16619+500
        location = self.ref.make_range("lambda", sub_start, sub_end)
        varg = bgsdk.VarGraph(location)
        for entry in self.variants:
            #VarGraph over a section of the reference
            #I need to ensure that the variant doesn't span across contigs
            #They will, though
            #overload add_variant so it can take a sequence OR a string
            if sub_start < entry.POS < sub_end:
                varg.add_variant(entry.start, entry.end, bgsdk.Sequence(str(entry.ALT[0])))
        #I guess 
        #total trace
        #TODO defaults
        varg.trace_sub(self.bg, sub_start, sub_end)
        cov = self.bg.seq_coverage(location.sequence)
        #would prefer to have an accessor to nodes[name]
        print "is_ref start end pair_upcov pair_dncov unpair_upcov unpair_dncov seq_cov"
        j = varg.get_nodes() #this is something
        for n in varg.get_nodes():
            s = cov[n.start:n.end]
            if len(s) == 0:
                mcov = 0
            else:
                mcov = sum(s)/float(len(s))
            print n.is_ref, n.start, n.end, n.paired.upstream_edge, n.paired.downstream_edge, \
                  n.unpaired.upstream_edge, n.unpaired.downstream_edge, mcov, s
            
 
        #varg.trace(b, sub_start, sub_end)
    
    def test_030_fullgraph_trio(self):
        """
        pass
        """
        pass
    #A different module will consume the coverages and make GTs
    #Along with joint genotyping...maybe.. Just a Genotyping module that consumes
    #a VCF and a BioGraph and optionally a pedigree - this is the pcmp tool

if __name__ == '__main__':
    unittest.main()
