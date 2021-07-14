""" First test - not successful"""
from biograph import *
import datetime
r = Reference("/reference/hs37d5/")
chrom = "4"

p1 = datetime.datetime.now()
b = BioGraph("/home/english/ajtrio/biographs/HG002_v3Biograph.bg/", CacheStrategy.MMAPCACHE)
p2 = datetime.datetime.now()
print p2-p1

def add_vcf(self, vcf_record):
    for allele in vcf_record.alleles[1:]:
        if allele.type == "BND":
            self.add_bnd(vcf_record, allele)
            continue
        if allele.type == "DEL":#FOR <DEL> alts
            if "END" not in vcf_record.INFO:
                raise KeyError("END not found in vcf entry's INFO")
            end = vcf_record.INFO["END"]
            alt_seq = vcf_record.REF
        else:
            end = vcf_record.end
            alt_seq = allele.sequence
        self.add_var(vcf_record.CHROM, vcf_record.start, end,
                        alt_seq, var_data=vcf_record)

#Okay - the above might work.. except for vt normalized <DEL>
#anchor base and all that jazz puts it back
ref_seq = "CCAGGGAAGCATTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTATACTCTAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAATGTGTCATCTAGCATTAGGTATATCTCCCAATGCTATCCCTCCCCCCTCCCCCGACCCCACCACAGTCCCCAGAGTGTGATATTCCCCTTCCTGTGTCCATGTGATCTCATTGTTCAATTCCCACCTATGAGTGAGAATATGCGGTGTTTGGTTTTTTGTTCTTGCGATAGTTTACTGAGAATGATGGTTTCCAATTTCATCCATGTCCCTACAAAGGATATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACAGTCCCACCAACAGTGTAAAAGTGTTCCTATTTCTCCGCATCCTCTCCAGCACCTGTTGTTTCCTGACTTCATGTCCAAAACACCAAAAGCAATGGCAACAAAAGACAAAATTGACAAATGGGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCAGAGTGAACAGGCAACCTACAACATGGGAGAAAATTTTTGCAACCTACTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAACTCAAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCGAAGGACATGAACAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAAACACATGAAGAAATGCTCATCATCACTGGCCATCAGAGAAATGCAAATCAAAACCACTATGAGATATCATCTCACACCAGTTAGAATGGCAATCATTAAAAAGT"
start = 7995794 - 1
end = start + len(ref_seq)
seq = "C"
seq2 = r.make_range("4", start, end)
#print ref_seq == str(seq2.sequence)
#print end - start
#print ref_seq
print str(seq2.sequence)
ctg = r.get_supercontig(chrom, start)
v = VarGraph(ctg, 1000)
v.add_variant(start, end, Sequence(seq))

p1 = datetime.datetime.now()
v.trace_sub(b, start - 500, end + 500)
p2 = datetime.datetime.now()
print p2-p1

print "is_ref start end pair_upcov pair_dncov unpair_upcov unpair_dncov"
all_collapsed_nodes = v.get_nodes()
for x in all_collapsed_nodes:
    for n in x:
        print n.is_ref, n.start, n.end, n.paired.upstream_edge, n.paired.downstream_edge, n.unpaired.upstream_edge, n.unpaired.downstream_edge


#grab straight up number of reads
seq = r.make_range("4", end-25, end+25)

print "number of reads", len(b.read_search(seq.sequence))

seq = r.make_range("4", end-500, end+500)
print b.seq_coverage(seq.sequence)

