import pysam
import biographsdk as bgsdk
import itertools
from string import maketrans
import pysam
from collections import Counter
import json

baseName = "/EnglishVolume/data/lambda_toys/weirdLambda/imperfect"
bg = bgsdk.Biograph(baseName + ".seqset", baseName + ".readmap")
print("truth for someshort = 13697 / 13703 reads & 1369695 bases")
print("Number of reads in biograph: %d" % bg.num_reads())
print("Number of bases in biograph: %d" % bg.num_bases())

bam = pysam.Samfile(baseName + ".bam", 'rb')
num_bases = 0
num_reads = 0
for read in bam:
    num_bases += len(read.seq)
    num_reads += 1;

print("Number of reads in bam: %d" % num_reads)
print("Number of bases in bam: %d" % num_bases)

#19 occurances of this kmer
nineteen = "ACATACATGTGCCTGCAAACAGTACC"
print("Finding the exact kmer")
x = bg.find(nineteen)
print(x)

print("Finding kmer with some wiggle room")
for i in range(3):
    print(i)
    y = bg.find_near(nineteen, i)
    print(y)
print(bg.find_near(nineteen))
print("Finding reads with the exact kmer")

"""
>> myIndiv = bg.samples["A"] # '09100-1094011
>> edge = bg.assemble(chr, start, end).edges[0]
>> myLookup = [bg.samples["A"], bg.samples["B"]]

>> edge.coverage()
12
>> edge.coverage(myIndiv)
1
>> edge.coverage(bg.samples.values()[:3])
3
>> edge.coverage(myLookup)
2

for key in bg.metadata.samples.vlaues():
    bg.readmaps[key].anfoinawef

my_interst = bg.unite_samples('a')

z = bg.find_reads(kmer, 1000, mask=bg.metada['a'])

z.extedn(bg.findrawe, indv=B)
z bg.fore_call( 
print(z)
print(len(z))
"""

print("Finding all permutations upto 4mers")
j = []
for i in range(1, 3):
    j.extend(itertools.combinations_with_replacement("ATCG", i))

j.sort()
print("kmer\tctx")
import sys
for i in j:
    i = "".join(i)
    ctx = bg.find(i)
    
    print("%s\t%s" % (i, ctx))

revcmp = maketrans("ATCG", "TAGC")
def orientation_test(kmer):
    fwd = kmer
    rev = fwd.translate(revcmp)[::-1]
    
    print("orientation test on %s, %s" % (fwd, rev) )
    
    bgFcnt = bg.find_reads(fwd, 1000)
    bgRcnt = bg.find_reads(rev, 1000)
    bgFnearCnt = len(bg.find_reads_near(fwd, 3, 1000))
    bgRnearCnt = len(bg.find_reads_near(rev, 3, 1000))
    bamFcnt, bamRcnt, bothCnt = 0, 0, 0
    bamFncnt, bamRncnt = 0, 0
    
    posCnt = Counter()
    bam = pysam.Samfile(baseName + '.bam', 'rb')
    for read in bam:
        if read.is_reverse:
            query = read.query.translate(revcmp)[::-1]
        else:
            query = read.query
        
        if fwd in query:
            bamFcnt += 1
            posCnt[query.index(fwd)] += 1
        
        if rev in query:
            bamRcnt += 1
            p = query.index(rev)
            posCnt[query.index(rev)] += 1
        
        if fwd in query and rev in query:
            bothCnt += 1
        
        bamFncnt += query.count(fwd)
        bamRncnt += query.count(rev)
        
    x = set([str(i.sequence) for i in bgFcnt])
    y = set([str(i.sequence) for i in bgRcnt])
    z = x.intersection(y)
    print("bg fwd reads (unq): %d (%d)" % (len(bgFcnt), len(x)))
    print("bg rev reads (unq): %d (%d)" % (len(bgRcnt), len(y)))
    print("bg intersection reads: %d" % (len(x.intersection(y))))
    print("bg union reads: %d" % (len(x.union(y))))
    print("bg reads near fwd (3): %d" % (bgFnearCnt))
    print("bg reads near rev (3): %d" % (bgRnearCnt))
    print("bam fwd reads containing: %d" % bamFcnt)
    print("bam rev reads containing: %d" % bamRcnt)
    print("bam fwd reads count: %d" % bamFncnt)
    print("bam rev reads count: %d" % bamRncnt)
    print("bam fwd+rev reads: %d" % bothCnt)
    print("Positions: %s" % (json.dumps(posCnt)))
    print("sum poss: %s" % (sum(posCnt.values())))
    #for i in bgFcnt:
        #print("bbbbf %s" % i)
    #for i in bgRcnt:
        #print("bbbbr %s" % i)
    
print("same cnt test")
orientation_test("AAACATAGCA")
print("diff cnt test")
orientation_test('AATAATTGCA')
                        
x = bg.find_reads_near(nineteen, 2, 1000)
print(len(x))

##trying to do an assembly

#fwd = bgsdk.find_anchors
#results = bgsdk.assemble(fwd, rev, min_overlap, max_steps, skip_ambiguous, read_map)

def multiPosTest(kmer):
    """
    Does this kmer have multiple reads possesing the kmer in the same position
    """
    fwd = kmer
    rev = fwd.translate(revcmp)[::-1]
    bgFcnt = bg.find_reads(fwd, 1000)
    bgRcnt = bg.find_reads(rev, 1000)
    bamFcnt, bamRcnt, bothCnt = 0, 0, 0
    bamFncnt, bamRncnt = 0, 0
    
    posCnt = Counter()
    bam = pysam.Samfile(baseName + '.bam', 'rb')
    for read in bam:
        if read.is_reverse:
            query = read.query.translate(revcmp)[::-1]
        else:
            query = read.query
        
        if fwd in query:
            bamFcnt += 1
            posCnt[query.index(fwd)] += 1
        
        if rev in query:
            bamRcnt += 1
            p = query.index(rev)
            posCnt[query.index(rev)] += 1
        
        if fwd in query and rev in query:
            bothCnt += 1
        
        bamFncnt += query.count(fwd)
        bamRncnt += query.count(rev)
    return (len(bgFcnt) == len(bgRcnt), len([x for x in posCnt if posCnt[x] > 1])) 

#I didn't go all the way through this huge brute force, but 
#I didn't see anythin for the first 20~minutes come through.
#It's possible that this redundant read fetch problem doesn't happen when
#reads don't have reads in the same
#with open("/share/functest/libspiralapi/toydata/lambda/ref_weirdLambda/source.fasta", 'r') as fh:
    #fh.readline()
    #ref = "".join(fh.readlines()).replace('\n','')
    #for i in itertools.product("ATCG", repeat=10):
        #i = "".join(i)
        #if ref.count(i) ==  1:
            #print(i, ref.count(i), multiPosTest(i))
        #if not same and cnt == 0:
            #print(i, same, cnt) #multiPosTest(i))

def weird_test():
    print("Finding all read 20mers present")
    
    print("kmer\tnum_reads")
    rev = maketrans("ATCG", "TAGC")
    for i in itertools.product("ATCG", repeat=10):
        i = "".join(i)
        rc = i.translate(rev)[::-1]
        reads = bg.find_reads(i, 1000)
        readsrc = bg.find_reads(rc, 1000)
        if len(reads) == 0 and len(readsrc) == 0:
            continue
        if len(reads) != len(readsrc):
            print(i, rc, len(reads), len(readsrc))
        #if len(reads) > 0:
                #print("%s\t%s" % (i, len(reads)))
