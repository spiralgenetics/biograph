
# coding: utf-8

# In[195]:

from biograph import BioGraph, Sequence
cseq = "GGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATG"
diff = "         dd        *      i"
qseq = "GGTTTAAGGTTTCCGTTTTTCTTCAGTCATAACTTAATG"
diff = "         dd        *      i        ****"
qseq = "GGTTTAAGGTTTCCGTTTTTCTTCAGTCATAACTTTTTT"
qseq = "GGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTCCCC"
diff = "                                    ****"
qseq = "acattttaaacaatttatactagccccccccgcccttttttaacagctctcataaatatgtagttcatttctctTgaatacttaagtcaactgtatga".upper()
diff = "acattttaaacaatttatactagcccccccc..gcccttttttaacagctctcataaatatgtagttcatttctcttgaatacttaagtcaactgtatga"
bases = ["A", "T", "C", "G"]
bg_file = "/home/english/temp/lambdaToyData/proband.bg"

bg_file = "/home/english/ajtrio/biographs/HG002-NA24385-50x.bg/"
max_diff = 3
find1 = "ACATTTTAAACAATTTATACTAGCCCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTTGAATACTTAAGTCAACTGTATGA"
find2 = "ACATTTTAAACAATTTATACTAGCCCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTGGAATACTTAAGTCAACTGTATGA"
""""
>ref
ACATTTTAAACAATTTATACTAGCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTTGAATACTTAAGTCAACTGTATGA
>query1
ACATTTTAAACAATTTATACTAGCCCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTTGAATACTTAAGTCAACTGTATGA
>query2
ACATTTTAAACAATTTATACTAGCCCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTGGAATACTTAAGTCAACTGTATGA

   ACTIONS      QUERY           SCORE START  END QSIZE IDENTITY CHRO STRAND  START    END      SPAN
   ---------------------------------------------------------------------------------------------------
   browser details query1            97     1   100   100  96.0%    17   +   57832317  57832414     98
   browser details query2            95     1   100   100  94.9%    17   +   57832317  57832414     98
   browser details ref               98     1    98    98 100.0%    17   +   57832317  57832414     98

query1 - details
00000001 acattttaaacaatttatactagccccccccccgcccttttttaacagct 00000050
>>>>>>>> |||||||||||||||||||||||||||||||  ||||||||||||||||| >>>>>>>>
57832317 acattttaaacaatttatactagcccccccc..gcccttttttaacagct 57832364
                                        ^
                                    57832348
00000051 ctcataaatatgtagttcatttctcttgaatacttaagtcaactgtatga 00000100
>>>>>>>> |||||||||||||||||||||||||||||||||||||||||||||||||| >>>>>>>>
57832365 ctcataaatatgtagttcatttctcttgaatacttaagtcaactgtatga 57832414

query2 - details
00000001 acattttaaacaatttatactagccccccccccgcccttttttaacagct 00000050
>>>>>>>> |||||||||||||||||||||||||||||||  ||||||||||||||||| >>>>>>>>
57832317 acattttaaacaatttatactagcccccccc..gcccttttttaacagct 57832364
                                        ^
                                    57832348
00000051 ctcataaatatgtagttcatttctctggaatacttaagtcaactgtatga 00000100
>>>>>>>> |||||||||||||||||||||||||| ||||||||||||||||||||||| >>>>>>>>
57832365 ctcataaatatgtagttcatttctcttgaatacttaagtcaactgtatga 57832414
                                   ^
                                57832392
Both variants exist!
https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A57832317-57832414&hgsid=596652695_gh1Wcz1EruLV7qLCcogBY3CdjlV9
"""
# hotcache
import sys
# sys.stderr.write("hotcache")
# with open(bg_file + "/seqset", "rb") as fh:
    # for chunk in iter(lambda: fh.read(4096), b""):
        # pass
# sys.stderr.write(" fin\n")
# In[196]:

my_bg = BioGraph(bg_file)
print dir(my_bg)
print my_bg.entry_search(qseq)

# try to put on every base
# if ctx becomes invalid, or it uses up errors before the end
#   pop off and try an alternate


class mtemp(object):

    def __init__(self, ctx, err):
        self.ctx = ctx
        self.err = err
        self.key = "%d_%d" % (self.ctx.start, self.ctx.end)
        self.mhash = int(self.key.__hash__())

    def __hash__(self):
        return self.mhash

    def __eq__(self, other):
        return self.key == other.key

MAXOPS = 500000
MAXOPS = 20000


def fuzz_find(ctx, remain, error=0, max_error=5, max_results=2, results=None, max_depth=100, pindel=0):
    """
    returns [Found a match, [matches]]
    pindel = -1 if previously deleted
    pindel = 1 if previously deleted
    """
    global MAXOPS
    # print pindel, ctx
    MAXOPS -= 1
    if results is None:
        results = set()

    if MAXOPS <= 0:
        # print remain, 'max'
        return False, []

    # don't search anymore
    if len(results) >= max_results:
        return False, []

    # quit on this?
    if not ctx.valid:
        return False, []

    # print MAXOPS, ctx, len(ctx.sequence), len(remain), error, len(results)
    # only matches can save you now
    if error >= max_error:
        pos = 0
        while ctx.valid and pos < len(remain):
            ctx = ctx.push_front(remain[pos])
            pos += 1
        if ctx.valid and pos == len(remain):
            return [True, [mtemp(ctx, error)]]
        return [pos == len(remain), []]

    # print len(remain), len(ctx.sequence), len(results)
    # we're finished
    if len(remain) == 0:
        return [True, [mtemp(ctx, error)]]

    # match
    new = fuzz_find(ctx.push_front(remain[0]), remain[
                    1:], error, max_error, max_results, results, max_depth - 1, pindel=0)

    # we've found a match and it won't get any better because this ctx has a single path
    # so when we pop all the way back,
    # it'll be at a context with >1 sequences, and we'll explore the next possibility
    if new[0] and ctx.end - ctx.start == 1:
        results.update(new[1])
        return [False, results]

    # mis
    any_match = False
    for i in [x for x in bases if x != remain[0]]:
        new = fuzz_find(ctx.push_front(i), remain[1:], 
                    error + 1, max_error, max_results, results, max_depth - 1, pindel=0)
        if new[0]:
            results.update(new[1])
            any_match = True
    if any_match and ctx.end - ctx.start == 1:
        return [False, results]

    # del in query so we'll look for an extra base
    if pindel != 1:
        any_match = False
        for i in [x for x in bases if x != remain[0]]:
            new = fuzz_find(ctx.push_front(i), remain, error + 1,
                            max_error, max_results, results, max_depth - 1, pindel=-1)
            # results.update(new[1])
            if new[0]:
                results.update(new[1])
                any_match = True
                # print len(results), max_depth
                # return new
        if any_match and ctx.end - ctx.start == 1:
            return [False, results]

    # ins in query so we'll query without this next base
    if pindel != -1:
        new = fuzz_find(ctx, remain[1:], error + 1, max_error, max_results, results, max_depth - 1, pindel=1)
        # results.update(new[1])
        if new[0] and ctx.end - ctx.start == 1:
            results.update(new[1])
            return [False, results]
    # print 'backup'
    return [True, results]

ctx = my_bg.entry_search("")
print fuzz_find(ctx[0], find1[::-1])
print fuzz_find(ctx[0], find2[::-1])
# flag, ret = search_further(ctx, qseq[::-1], 0, max_results=5)
# uniquify
"""
Remove redundant candidates.
"""
# for j in ret:
    # print j.ctx.sequence, j.ctx, j.err

# print len(my_bg.find_reads(find1))
# print len(my_bg.find_reads(find2))
