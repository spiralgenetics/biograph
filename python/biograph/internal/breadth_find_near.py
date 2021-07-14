
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

bg_file = "/home/english/ajtrio/HG002-NA24385-50x.bg/"
my_bg = BioGraph(bg_file)
max_diff = 3

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
sys.stderr.write("hotcache")
with open(bg_file + "/seqset", "rb") as fh:
    for chunk in iter(lambda: fh.read(4096), b""):
        pass
sys.stderr.write(" fin\n")
# In[196]:

print my_bg.find(qseq)


def find_diff(my_bg, search_seq, max_diff=3):
    """
    search for a string with up to max_diff mis/ins/del
    returns a list of seqset_ranges
    """
    search_seq = search_seq[::-1]  # Reverse it for pop_front
    candidates = []
    cur_diff = 0
    if max_diff == 0:
        ret = my_bg.find(search_seq)
        if ret.valid:
            return [ret]
        return []
    for base in bases:
        new = my_bg.find(base[0])
        if new.valid:
            # print base, cur_diff, search_seq[0], base == search_seq[0]
            if base == search_seq[0]:
                candidates.append((new, 0, search_seq[1:], False))  # mat
            else:
                candidates.append((new, 1, search_seq, False))  # ins err
                candidates.append((new, 1, search_seq[1:], True))  # mis err

    candidates.append((my_bg.find(""), 1, search_seq[1:], True))  # del err
    # print "start" redundancy filtering
    ret = find_diff_recursive(candidates, max_diff, len(search_seq))
    for i in ret:
        print i[0].sequence, i
    return candidate_unique(ret, False)


def candidate_unique(candidates, with_remain=True):
    """
    Remove redundant candidates.
    if with_remain: we're still building so we need to take remaining sequence into account
    """
    uhash = {}
    for i in candidates:
        if with_remain:
            key = "%d_%d_%s" % (i[0].start, i[0].end, i[2])
        else:
            key = "%d_%d" % (i[0].start, i[0].end)
        if key not in uhash:
            uhash[key] = i
        else:
            if uhash[key][1] > i[1]:
                uhash[key] = i
    print 'filt', len(candidates), '->', len(uhash)
    return uhash.values()


def find_diff_recursive(candidates, max_diff, query_len):
    """
    Candidate = (sequence being build, it's difference, sequence remaining to explore, previous base indel)
    #can't have side by side insertion/deletions
    """
    if len(candidates) == 0:
        return []
    # reduce redundancy -- hardly worth it?
    candidates = candidate_unique(candidates)

    new_candidates = []
    finished = []
    for cand, diff_cnt, remain_seq, pindel in candidates:  # for each of the candidates
        # print cand.sequence, diff_cnt, remain_seq
        if len(remain_seq) == 0:  # nothing left to do
            # print('completed')
            finished.append((cand, diff_cnt))
        elif diff_cnt >= max_diff:  # can only try matching
            # print("maxdiff'd")
            pos = 0
            while cand.valid and pos < len(remain_seq):
                cand = cand.push_front(remain_seq[pos])
                pos += 1
            if cand.valid and pos == len(remain_seq):
                # print('but good', pos, len(remain_seq), cand.sequence)
                finished.append((cand, diff_cnt))
        else:
            has_match = False
            for base in bases:  # try pushing on a new base
                new = cand.push_front(base[0])
                if new.valid:
                    if base == remain_seq[0]:
                        has_match = True
                        new_candidates.append((new, diff_cnt, remain_seq[1:], False))  # mat
                        # you can't ins or del a match
                        # left-aligned
                    else:
                        new_candidates.append((new, diff_cnt + 1, remain_seq[1:], False))  # mis
                        if not pindel:
                            new_candidates.append((new, diff_cnt + 1, remain_seq, True))  # ins

            # didn't find this base, give a candidate with a del
            if not has_match and not pindel:
                new_candidates.append((cand, diff_cnt + 1, remain_seq[1:], True))  # del
    return finished + find_diff_recursive(new_candidates, max_diff, query_len)


# In[199]:
import time
sys.stderr.write(time.asctime() + ' first\n')
j = find_diff(my_bg, qseq, 5)

sys.stderr.write(time.asctime() + ' second\n')
j = find_diff(my_bg, qseq, 5)
sys.stderr.write(time.asctime() + ' finish\n')

# In[200]:
print cseq
print diff
print qseq, len(j)
for x in j:
    print x[0].sequence, x


# In[ ]:
exit()


# In[136]:

def ne(self, other):
    return not self == other
Sequence.__ne__ = ne


def hash(self):
    return str(self).__hash__()
Sequence.__hash__ = hash


# In[137]:

print(len(j))


# In[135]:

"j".__hash__


# In[ ]:
