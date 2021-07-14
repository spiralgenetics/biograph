# coding: utf-8
import sys
import copy
import itertools
import logging
from StringIO import StringIO

from biograph import BioGraph, Sequence


class ContextState(object):

    """
    Holds a context, and it's state.
    State Properties:
        Operations to be performed
        Seen Paths (ctx_keys, score, qpos)

    This will allow us to:
        Resume query having explored part of the operations
        Short circuit when divergence gets us to somewhere we've already seen

    """

    def __init__(self, context, query, mask=None, error_cnt=0, pindel=0):
        self.context = context
        self.oquery = query
        self.query = buffer(query)
        self.cur_pos = 0
        self.error_cnt = error_cnt
        self.pindel = 0
        self.result_cnt = 0
        # Need to keep where we've been for back_tracking
        self.__query_stack = []
        self.__context_stack = []
        self.__hide_results = set()

    def get_results(self):
        """
        Return the results
        """
        return self.__hide_results

    def record_result(self):
        """
        Save the state of the current alignment and be done
        """
        # take this out for debugging the - you're back to the same
        self.__hide_results.add(self.context.sequence)

    def roll_back(self, to_step):
        """
        Restore to a previous context being the head.
        Remove from context stack and alignment where needed
        """
        logging.debug('rolling %d %d', to_step, len(self.__context_stack))
        self.query = self.__query_stack[to_step]
        self.__query_stack[:to_step]
        self.context = self.__context_stack[to_step]
        self.__context_stack = self.__context_stack[:to_step]

    # shortcuts to the alignment
    def is_aligned(self):
        """
        We've used up all the query bases
        """
        return len(self.query) == 0

    def next_base(self):
        """
        gimme the next base
        """
        return self.query[0]

    def next_mask(self):
        """
        Is the next  base masked
        """
        return False

    def add_match(self, base=None):
        """
        Align a match
        """
        if base is None:
            base = self.next_base()
        self.__context_stack.append(self.context)
        logging.debug('pushing %s', base)
        self.context = self.context.push_front(base)
        self.__query_stack.append(self.query)
        self.query = buffer(self.query, 1)

        # self.push(base)
        # self.alignment.add_match(base)

    def add_ins(self):
        """
        Insertion in query relative to sequence
        """
        self.__context_stack.append(self.context)
        self.__query_stack.append(self.query)

        # self.alignment.add_ins()

    def add_del(self, base):
        """
        Deletion in query relative to sequence
        """
        # self.push(base)
        self.__context_stack.append(self.context)
        self.context = self.context.push_front(base)
        self.__query_stack.append(self.query)
        self.query = buffer(self.query, 1)


class GAligner(object):

    """
    Object for moving a seqset_ctx around on a query and tracking what you're doing
    """

    def __init__(self, my_bg, max_error=5, max_results=5, max_ops=20000):
        """
        Object that will perform all the operations
        """
        self.my_bg = my_bg
        self.max_error = max_error
        self.max_results = max_results
        self.max_ops = max_ops
        self.__cur_ops = max_ops
        self.__bases = ["A", "T", "C", "G"]

    def query(self, sequence, mask=None):
        """
        Queries the biograph for the sequence given the parameters provided during the __init__
        mask is some iterable of length == sequence where elements are True if you can't edit

        Returns a list of Alignments
        """
        self.__cur_ops = self.max_ops
        ctx = ContextState(self.my_bg.find(""), sequence, mask)
        self.__rec_query(ctx)
        return ctx.get_results()

    def __rec_query(self, ctx, cur_err=0, cur_depth=0, pindel=0):
        """
        returns path_status which says if the path you just tried found a value
        if you did find a value and your context is 1 - no use in trying anymore,
        you found an answer

        pindel - was the previous base an insertion/deletion or not (1/-1 or 0)
        """
        # print pindel, ctx.context
        self.__cur_ops -= 1

        if self.__cur_ops <= 0:
            logging.debug('maxops')
            return False

        if len(ctx.get_results()) >= self.max_results:
            logging.debug('results')
            return False

        # quit on this?
        if not ctx.context.valid:
            logging.debug('invalid')
            return False
        print self.__cur_ops, ctx.context, len(ctx.context.sequence), len(ctx.query), cur_err, len(ctx.get_results())
        # logging.debug('depth %d', cur_depth)
        # logging.debug(ctx.alignment.query)
        # logging.debug('aln\n%s', ctx.alignment)
        # logging.debug(str(ctx.context.sequence[::-1]))
        # logging.debug(ctx.alignment.get_sequence() == str(ctx.context.sequence[::-1]))
        # only matches can save you now
        if cur_err >= self.max_error:
            logging.debug('exhaust')
            # exhaust the sequence
            is_valid = ctx.context.valid
            is_aln = ctx.is_aligned()
            while is_valid and not is_aln:
                is_valid = ctx.add_match()
                is_aln = ctx.is_aligned()

            ret = ctx.is_aligned()
            if is_valid and is_aln:
                ctx.record_result()
                ret = True
            return ret

        # we're finished
        if ctx.is_aligned():
            logging.debug('finished')
            ctx.record_result()
            return True

        # Get back to this position
        ctx.add_match()
        path_status = self.__rec_query(ctx, cur_err, cur_depth + 1, pindel=0)
        ctx.roll_back(cur_depth)

        # we have to have used this masked base and that's our only option
        # if ctx.next_mask():
            # ctx.roll_back(cur_depth)
            # return False
        # there's only one path forward and we tried it
        if path_status and ctx.context.end - ctx.context.start == 1:
            logging.debug('short_mat')
            return False

        # mis
        next_base = ctx.next_base()
        any_match = False
        logging.debug('trying mis')
        for i in [x for x in self.__bases if x != next_base]:
            # if pindel=-1 and next_base = ctx.prev
            ctx.add_match(i)
            path_status = self.__rec_query(ctx, cur_err + 1, cur_depth + 1, pindel=0)
            any_match = any_match or path_status
            ctx.roll_back(cur_depth)

        if any_match and ctx.context.end - ctx.context.start == 1:
            logging.debug('short_mis')
            return False

        # del in query so we'll look for an extra base in sequence
        if pindel != 1:  # feel like if there were any-match..nope.. check valids shortcut..
            next_base = ctx.next_base()
            logging.debug('trying_del')
            any_match = False
            # homopolymer
            # can't insert something that should be a match
            for i in [x for x in self.__bases if x != next_base]:
                ctx.add_del(i)
                path_status = self.__rec_query(ctx, cur_err + 1, cur_depth + 1, pindel=-1)
                any_match = any_match or path_status
                ctx.roll_back(cur_depth)
            if any_match and ctx.context.end - ctx.context.start == 1:
                logging.debug('short_del')
                return False

        # ins in query so we'll skip to the next round for sequence
        if pindel != -1:
            logging.debug('trying_ins')
            ctx.add_ins()
            path_status = self.__rec_query(ctx, cur_err + 1, cur_depth + 1, pindel=1)
            ctx.roll_back(cur_depth)
            if path_status and ctx.context.end - ctx.context.start == 1:
                logging.debug('short_ins')
                return False

        # no need to roll back anymore... whatever I return to will do that work
        # And tell the previous guy that this all worked out just fine
        logging.debug('end_of_loop')
        return True


import unittest


class GAlnTestCases(unittest.TestCase):

    def test_simplealn(self):
        """
        aln = Alignment("ATCGAATCGA")
        ans = "ATCGA-ATCGA\n||*|| || ||\nATGGAAAT-GA"
        aln.add_match()
        aln.add_match()
        aln.add_match("G")
        aln.add_match()
        aln.add_match()
        aln.add_del("A")
        aln.add_match()
        aln.add_match()
        aln.add_ins()
        aln.add_match()
        aln.add_match()
        self.assertEqual(str(aln), ans)
        """
        pass

    def test_query1(self):
        bg_file = "/home/english/ajtrio/HG002-NA24385-50x.bg/"
        max_diff = 3
        find1 = "ACATTTTAAACAATTTATACTAGCCCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTTGAATACTTAAGTCAACTGTATGA"
        find2 = "ACATTTTAAACAATTTATACTAGCCCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTGGAATACTTAAGTCAACTGTATGA"
        qseq = "ACATTTTAAACAATTTATACTAGCCCCCCCCGCCCTTTTTTAACAGCTCTCATAAATATGTAGTTCATTTCTCTTGAATACTTAAGTCAACTGTATGA"
        # reversing
        my_bg = BioGraph(bg_file)
        qseq = qseq[::-1]
        alner = GAligner(my_bg)
        ret = alner.query(qseq)
        print '>>>', ret
        for x in ret:
            # x.query_cmp = x.query_cmp[::-1]
            # x.cmp_seq = x.cmp_seq[::-1]
            # x.sequence_cmp = x.sequence_cmp[::-1]
            print str(x)

if __name__ == '__main__':
    unittest.main()


def off_test():
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
