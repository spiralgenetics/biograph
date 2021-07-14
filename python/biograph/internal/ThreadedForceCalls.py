import sys
import copy
import signal
import traceback
from collections import Counter, OrderedDict
import multiprocessing
import random
import progressbar
import logging

import biograph as bgsdk


class Alarm(Exception):
    pass


def alarm_handler(signum, frame):
    raise Alarm


class Consumer(multiprocessing.Process):

    """
    Basic Consumer. Follow the two queues with your *args and **kwargs that should be sent
    to the task when __call__ 'd

    NOTE! args can't hold anything that isn't pickle-able for the subprocess
    task_queue, result_que, biograph_file_path?
    """

    def __init__(self, task_queue, result_queue, force_indiv):
        multiprocessing.Process.__init__(self)
        my_sg.queryRegion("chr1", 400000, 405000)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        try:
            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    # Poison pill means shutdown
                    self.task_queue.task_done()
                    break
                try:
                    next_task()
                except Exception as e:
                    logging.error("Exception raised in task - %s" % (str(e)))

                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    logging.error("Dumping Traceback:")
                    traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stderr)

                    next_task.failed = True
                    next_task.errMessage = str(e)

                self.result_queue.put(next_task)
                self.task_queue.task_done()

            return
        except Exception as e:
            logging.error("Consumer %s Died\nERROR: %s" % (self.name, e))
            return


class ForceCall():

    """
    Counts all the errors in the region
    """

    def __init__(self, svcall):
        self.svcall = svcall
        self.failed = False
        self.errMessage = ""

    def __call__(self):
        """
        Takes a pysam.Samfile
        """
        # logging.info("Starting %s" % (self.name))
        # regName =  "%s:%d-%d" % (self.chrom, self.start, self.end)
        my_new = copy.copy(newFields)
        chrom, start, end, svid, svtype, size, pforce, gt, gap = self.svcall[:9]

        # Ignore big spans
        if abs(start - end) > 5000 or gap != "." \
           or svtype == "MIS" or start >= end:
            my_new["BG_force"] = 'skip'
            self.svcall.append(toTab(my_new))
            return

        ASM_BUFFER = 300
        bufStart = max(1, start - ASM_BUFFER)
        bufEnd = end + ASM_BUFFER
        # I'm only doing subset for testing
        # if random.random() > .05:
            # finisher(input_queue, output_queue, svcall, my_new)
            # continue
        # if chrom != "chr1":
            # finisher(input_queue, output_queue, svcall, my_new)
            # continue

        try:
            # time it out because we don't want to clog the pool
            signal.signal(signal.SIGALRM, alarm_handler)
            signal.alarm(90)  # 90 seconds MAX
            variants = my_sg.queryRegion(chrom, bufStart, bufEnd)
            signal.alarm(0)  # reset the alarm
        except:  # Can't catch custom, yet, should probably get the alarm
            # need to watch out for this
            my_new["BG_force"] = 'error'
            self.svcall.append(toTab(my_new))
            return

        # Then figure out which one of these variants to keep and report on
        keeps = []
        for v in variants:
            if v.anno == svtype and v.length >= 50:
                keeps.append(v)

        if len(keeps) == 0:
            if len(variants) == 1:
                my_new["BG_force"] = "asm"
                my_new["BG_svtype"] = "ref"
                # my_new["BG_mu_refCov"] = variants[0].mRef # Why not include this?
                my_new["BG_gt"] = "0/0"
                self.svcall.append(toTab(my_new))
            else:  # Got some noise I don't want to deal with, yet
                my_new["BG_force"] = "mix"
                self.svcall.append(toTab(my_new))
        else:
            # just the biggest
            keeps.sort(cmp=lambda x, y: x.length - y.length, reverse=True)
            my_new["BG_force"] = 'asm'
            my_new["BG_chrom"] = keeps[0].chromosome
            my_new["BG_start"] = keeps[0].refStart
            my_new["BG_end"] = keeps[0].refEnd
            my_new["BG_svtype"] = keeps[0].anno
            my_new["BG_length"] = keeps[0].length
            my_new["BG_gt"] = keeps[0].genotype
            my_new["BG_mu_refCov"] = keeps[0].mRef
            my_new["BG_mu_altCov"] = keeps[0].mAlt
            my_new["BG_seq"] = keeps[0].altSeq if keeps[0].anno == "INS" else "."
            self.svcall.append(toTab(my_new))
        return


"""
familyId = 14509

individuals = [("SSC12655", "fa"),
               ("SSC12647", "mo"),
               ("SSC12643", "p1"),
               ("SSC12656", "s1")]
"""
familyId = 14497
individuals = [("SSC11552", "fa"),
               ("SSC11550", "mo"),
               ("SSC11549", "p1"),
               ("SSC11553", "s1")]

newFields = OrderedDict()
# Notes on BG_force
# skip = skipped it; fail = internal error; asm = graph reference path found*; kmer = KmerQueried"; both = asm & kmer
# mix = buncha snps around
# I'll manually put in the header
newFields["BG_force"] = '.'
newFields["BG_chrom"] = '.'
newFields["BG_start"] = '.'
newFields["BG_end"] = '.'
newFields["BG_svtype"] = '.'
newFields["BG_gt"] = '.'
newFields["BG_length"] = '.'
newFields["BG_mu_refCov"] = '.'
newFields["BG_mu_altCov"] = '.'
newFields["BG_seq"] = "."
# newFields["BG_kmer_ref"] = '.'
# newFields["BG_kmer_alt"] = '.' #Only done on deletions for now


def toTab(data):
    """newFields helper method"""
    return "\t".join([str(data[x]) for x in data])


def process_svcall(thread_id, force_indiv, input_queue, output_queue):
    my_bg = bgsdk.BioGraph(force_indiv + ".hgv.seqset", force_indiv + ".hgv.readmap")
    my_sg = bgsdk.SpiralGenome(my_bg, "/share/reference/human_grch38_nospace/")

    def finisher(input_queue, output_queue, svcall, my_new):
        output_queue.put(svcall)
        try:
            input_queue.task_done()
        except ValueError:
            pass
        updateprogress()
    while True:
        svcall = input_queue.get()
        if svcall is None:
            input_queue.task_done()
            return
        my_new = copy.copy(newFields)
        chrom, start, end, svid, svtype, size, pforce, gt, gap = svcall[:9]

        # Ignore big spans
        if abs(start - end) > 5000 or gap != "." \
           or svtype == "MIS" or start >= end:
            my_new["BG_force"] = 'skip'
            finisher(input_queue, output_queue, svcall, my_new)
            continue

        ASM_BUFFER = 300
        bufStart = max(1, start - ASM_BUFFER)
        bufEnd = end + ASM_BUFFER
        # I'm only doing subset for testing
        # if random.random() > .05:
            # finisher(input_queue, output_queue, svcall, my_new)
            # continue
        # if chrom != "chr1":
            # finisher(input_queue, output_queue, svcall, my_new)
            # continue

        try:
            # time it out because we don't want to clog the pool
            signal.signal(signal.SIGALRM, alarm_handler)
            signal.alarm(90)  # 90 seconds MAX
            variants = my_sg.queryRegion(chrom, bufStart, bufEnd)
            signal.alarm(0)  # reset the alarm
        except:  # Can't catch custom, yet, should probably get the alarm
            # need to watch out for this
            my_new["BG_force"] = 'error'
            finisher(input_queue, output_queue, svcall, my_new)
            print('error on ', svcall)
            # exit(1)
            continue

        # Then figure out which one of these variants to keep and report on
        keeps = []
        for v in variants:
            if v.anno == svtype and v.length >= 50:
                keeps.append(v)
        if len(keeps) == 0:
            if len(variants) == 1:
                my_new["BG_force"] = "asm"
                my_new["BG_svtype"] = "ref"
                my_new["BG_gt"] = "0/0"
                finisher(input_queue, output_queue, svcall, my_new)
                continue
            else:
                # Got some noise I don't want to deal with, yet
                my_new["BG_force"] = "mix"
                finisher(input_queue, output_queue, svcall, my_new)
                continue
        # just the biggest
        keeps.sort(cmp=lambda x, y: x.length - y.length, reverse=True)
        # keeps = keeps[0]
        my_new["BG_force"] = 'asm'
        my_new["BG_chrom"] = keeps[0].chromosome
        my_new["BG_start"] = keeps[0].refStart
        my_new["BG_end"] = keeps[0].refEnd
        my_new["BG_svtype"] = keeps[0].anno
        my_new["BG_length"] = keeps[0].length
        my_new["BG_gt"] = keeps[0].genotype
        my_new["BG_mu_refCov"] = keeps[0].mRef
        my_new["BG_mu_altCov"] = keeps[0].mAlt
        my_new["BG_seq"] = keeps[0].altSeq if keeps[0].anno == "INS" else "."
        finisher(input_queue, output_queue, svcall, my_new)

        # Then let's do the simple kmer querying?

print familyId
MAX_THREAD_COUNT = 16

count_lock = multiprocessing.Lock()
count_value = multiprocessing.Value('i', 0)


def updateprogress():
    with count_lock:
        count_value.value += 1


def progress_bar_watcher(callCount):
    # You kill this when you're done with it
    bar = progressbar.ProgressBar(redirect_stdout=True, max_value=callCount, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(callCount), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])
    cont = True
    while cont:
        with count_lock:
            if count_value.value < bar.max_value - 1:
                bar.update(count_value.value)
            else:
                cont = False

cur_indiv = sys.argv[1]
force_indiv = sys.argv[2]

# reset the queues! I'm new(s)!
input_queue = multiprocessing.JoinableQueue()
output_queue = multiprocessing.Queue()
callCount = 0

sys.stderr.write("Loading BioGraph\n")
my_bg = bgsdk.BioGraph(force_indiv + ".hgv.seqset", force_indiv + ".hgv.readmap")
my_sg = bgsdk.SpiralGenome(my_bg, "/share/reference/human_grch38_nospace/")

sys.stderr.write("Forcing %s over %s calls\n" % (force_indiv, cur_indiv))
jobs = [Consumer(input_queue, output_queue, force_indiv) for i in range(MAX_THREAD_COUNT)]
for w in jobs:
    w.start()

with open(cur_indiv + '.hgv.parliament.calls.txt', 'r') as fh:
    fh.readline()
    for line in fh:
        data = line.strip().split('\t')
        data[1] = int(data[1])
        data[2] = int(data[2])
        data[3] = cur_indiv + "_" + data[3]
        data[5] = int(data[5])
        input_queue.put(ForceCall(data))
        callCount += 1


bar = progressbar.ProgressBar(redirect_stdout=True, max_value=callCount, widgets=[
    ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(callCount), '] ',
    progressbar.Bar(),
    ' (', progressbar.ETA(), ') ',
])

found = 0
with open(cur_indiv + 'calls_' + force_indiv + "data" + "_bgForce.txt", 'w') as fout:
    while callCount:
        data = output_queue.get()
        fout.write("\t".join([str(y) for y in data.svcall]) + '\n')
        callCount -= 1
        found += 1
        bar.update(found)
bar.finish()

for j in jobs:
    callCount += 1
    input_queue.put(None)
