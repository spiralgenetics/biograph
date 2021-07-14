"""
Toy script to take a bunch of force-called results and unite them based on the
find_variants results.
NOT a production script
"""
import sys
import json
import random
from collections import Counter, defaultdict

import progressbar
from biograph import libspiralapi as spapi

import biograph as bgsdk

familyId = 14497
individuals = [("SSC11549", "p1"),
               ("SSC11552", "fa"),
               ("SSC11550", "mo"),
               ("SSC11553", "s1")]  # ignore for now


def uniteCalls():
    """
    1) Load in and unite based on the original coordinates.
    I'm taking out MIS and making their columns consistent
    """
    # Put them in my two piles
    # par - parliament calls
    # bgfv - findvariant calls
    # bggt - queryKmer over bgfv coordinates
    allCalls = {"par": defaultdict(list),
                "bgfv": defaultdict(list),
                "bggt": defaultdict(list)}
    for called, idc in individuals:
        for forced, idf in individuals:
            with open("%scalls_%sdata_bgForce.txt" % (called, forced), 'r') as fh:
                for line in fh:
                    data = line.strip().split('\t')
                    key = "%s:%s-%s" % (data[0], data[1], data[2])
                    if data[4] == 'MIS':
                        allCalls["par"][key].append(data)
                    elif data[9] == 'skip':
                        allCalls["par"][key].append(data)
                    elif data[13] == 'ref':  # no new coords to work off of
                        allCalls["par"][key].append(data)
                    else:
                        key = "%s:%s-%s" % (data[10], data[11], data[12])
                        allCalls["bgfv"][key].append(data)
    sys.stderr.write("%d parliment pile calls\n" % (len(allCalls["par"])))
    sys.stderr.write("%d biograph pile calls\n" % (len(allCalls["bgfv"])))

    todo = []
    # Downsampling
    for i in allCalls["bgfv"]:
        if allCalls["bgfv"][i][0][4] == "DEL":
            todo.append(i)

    # todo = random.sample(todo, 10)
    callCount = len(todo) * len(individuals)
    bar = progressbar.ProgressBar(redirect_stdout=True, max_value=callCount, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(callCount), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])

    n = 0
    for individual, iid in individuals:
        my_bg = bgsdk.BioGraph(individual + ".hgv.seqset", individual + ".hgv.readmap")
        my_sg = spapi.SpiralGenome(my_bg, "/share/reference/human_grch38_nospace/")
        for callKey in todo:
            n += 1
            bar.update(n)
            data = allCalls["bgfv"][callKey][0]
            # Only deletions for now
            if data[4] != "DEL":
                continue

            chrom, start, end = data[10], int(data[11]), int(data[12])
            try:
                x = my_sg.queryKmer(chrom, start, end, 50)
            except:
                continue  # don't care
            allCalls["bggt"][callKey].append((iid, x[0][0], x[0][1], len(x[1]), len(x[2])))

    return allCalls["bggt"]
    # with open("results.json", 'w') as fout:
        # fout.write(json.dumps(allCalls["bggt"], indent=4))


def uniteCallsINS():
    """
    1) Load in and unite based on the original coordinates.
    I'm taking out MIS and making their columns consistent
    """
    # Put them in my two piles
    # par - parliament calls
    # bgfv - findvariant calls
    # bggt - queryKmer over bgfv coordinates
    allCalls = {"par": defaultdict(list),
                "bgfv": defaultdict(list),
                "bggt": defaultdict(list)}
    for called, idc in individuals:
        for forced, idf in individuals:
            with open("%scalls_%sdata_bgForce.txt" % (called, forced), 'r') as fh:
                for line in fh:
                    data = line.strip().split('\t')
                    key = "%s:%s-%s" % (data[0], data[1], data[2])
                    if data[4] == 'MIS':
                        allCalls["par"][key].append(data)
                    elif data[9] in ['skip', 'mix']:
                        allCalls["par"][key].append(data)
                    elif data[13] == 'ref':  # no new coords to work off of
                        allCalls["par"][key].append(data)
                    else:
                        key = "%s:%s-%s" % (data[10], data[11], data[12])
                        allCalls["bgfv"][key].append(data)
    sys.stderr.write("%d parliment pile calls\n" % (len(allCalls["par"])))
    sys.stderr.write("%d biograph pile calls\n" % (len(allCalls["bgfv"])))

    todo = []
    # Downsampling
    for i in allCalls["bgfv"]:
        if allCalls["bgfv"][i][0][4] == "INS":
            todo.append(i)

    # todo = random.sample(todo, 10)
    callCount = len(todo) * len(individuals)
    bar = progressbar.ProgressBar(redirect_stdout=True, max_value=callCount, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(callCount), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])

    n = 0
    kmer = 50
    for individual, iid in individuals:
        my_bg = bgsdk.BioGraph(individual + ".hgv.seqset", individual + ".hgv.readmap")
        my_sg = spapi.SpiralGenome(my_bg, "/share/reference/human_grch38_nospace/")
        for callKey in todo:
            n += 1
            bar.update(n)
            data = allCalls["bgfv"][callKey][0]
            # Only deletions for now
            if data[4] != "INS":
                continue
            chrom, start, end = data[10], int(data[11]), int(data[12])
            # try:
            if True:
                altSeq1 = str(my_sg.reference.make_range(
                    chrom, start - kmer / 2, start + 1, True).sequence) + data[-1][:kmer / 2]
                altSeq2 = data[-1][-(kmer / 2):] + str(
                    my_sg.reference.make_range(chrom, end, end + kmer / 2, True).sequence)
                refSeq1 = str(my_sg.reference.make_range(chrom, start - kmer / 2, start + kmer / 2, True).sequence)
                # refSeq2 = str(my_sg.reference.make_range(chrom, end - kmer/2, end + kmer/2, True).sequence)
                altCov = set()
                refCov = set()
                altCov.update(my_sg.my_bg.find_reads(altSeq1, 100))
                altCov.update(my_sg.my_bg.find_reads(altSeq2, 100))
                refCov.update(my_sg.my_bg.find_reads(refSeq1, 100))
                # refCov.update(my_sg.my_bg.find_reads(altSeq1, 100))
                x = spapi.genotyper(len(refCov) + len(altCov), len(altCov)), refCov, altCov
            # except:
                # continue#don't care
            allCalls["bggt"][callKey].append([iid, x[0][0], x[0][1], len(x[1]), len(x[2]), data[-1]])

    with open("results_ins.json", 'w') as fout:
        fout.write(json.dumps(allCalls["bggt"], indent=4))


def resultsToTable(data):
    # with open("results_ins.json", 'r') as fh:
    # with open("results.json", 'r') as fh:
        # data = json.load(fh)

    def regenotype(data):
        THRESHOLD = 1
        data[0] = str(data[0])
        data[1] = str(data[1])
        refcnt = data[3]
        altcnt = data[4]
        if refcnt >= THRESHOLD:
            if altcnt >= THRESHOLD:
                data[1] = "0/1"
            else:
                data[1] = "0/0"
        else:
            if altcnt >= THRESHOLD:
                data[1] = "1/1"
            else:
                data[1] = "./."
        del(data[2])
        # del(data[-1])

    numCalls = 0

    genotypes = defaultdict(Counter)
    mendelGood = 0
    mendelBad = 0
    coverageProblem = 0
    noCalls = 0
    low_coverage = 0
    callSizes = []

    for key in data:
        numCalls += 1
        p, f, m, s = data[key]
        coords = key.replace(':', '\t').replace('-', '\t').split('\t')

        chrom, start, end = coords[0], int(coords[1]), int(coords[2])
        # seq = p[-1]
        coords.append(str(end - start))
        # coords.append(len(seq)) ins
        regenotype(p)
        regenotype(f)
        regenotype(m)
        regenotype(s)
        genotypes[p[0]][p[1]] += 1
        genotypes[f[0]][f[1]] += 1
        genotypes[m[0]][m[1]] += 1
        if p[1] in ['./.', '0/0'] and f[1] in ['./.', '0/0'] and m[1] in ['./.', '0/0']:
            noCalls += 1
            coords.append("low_quality")
        elif p[1] == '0/0':
            if f[1].count('0') and m[1].count('0'):
                coords.append("consistent")
                mendelGood += 1
            else:
                coords.append("inconsistent")
                mendelBad += 1
        elif p[1] == '0/1':
            if (f[1] == '1/1' and m[1] == '1/1'):
                coords.append("consistent")
                mendelBad += 1
            elif (f[1] == '0/0' and m[1] == '0/0'):
                mendelBad += 1
                if p[3] >= 5:
                    coords.append("denovo")
                else:
                    coords.append("low_coverage")
            elif f[1].count('1') + m[1].count('1'):
                coords.append("consistent")
                mendelGood += 1
            else:
                coords.append("low_coverage")
                low_coverage += 1
        elif p[1] == '1/1':
            if (f[1] == '0/0' or m[1] == '0/0'):
                coords.append("inconsistent")
                mendelBad += 1
            elif (f[1] == './.' and m[1] == './.'):
                coords.append("low_coverage")
                low_coverage += 1
            else:
                coords.append("consistent")
                mendelGood += 1
        elif p[1] == './.':
            coords.append("low_coverage")
            coverageProblem += 1

        print("\t".join([str(x) for x in coords + data[key]]).replace('[', '')
              .replace(']', '').replace("'", '').replace(', ', ':'))
    sys.stsderr.write("numcalls %d\n" % numCalls)
    """"
    print("numcalls %d" % numCalls)
    print("noCalls %d" % noCalls)
    print("genotypes %s" % (json.dumps(genotypes)))
    print("mendelGood %d" % mendelGood)
    print("mendelBad %d" % mendelBad)
    print("coverageProblem %d" % coverageProblem)
    print("low_coverage %d" % low_coverage)
    """

if __name__ == '__main__':
    calls = uniteCalls()
    resultsToTable(calls)
    # uniteCallsINS()
