from typing import Dict, IO, List, Set, OrderedDict, Any, Tuple
from collections import defaultdict
import csv
import operator
from sklearn.cluster import DBSCAN
import numpy as np
import functools
import pickle
from Bio import SearchIO
import sys
from read import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import subprocess

basename = sys.argv[1]

gff_file = "/home/jiangx6/data/database/UHGG/MGYG-HGUT/gff/"+basename+".gff"
in_seq_file = "/home/jiangx6/data/database/UHGG/MGYG-HGUT/fna/"+basename+".fna"
in_xml_file = "/home/jiangx6/data/database/UHGG/MGYG-HGUT/ipr_annot/" + \
    basename+"/" + basename+".xml"

g1 = GetGenomeFromGFF(gff3file=gff_file, fnafile=in_seq_file)
[ct.sort() for ct in g1.contigs.values()]

for qresult in SearchIO.parse(in_xml_file, 'interproscan-xml'):
    setattr(g1.genes[qresult.id], qresult.program, qresult)

handle = open("output/"+basename+".pkl", "wb")
pickle.dump(g1, file=handle)
handle.close()


OrthScore_dict = defaultdict(dict)


with open('/home/jiangx6/data/test/GutFunctionalProfile/Code/data/Flagellar.input.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        OrthScore_dict[row[0]][row[1]] = float(row[2])


def ScoreKO(Query, annotdict) -> float:
    try:
        domScore = [annotdict[i]
                    for i in Query.hit_keys if i in annotdict.keys()]
        index, value = max(enumerate(domScore), key=operator.itemgetter(1))
        # return (1 - functools.reduce(lambda a,b : (1-a)*(1-b), [ annotdict[i] for i in Query.hit_keys if i in annotdict.keys()]))
        return (value)
    except:
        return(0)


def AnnotateGene(gene: Gene, KOAnnotdict: OrderedDict[str, Any], exe: str = "InterProScan") -> None:
    if hasattr(gene, exe):
        Query = getattr(gene, exe)
        ScoreList = [ScoreKO(Query, KOAnnotdict[KO])
                     for KO in KOAnnotdict.keys()]
        index, value = max(enumerate(ScoreList), key=operator.itemgetter(1))
        annot = (list(KOAnnotdict.keys())[
                 index], value) if value > 0 else (None, 0)
        #gene.annot = annot
        return(annot)
    else:
        annot = (None, 0)
        #gene.annot = annot
        return(annot)


def CheckGeneQualifiers(ct: Contig, scoredic: str) -> Tuple[List[int], List[Gene], List[str], List[float]]:
    poslist, genelist, termlist, weightlist = [], [], [], []
    for i, gene in enumerate(ct.genes.values()):
        term, weight = AnnotateGene(gene, KOAnnotdict=scoredic)
        if weight > 0:
            poslist.append(i)
            genelist.append(gene)
            termlist.append(term)
            weightlist.append(weight)
    return(poslist, genelist, termlist, weightlist)


def DBSCANCluster(pos: List[int], weight: List[float] = None, eps: float = 2, min_samples: float = 0.9) -> List[int]:
    pos = np.array(pos).reshape(-1, 1)
    cluster = DBSCAN(eps=eps, min_samples=min_samples).fit(
        pos, sample_weight=weight)
    return(list(cluster.labels_))


for ct in g1.contigs.values():
    poslist, genelist, termlist, weightlist = CheckGeneQualifiers(
        ct=ct, scoredic=OrthScore_dict)

    if len(poslist) > 1:
        c = DBSCANCluster(pos=poslist, weight=weightlist)
        for i in range(0, len(c)):
            if c[i] >= 0:
                print(ct.id, genelist[i].id, c[i], termlist[i], weightlist[i])
