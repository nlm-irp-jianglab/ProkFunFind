from configparser import ConfigParser
import csv
from typing import Dict, IO, List, Set, OrderedDict, Callable
from Bio.SearchIO._model.query import QueryResult


def blast_filter(config: ConfigParser, qres: QueryResult) -> QueryResult:
    cf = config

    ##################################################################
    #  User can customize the filter function to remove QueryResult  #
    ##################################################################
    global_evalue = cf["filter.global"]["evalue"]
    if global_evalue:
        def hsp_filter_func(hsp): return hsp.evalue < float(global_evalue)
        # hit_filter_func
        return qres.hsp_filter(hsp_filter_func)
    else:
        return qres


def blast_ortho(qres: QueryResult, ortho_pair_file: str) -> QueryResult:

    # Create a dict of dict[bait][orthoID] = specificity
    from collections import defaultdict
    OrthScore_dict = defaultdict(dict)
    #############################################################
    #  User can change the last column to indicate specificity  #
    #############################################################
    with open(ortho_pair_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        header = next(csv_reader)
        col_num = len(header)
        if col_num == 2:
            OrthScore_dict[header[1]][header[0]] = 1
            for row in csv_reader:
                OrthScore_dict[row[1]][row[0]] = 1
        else:
            OrthScore_dict[header[1]][header[0]] = float(row[2])
            for row in csv_reader:
                OrthScore_dict[row[1]][row[0]] = float(row[2])

    ######################################
    #  Sort method can be defined later  #
    ######################################
    # sort by the QueryResult by hit length
    def sort_key(hit): return sum([hsp.aln_span for hsp in hit.hsps])
    qres.sort(key=sort_key, reverse=True, in_place=False)

    # Use the top matched hit to assgin orthoID to gene
    hit_key = qres.hit_keys[0]
    ortho_dict = OrthScore_dict[hit_key]
    setattr(qres, "orthoID", list(ortho_dict.keys())[0])
    setattr(qres, "orthoID_weight", list(ortho_dict.values())[0])

    return qres
