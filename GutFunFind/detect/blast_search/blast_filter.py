import csv
import os
import operator
from collections import defaultdict

from typing import IO, Union
from Bio.SearchIO._model.query import QueryResult
from GutFunFind.toolkit.base import read_config, check_path_existence, read2orthoDict

# def blast_filter(config: ConfigParser, qres: QueryResult) -> QueryResult:
def blast_filter(config_file: Union[str, IO], qres: QueryResult) -> QueryResult:

    #cf = config
    cf = read_config(config_file)
    basedir = os.path.dirname(os.path.abspath(config_file))+"/"

    global_evalue = cf["filter.global"]["evalue"]
    global_ident = cf["filter.global"]["ident_pct"]

    ##################################################################
    #  User can customize the filter function to remove QueryResult  #
    ##################################################################
    ops = {
        "<=": operator.le,
        ">=": operator.ge,
        ">": operator.gt,
        "<": operator.lt,
        "==": operator.eq,
        "!=": operator.ne
    }

    filter_dict = defaultdict(list)

    # if there is section filter.local and that section has filter_file, the info of the filter_file will be passed to  filter_dict; Otherwise filter_dict remains empty
    if cf.has_option("filter.local", "filter_file"): 
        hit_filter_file = check_path_existence(basedir + cf["filter.local"]["filter_file"])
        # check if file exist or empty
        with open(hit_filter_file) as filter_file:
            for row in csv.reader(filter_file, delimiter='\t'):
                filter_dict[row[0]].append(
                    {"attr": row[1], "cpfun": ops[row[2]], "value": float(row[3])})

    def hsp_filter_func(hsp):
        status = True
        if hsp.hit_id in filter_dict:
            for one in filter_dict[hsp.hit_id]:
                if one["cpfun"](getattr(hsp, one['attr']), one["value"]):
                    pass
                else:
                    status = False
                    break
        else:
            if hsp.evalue > float(global_evalue) or hsp.ident_pct < float(global_ident):
                status = False
        return status

    return qres.hsp_filter(hsp_filter_func)


def blast_ortho(qres: QueryResult, ortho_pair_file: str) -> QueryResult:

    OrthScore_dict = read2orthoDict(ortho_pair_file=ortho_pair_file)

    ######################################
    #  Sort method can be defined later  #
    ######################################
    # sort by the QueryResult by hit length
    def sort_key(hit):
        return sum([hsp.aln_span for hsp in hit.hsps])

    qres.sort(key=sort_key, reverse=True, in_place=False)

    # Use the top matched hit to assgin orthoID to gene
    hit_key = qres.hit_keys[0]
    max_dict = OrthScore_dict[hit_key]

    setattr(qres, "orthoID", max_dict["orthoID"])
    setattr(qres, "orthoID_weight", max_dict["precision"])

    return qres
