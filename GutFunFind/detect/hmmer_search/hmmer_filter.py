import csv
import os
import operator
from collections import defaultdict

from typing import IO, Union
from Bio.SearchIO._model.query import QueryResult
from GutFunFind.toolkit.utility import read_config, check_path_existence, read2orthoDict

def hmmer_filter(config: dict, qres: QueryResult, basedir=str) -> QueryResult:
    """Function to filetr HMMER search results

       Arguments
           config: A configuration dictionary
           qres: A Bio.SearchIO._model.query.QueryResult object
           basedir: Path to the configuration files directory

       Returns:
           status: True if query passes filtering.
    """
    global_evalue = config['filter']['evalue'] if  config['filter'].get('evalue') else 10
    global_bitscore = config['filter']['bitscore'] if config['filter'].get('bitscore') else 0

    ##################################################################
    #  User can customize the filter function to remove QueryResult  #
    ##################################################################
    ops = {
        '<=': operator.le,
        '>=': operator.ge,
        '>': operator.gt,
        '<': operator.lt,
        '==': operator.eq,
        '!=': operator.ne
    }

    filter_dict = defaultdict(list)

    query_len = qres.seq_len;

    # if there is section filter.local and that section has filter_file, the info of the filter_file will be passed to  filter_dict; Otherwise filter_dict remains empty
    if config.has_option("filter", "filter_file"):
        hit_filter_file = check_path_existence(basedir + config['filter']['filter_file'])
        with open(hit_filter_file) as filter_file:
            for row in csv.reader(filter_file, delimiter="\t"):
                filter_dict[row[0]].append({'attr': row[1], 'cpfun': ops[row[2]], 'value': float(row[3])})

    def hit_filter_func(hit):
        status = True
        if hit.id in filter_dict:
            hit_len = hit.seq_len
            for one in filter_dict[hit.id]:
                if one['attr'] == "query_match_pct":
                    check1 = [ one['cpfun']((hsp.query_span)/query_len, one['value']) for hsp in hit.hsps ]
                    if not sum(check1) :
                        status = False
                elif one['attr'] == "hit_match_pct":
                   check2 = [ one['cpfun']((hsp.hit_span)/hit_len, one['value']) for hsp in hit.hsps ]
                   if not sum(check2) :
                       status = False
                elif one['cpfun'](getattr(hit, one['attr']), one['value']):
                    pass
                else:
                    status = False
                    break
        else:
            if hit.evalue > float(global_evalue) or hit.bitscore < float(global_bitscore):
                status = False
        return status

    return qres.hit_filter(hit_filter_func)
