import csv
import operator
from collections import defaultdict

from Bio.SearchIO._model.query import QueryResult
from ProkFunFind.toolkit.utility import check_path_existence


def hmmer_filter(config: dict, qres: QueryResult, basedir: str, filter_dict: dict) -> QueryResult:
    """Function to filetr HMMER search results

       Arguments
           config: A configuration dictionary
           qres: A Bio.SearchIO._model.query.QueryResult object
           basedir: Path to the configuration files directory

       Returns:
           status: True if query passes filtering.
    """
    global_evalue = config['hmmer'].get('evalue', 0.01)
    global_bitscore = config['hmmer'].get('bitscore', 0)

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

    query_len = qres.seq_len

    def hit_filter_func(hit):
        status = True
        if hit.id in filter_dict:
            hit_len = hit.seq_len
            for one in filter_dict[hit.id]:
                if one['attr'] == "query_match_pct":
                    check1 = [one['cpfun']((hsp.query_span)/query_len,
                              one['value']) for hsp in hit.hsps]
                    if not sum(check1):
                        status = False
                elif one['attr'] == "hit_match_pct":
                    check2 = [one['cpfun']((hsp.hit_span)/hit_len,
                              one['value']) for hsp in hit.hsps]
                    if not sum(check2):
                        status = False
                elif one['cpfun'](getattr(hit, one['attr']), one['value']):
                    pass
                else:
                    status = False
                    break
        else:
            if hit.evalue > float(global_evalue) or \
                    hit.bitscore < float(global_bitscore):
                status = False
        return status

    return qres.hit_filter(hit_filter_func)
