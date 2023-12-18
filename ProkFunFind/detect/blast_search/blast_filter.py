import csv
import operator
from collections import defaultdict

from Bio.SearchIO._model.query import QueryResult
from ProkFunFind.toolkit.utility import check_path_existence


def blast_filter(config: dict, qres: QueryResult, basedir: str, filter_dict: dict) -> QueryResult:
    """Function to filter blast query results

       Arguments:
           config: A configuration dictionary
           qres: A Bio.SearchIO._model.query.QueryResult object
           basedir: A path to the directory containing the config files

      Returns:
           Bool. If true then the hit has passed the filtering steps

    """

    global_evalue = config['blast'].get('evalue', 0.01)
    global_ident = config['blast'].get('ident_pct', 30)

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

    def hsp_filter_func(hsp):
        """Helper function to filter high scoring pairs

           Arguments:
              hsp: A QueryResult HSP object

           Returns:
              Bool, true if hit passes filtering.

        """
        status = True

        if hsp.hit_id in filter_dict:
            for one in filter_dict[hsp.hit_id]:
                if one['cpfun'](getattr(hsp, one['attr']), one['value']):
                    pass
                else:
                    status = False
                    break
        else:
            if hsp.evalue > float(global_evalue) or \
                    hsp.ident_pct < float(global_ident):
                status = False
        return status

    return qres.hsp_filter(hsp_filter_func)


def blast_ortho(qres: QueryResult, OrthScore_dict: dict) -> QueryResult:
    """Function to sort and assign properties to blast results

       Arguments:
          qres: Bio.SearchIO._model.query.QueryResult object
          OrthScore_dict:

       Returns:
          qres: Updated QueryResult object

    """
    ######################################
    #  Sort method can be defined later  #
    ######################################
    # sort by the QueryResult by hit length
    def sort_key(hit):
        return sum([hsp.aln_span for hsp in hit.hsps])

    qres.sort(key=sort_key, reverse=True, in_place=False)

    # Use the top matched hit to assgin geneID to gene
    hit_key = qres.hit_keys[0]
    max_dict = OrthScore_dict[hit_key]

    setattr(qres, "geneID", max_dict['geneID'])
    setattr(qres, "geneID_weight", max_dict['precision'])
    setattr(qres, "detect_tool", "blast")

    return qres
