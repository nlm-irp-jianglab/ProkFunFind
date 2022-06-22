from typing import IO, List, Union
from configparser import ConfigParser

from Bio import SearchIO
from Bio.SearchIO._model.query import QueryResult

from GutFunFind.toolkit.utility import *
from GutFunFind.detect.kofam_search.kofam_filter import *


def pipeline(config: dict,
             in_file: Union[str, IO],
             basedir, OrthScore_dict: dict, q_list: dict,
             ) -> List[QueryResult]:
    """Run kofamscan analysis"""
    # 1. Read kofamscan tsv file and parse results
    qresults = kofam_tab_parse(in_file)

    # 2. Process all QueryResult
    tmp_list = []
    for qres in qresults:
        # remove and query results that do not hit to searched KOs
        i = qres.hsp_filter(lambda hsp: hsp.hit_id in OrthScore_dict.keys())

        def sort_key(hit):
            return hit.hsps[0].evalue

        if len(i) > 0:
            # sort KO hits by evalue.
            i.sort(key=sort_key, reverse=False, in_place=True)
            max_dict = OrthScore_dict[i.hits[0].id]

            # set the QueryResult attributes
            setattr(i, "orthoID", max_dict['orthoID'])
            setattr(i, "orthoID_weight", max_dict['precision'])
            setattr(i, "detect_tool", "kofamscan")
            tmp_list.append(i)

    # 3. filter results based on evalue and thresholds
    if config['filter']:
        filter_res = [kofam_filter(config=config, qres=i, basedir=basedir) for i in tmp_list]
    else:
        filter_res = tmp_list

    # 4. Append significant hits to q_list
    for i in filter_res:
        if len(i) > 0:
            q_list.append(i)
    return q_list
