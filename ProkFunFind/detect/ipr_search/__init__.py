from typing import IO, List, Union
from configparser import ConfigParser

from Bio import SearchIO
from Bio.SearchIO._model.query import QueryResult

from ProkFunFind.toolkit.utility import *
from ProkFunFind.detect.ipr_search.interproscan_tab import InterproscanTabParser, ipr_tab_parse


def pipeline(config: dict,
             in_file: Union[str, IO],
             fmt: str,
             basedir: str, OrthScore_dict: dict, q_list: dict) -> List[QueryResult]:
    """Main function for running IPR based search

       Arguments:
           config: A configuration dictionary
           in_file: Path to interproscan output file
           fmt: format of input file (xml or tsv)
           basedir: path to configuration file directory
           OrthScore_dict: OrthScore_dict dictionary
           q_list: list of query results or empty list

       Returns:
           q_list: Updated list of QueryResults
    """
    # 1. Read the input files
    if fmt == "xml":
        qresults = SearchIO.parse(in_file, "interproscan-xml")
    elif fmt == "tsv":
        qresults = ipr_tab_parse(in_file)

    # 2. Process all QueryResult
    for qres in qresults:

        # remove QueryResult that doest not hit any domain in function-related domain list
        i = qres.hsp_filter(lambda hsp: hsp.hit_id in OrthScore_dict.keys())

        # for those without any hit match to the domain
        if len(i) > 0:
            # sort the hits by precision
            i.sort(
                key=lambda hit: OrthScore_dict[hit.id]['precision'], reverse=True)

            max_dict = OrthScore_dict[i.hits[0].id]

            # set the QueryResult attribution
            setattr(i, "orthoID", max_dict['orthoID'])
            setattr(i, "orthoID_weight", max_dict['precision'])
            setattr(i, "detect_tool", "interproscan")
            q_list.append(i)

    return q_list
