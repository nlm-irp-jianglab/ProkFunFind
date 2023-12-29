from typing import IO, List, Union
import os

from Bio.SearchIO._model.query import QueryResult

from ProkFunFind.detect.bakta_search.bakta_parse import \
    bakta_tab_parse


def pipeline(config: dict,
             in_file: Union[str, IO],
             basedir: Union[str, IO], OrthScore_dict: dict, q_list: dict,
             filter_dict: dict
             ) -> List[QueryResult]:
    """Run Bakta annotation analysis

       Arguments:
           config: A configuration dictionary
           in_file: Input bakta tab file
           basedir: path to directory containing config files
           OrthScore_dict: parsed ortholog table dictionary
           q_list: list of QueryResult objects or empty list

       Returns:
           q_list: An updated list of QueryResults
    """
    # 1. Read baktascan tsv file and parse results
    basedir = os.path.abspath(basedir)+"/"
    qresults = bakta_tab_parse(in_file)

    # 2. Process all QueryResult
    results = []
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
            setattr(i, "geneID", max_dict['geneID'])
            setattr(i, "geneID_weight", max_dict['precision'])
            setattr(i, "detect_tool", "bakta")
            results.append(i)

    # 3. Append hits to overall q_list
    for i in results:
        if len(i) > 0:
            q_list.append(i)
    return q_list
