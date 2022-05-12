from typing import IO, List, Union
from configparser import ConfigParser

from Bio import SearchIO
from Bio.SearchIO._model.query import QueryResult

from GutFunFind.toolkit.base import *
from GutFunFind.detect.emap_search.emap_filter import *


def pipeline(config_file: Union[str, IO],
             in_file: Union[str, IO]
             ) -> List[QueryResult]:
    """Run emappperscan analysis"""
    # 1. Read the configuration file into configuration object
    cf = read_config(config_file)["emapper"]
    basedir = os.path.dirname(os.path.abspath(config_file))+"/"

    # 2. Read emappperscan tsv file and parse results
    qresults = emappper_tab_parse(in_file)

    # 3. Read orthoID info into dictionary
    ortho_file = check_path_existence(basedir + cf["map.ortho_pair"])
    OrthScore_dict = read2orthoDict(ortho_pair_file=ortho_file)

    # 4. Process all QueryResult
    q_list = []
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
            setattr(i, "orthoID", max_dict["orthoID"])
            setattr(i, "orthoID_weight", max_dict["precision"])
            q_list.append(i)

    # filter results based on evalue and thresholds
    if cf.get("filter.config") and cf["filter.config"]:
        filter_path = check_path_existence(basedir + cf["filter.config"])
        filter_res = [emappper_filter(config_file=filter_path, qres=i) for i in q_list]

    # generate final list
    q_list = [i for i in filter_res if len(i) > 0]
    return q_list
