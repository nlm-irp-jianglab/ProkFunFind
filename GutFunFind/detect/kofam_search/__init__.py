from typing import IO, List, Union
from configparser import ConfigParser

from Bio import SearchIO
from Bio.SearchIO._model.query import QueryResult

from GutFunFind.toolkit.base import *
from GutFunFind.detect.kofam_search.kofam_filter import *


def pipeline(config: dict,
             in_file: Union[str, IO],
             basedir
             ) -> List[QueryResult]:
    """Run kofamscan analysis"""
    # 1. Read the configuration file into configuration object
    # cf = read_config(config_file)["kofamscan"]
    # basedir = os.path.dirname(os.path.abspath(config_file))+"/"

    # 2. Read kofamscan tsv file and parse results
    qresults = kofam_tab_parse(in_file)

    # 3. Read orthoID info into dictionary
    ortho_file = check_path_existence(basedir + config['kofamscan']["map.ortho_pair"])
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
    if config['filter']:
        # filter_path = check_path_existence(basedir + cf["filter.config"])
        filter_res = [kofam_filter(config=config, qres=i, basedir=basedir) for i in q_list]

    # generate final list
    q_list = [i for i in filter_res if len(i) > 0]
    return q_list
