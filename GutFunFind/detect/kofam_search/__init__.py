from typing import IO, List, Union
from configparser import ConfigParser

from Bio import SearchIO
from Bio.SearchIO._model.query import QueryResult

from GutFunFind.toolkit.base import *
from GutFunFind.detect.kofam_search.kofam_filter import KofamscanTabParser, kofam_tab_parse, kofam_filter


def pipeline(config_file: Union[str, IO],
             in_file: Union[str, IO]
             ) -> List[QueryResult]:
    # 1. Read the configuration file into configuration object
    cf = read_config(config_file)["kofamscan"]
    basedir = os.path.dirname(os.path.abspath(config_file))+"/"

    # 2. Read kofamscan xml file
    qresults = kofam_tab_parse(in_file)

    # 3. Read score and orthoID info into dictionary
    ortho_file = check_path_existence(basedir + cf["orthoID_domain_precision"])
    OrthScore_dict = read2orthoDict(ortho_pair_file=ortho_file)

    # 4. Process all QueryResult
    q_list = []
    for qres in qresults:

        # remove QueryResult that doest not hit any domain in function-related domain list
        i = qres.hsp_filter(lambda hsp: hsp.hit_id in OrthScore_dict.keys())

        # for those without any hit match to the domain
        if len(i) > 0:
            i.sort(
                key=lambda hit: OrthScore_dict[hit.id]["precision"], reverse=True)

            max_dict = OrthScore_dict[i.hits[0].id]

            # set the QueryResult attribution
            setattr(i, "orthoID", max_dict["orthoID"])
            setattr(i, "orthoID_weight", max_dict["precision"])
            q_list.append(i)
    if cf.get("filter.config") and cf["filter.config"]:
        filter_path = check_path_existence(basedir + cf["filter.config"])
        filter_res = [kofam_filter(config_file=filter_path, qres=i) for i in q_list]
        
    q_list = [i for i in filter_res if len(i) > 0]
    return q_list
