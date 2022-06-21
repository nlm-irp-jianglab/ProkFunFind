from typing import IO, List, Union
from configparser import ConfigParser

from Bio import SearchIO
from Bio.SearchIO._model.query import QueryResult

from GutFunFind.toolkit.utility import *
from GutFunFind.detect.ipr_search.interproscan_tab import InterproscanTabParser, ipr_tab_parse


def pipeline(config: dict,
             in_file: Union[str, IO],
             fmt: str,
             basedir: str, OrthScore_dict: dict, q_list: dict) -> List[QueryResult]:
    # 1. Read the configuration file into configuration object
    # cf = read_config(config_file)["interproscan"]
    # basedir = os.path.dirname(os.path.abspath(config_file))+"/"

    # 2. Read interproscan xml file
    if fmt == "xml":
        qresults = SearchIO.parse(in_file, "interproscan-xml")
    elif fmt == "tsv":
        qresults = ipr_tab_parse(in_file)

    # 3. Read score and orthoID info into dictionary
    # ortho_file = check_path_existence(basedir + config['interproscan']['orthoID_domain_precision'])
    # OrthScore_dict = read2orthoDict(ortho_pair_file=ortho_file)

    # 4. Process all QueryResult
    # q_list = []
    for qres in qresults:

        # remove QueryResult that doest not hit any domain in function-related domain list
        i = qres.hsp_filter(lambda hsp: hsp.hit_id in OrthScore_dict.keys())

        # for those without any hit match to the domain
        if len(i) > 0:
            # sort the hits by precision
            # max_dict = sorted([OrthScore_dict[x.id] for x in i], key = lambda i: i['precision'],reverse=True)[0]
            i.sort(
                key=lambda hit: OrthScore_dict[hit.id]['precision'], reverse=True)

            max_dict = OrthScore_dict[i.hits[0].id]

            # set the QueryResult attribution
            setattr(i, "orthoID", max_dict['orthoID'])
            setattr(i, "orthoID_weight", max_dict['precision'])
            setattr(i, "detect_tool", "interproscan")
            q_list.append(i)

    return q_list
