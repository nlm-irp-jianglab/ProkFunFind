import os
from configparser import ConfigParser
from typing import IO, List, Union

from Bio.SearchIO._model.query import QueryResult
from Bio import SearchIO

from GutFunFind.toolkit.base import *
from GutFunFind.detect.hmmer_search import *
from GutFunFind.detect.hmmer_search.hmmer_filter import hmmer_filter


def pipeline(config_file: Union[str, IO],
             protein_file: Union[str, IO],
             outprefix: str) -> List[QueryResult]:
    # 1. Read the configuration file into configuration object
    cf = read_config(config_file)["hmmer"]
    basedir = os.path.dirname(os.path.abspath(config_file))+"/"

    # 2 generate command line to run

    query_path = check_path_existence(basedir+cf["hmmer.query"])

    tool_format_dict = dict({"hmmsearch":"hmmsearch3-domtab","hmmscan":"hmmscan3-domtab","phmmer":"phmmer3-domtab"});
    outfmt = tool_format_dict[cf["hmmer.exec"]]

    cmd = [
        cf["hmmer.exec"],
        "-o",
        outprefix + ".out",
        "--tblout",
        outprefix + ".tblout",
        "--domtblout",
        outprefix + ".domtblout",
        query_path,
        protein_file
        ]

    if cf.get("hmmer.threads"):
        cmd.insert(1,"--cpu")
        cmd.insert(2,cf["hmmer.threads"])

    # 3 run command line
    # check the protein_file+".blast.xml" existence
    res = execute(cmd)
    if res.return_code:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))

    # 4 read the output from previous step
    qresults = SearchIO.parse(outprefix + ".domtblout", outfmt)

    if cf.get("filter.config") and cf["filter.config"]:
        filter_path = check_path_existence(basedir + cf["filter.config"])
        filter_res = [hmmer_filter(config_file=filter_path, qres=i) for i in qresults]
        q_list = [i for i in filter_res if len(i) > 0]
    else:
        q_list = [i for i in qresults if len(i) > 0]


    ortho_file = check_path_existence(basedir + cf["orthoID_domain_precision"])
    OrthScore_dict = read2orthoDict(ortho_pair_file=ortho_file)

    # 4. Process all QueryResult
    a_list = []

    for i in q_list:
        # for those without any hit match to the domain
        if len(i) > 0:
            # sort the hits by precision
            # max_dict = sorted([OrthScore_dict[x.id] for x in i], key = lambda i: i['precision'],reverse=True)[0]
            i.sort(key=lambda hit: OrthScore_dict[hit.id]["precision"], reverse=True)

            max_dict = OrthScore_dict[i.hits[0].id]

            # set the QueryResult attribution
            setattr(i, "orthoID", max_dict["orthoID"])
            setattr(i, "orthoID_weight", max_dict["precision"])
            a_list.append(i)
    return a_list
