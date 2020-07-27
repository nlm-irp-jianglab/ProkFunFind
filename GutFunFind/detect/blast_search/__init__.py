from Bio.SearchIO._model.query import QueryResult
from Bio import SearchIO
from typing import Dict, IO, List, Set, OrderedDict, Any, Tuple, Union
from configparser import ConfigParser
from .blast_filter import blast_filter, blast_ortho
from GutFunFind.toolkit.base import *


def pipeline(config_file: Union[str,
                                IO],
             protein_file: Union[str,
                                 IO],
             outprefix: str) -> List[QueryResult]:
    # 1. Read the configuration file into configuration object
    cf = read_config(config_file)["blast"]

    # 2 generate command line to run
    # we can add "-num_threads" later
    # Right now we use the -m 6 format
    #cmd = [ cf["blast.exec"], "-query",  protein_file, "-db" ,cf["blast.query"], "-evalue", cf["blast.evalue"], "-outfmt" ,"5", "-out", outprefix+".blast.xml"]
    cmd = [
        cf["blast.exec"],
        "-query",
        protein_file,
        "-db",
        cf["blast.query"],
        "-evalue",
        cf["blast.evalue"],
        "-outfmt",
        "6",
        "-out",
        outprefix +
        ".blast.m6"]

    # 3 run command line
    # check the protein_file+".blast.xml" existence
    res = execute(cmd)
    if res.return_code:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))

    # 4 read the output from previous step
    qresults = SearchIO.parse(outprefix + ".blast.m6", 'blast-tab')
    #qresults = SearchIO.parse(protein_file+".blast.xml", 'blast-xml')
    q_list = [i for _, i in SearchIO.to_dict(qresults).items() if len(i) > 0]

    if cf["filter.config"]:
        filter_cf = read_config(cf["filter.config"])
        filter_res = [blast_filter(config=filter_cf, qres=i) for i in q_list]
        q_list = [i for i in filter_res if len(i) > 0]

    ortho_file = cf["map.ortho_pair"]
    q_list = [blast_ortho(ortho_pair_file=ortho_file, qres=i) for i in q_list]

    return q_list
