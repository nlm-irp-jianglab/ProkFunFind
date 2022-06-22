import os
from typing import IO, List, Union
import subprocess

from Bio.SearchIO._model.query import QueryResult
from Bio import SearchIO

from GutFunFind.toolkit.utility import *
from GutFunFind.detect.blast_search.blast_filter import blast_filter, blast_ortho


def pipeline(config: dict,
             protein_file: Union[str, IO],
             outprefix: str,
             basedir: str, q_list: dict, OrthScore_dict: dict) -> List[QueryResult]:

    # 1. Check for query file
    query_path = check_path_existence(basedir+config['blast']['blast.query'])

    # 2. Set up blast command.
    cmd = [
        config['blast']['blast.exec'],
        "-query",
        protein_file,
        "-db",
        query_path,
        "-evalue",
        config['blast']['blast.evalue'],
        "-outfmt",
        "6",
        "-out",
        outprefix + ".blast.m6"]

    if config['blast'].get('blast.threads'):
        cmd +=["-num_threads",config['blast']['blast.threads']]

    # 3. run command line
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))

    # 4. read the output from previous step
    qresults = SearchIO.parse(outprefix + ".blast.m6", "blast-tab")
    #qresults = SearchIO.parse(protein_file+".blast.xml", 'blast-xml')
    tmp_list = [i for _, i in SearchIO.to_dict(qresults).items() if len(i) > 0]

    # 5. Apply blast filtering.
    if config['filter']:
        filter_res = [blast_filter(
            config=config, qres=i, basedir=basedir) for i in tmp_list]
        res_list = [i for i in filter_res if len(i) > 0]
    else:
        res_list = tmp_list

    for i in res_list:
        q_list.append(blast_ortho(OrthScore_dict=OrthScore_dict, qres=i))

    return q_list
