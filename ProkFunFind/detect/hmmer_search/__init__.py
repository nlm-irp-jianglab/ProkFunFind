import os
from configparser import ConfigParser
from typing import IO, List, Union
import subprocess

from Bio.SearchIO._model.query import QueryResult
from Bio import SearchIO

from ProkFunFind.toolkit.utility import *
from ProkFunFind.detect.hmmer_search import *
from ProkFunFind.detect.hmmer_search.hmmer_filter import hmmer_filter


def pipeline(config: dict,
             protein_file: Union[str, IO],
             outprefix: str,
             basedir: str, OrthScore_dict: dict, q_list: dict) -> List[QueryResult]:
    """Run HMMER based search

       Arguments:
           config: A configuration dictionary
           protein_file: Input protein fasta file
           outprefix: path and prefix name of the output files
           basedir: path to config file directory
           OrthScore_dict: parsed ortholog table dictionary
           q_list: list of QueryResult objects or empty list

       Returns:
           q_list: An updated list of QueryResults
    """
    # 1. Read the query files
    query_path = check_path_existence(basedir+config['hmmer']['hmmer.query'])

    tool_format_dict = dict({'hmmsearch':"hmmsearch3-domtab",'hmmscan':"hmmscan3-domtab",'phmmer':"phmmer3-domtab"});
    outfmt = tool_format_dict[config['hmmer']['hmmer.exec']]

    # 2. Format the hmmer command
    cmd = [
        config['hmmer']['hmmer.exec'],
        "-o",
        outprefix + ".out",
        "--tblout",
        outprefix + ".tblout",
        "--domtblout",
        outprefix + ".domtblout",
        query_path,
        protein_file
        ]

    if config['hmmer'].get('hmmer.threads'):
        cmd.insert(1,"--cpu")
        cmd.insert(2,config['hmmer']['hmmer.threads'])

    # 3. run the hmmer command
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))

    # 4. read the output and apply the filtering
    qresults = SearchIO.parse(outprefix + ".domtblout", outfmt)

    if config['filter']:
        filter_res = [hmmer_filter(config=config, qres=i, basedir=basedir) for i in qresults]
        tmp_list = [i for i in filter_res if len(i) > 0]
    else:
        tmp_list = [i for i in qresults if len(i) > 0]

    # 5. Add the queries to the overall q_list
    for i in tmp_list:
        # for those without any hit match to the domain
        if len(i) > 0:
            # sort the hits by precision
            i.sort(key=lambda hit: OrthScore_dict[hit.id]['precision'], reverse=True)

            max_dict = OrthScore_dict[i.hits[0].id]
            # set the QueryResult attributes
            setattr(i, "orthoID", max_dict['orthoID'])
            setattr(i, "orthoID_weight", max_dict['precision'])
            setattr(i, "detect_tool", "hmmer")
            q_list.append(i)
    return q_list
