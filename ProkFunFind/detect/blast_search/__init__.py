from typing import IO, List, Union
import subprocess

from Bio.SearchIO._model.query import QueryResult
from Bio import SearchIO

from ProkFunFind.toolkit.utility import check_path_existence
from ProkFunFind.detect.blast_search.blast_filter import \
    blast_filter, blast_ortho


def pipeline(config: dict,
             protein_file: Union[str, IO],
             outprefix: str,
             basedir: str,
             q_list: dict,
             OrthScore_dict: dict) -> List[QueryResult]:
    """Main pipeline function for performing blast-based searches

       Arguments:
           config:
           protein_file:
           outprefix:
           basedir:
           q_list:
           OrthScore_dict:

       Returns:
           q_list: updated list of QueryResult objects.

    """

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
        cmd += ["-num_threads", config['blast']['blast.threads']]

    # 3. run command line
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))

    # 4. read the output from previous step
    qresults = SearchIO.parse(outprefix + ".blast.m6", "blast-tab")
    tmp_list = [i for _, i in SearchIO.to_dict(qresults).items() if len(i) > 0]

    # 5. Apply blast filtering.
    filter_res = [blast_filter(
        config=config, qres=i, basedir=basedir) for i in tmp_list]
    res_list = [i for i in filter_res if len(i) > 0]

    for i in res_list:
        q_list.append(blast_ortho(OrthScore_dict=OrthScore_dict, qres=i))

    return q_list
