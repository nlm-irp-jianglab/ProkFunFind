#!/usr/bin/env python
# -*- coding: utf-8 -*-

import GutFunFind
from GutFunFind import examine
from GutFunFind import cluster
from GutFunFind import detect
from GutFunFind import report
from GutFunFind.read import *

import os
import sys
from Bio import SeqIO
from typing import Dict, IO, List, Set, OrderedDict, Callable
from Bio.SeqRecord import SeqRecord
import subprocess
from configparser import ConfigParser


switcher = {
    'blast': 'blast_search',
    'interproscan': 'ipr_search',
    'hmmer': 'hmmer_search',
    'DBSCAN': 'DBSCAN'
}


def module_name(arg: str) -> str:
    """TODO: return the module name for name

    :arg: TODO
    :returns: TODO

    """
    func = switcher.get(arg)
    return(func)


# Write a funtion pipeline for function of interest
def retrieve_function_pipeline(database: str, fun_name: str) -> Callable:

    # 1. Obtain the configuration and check the exec as well as the database
    # required
    path_to_fun = database + fun_name

    # 1.1 check the existance of the function
    if not os.path.exists(path_to_fun):
        print("Function {} doesn't exist in GutFun".format(fun_name))
        exit()

    # 1.2 check the existance of function configuration
    config = ConfigParser()
    config.read(path_to_fun + "/config.ini")

    import importlib
    # 2. Run the correct detect tool
    detect_tool = config.get('main', 'detect.tool')
    detect_cf_path = config.get('main', 'detect.config')
    detect_module = importlib.import_module(
        "GutFunFind.detect" + "." + module_name(detect_tool), package=None)

    # 3. Run the cluster method
    cluster_tool = config.get('main', 'cluster.tool')
    cluster_cf_path = config.get('main', 'cluster.config')
    cluster_module = importlib.import_module(
        "GutFunFind.cluster" + "." + module_name(cluster_tool), package=None)

    system_file = config.get('main', 'system.file')

    def function_analysis(genome_prefix: str, outprefix: str) -> Genome:

        # varaible for the path for genome(fna),annotations(gff),proteins(faa)
        gff_file = genome_prefix + ".gff"
        fna_file = genome_prefix + ".fna"
        faa_file = genome_prefix + ".faa"

        # create genome object from genome and annotations file and sort the
        # genes
        genomeObj = GetGenomeFromGFF(gff3file=gff_file, fnafile=fna_file)
        [ct.sort() for ct in genomeObj.contigs.values()]

        # detect function related gene
        detect_list = detect_module.pipeline(
            config_file=detect_cf_path,
            protein_file=faa_file,
            outprefix=outprefix)

        # attache the detect result to genome object
        for query in detect_list:
            setattr(genomeObj.genes[query.id], detect_tool, query)

        # identify gene cluster at genome object
        genomeObj = cluster_module.pipeline(
            config_file=cluster_cf_path,
            genome_object=genomeObj,
            detect_tool=detect_tool)

        # examine the genome to see if it contain the function and annotate
        # related genes
        (system_dict, status, genomeObj) = examine.pipeline(
            system_file=system_file, genome_object=genomeObj, detect_tool=detect_tool)

        report.report_all(
            system_dict=system_dict,
            status=status,
            genomeObj=genomeObj,
            outprefix=outprefix,
            detect_tool=detect_tool,
            cluster_tool=cluster_tool)

        return (system_dict, status, genomeObj)

    return function_analysis


# if __name__ == "__main__":

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Identify genes related function of interest in genome')
    parser.add_argument(
        '-g',
        '--genomeprefix',
        help='The prefix of genome',
        required=True,
        dest='genome_prefix',
        metavar='')
    parser.add_argument(
        '-b',
        '--databasedir',
        help='The base dir of function',
        required=True,
        dest='database',
        metavar='')
    parser.add_argument(
        '-f',
        '--function',
        help='Name of the function',
        required=True,
        dest='fun_name',
        metavar='')
    parser.add_argument(
        '-o',
        '--outputprefix',
        help='The output prefix',
        required=True,
        dest='outprefix',
        metavar='')

    args = parser.parse_args()

    #database = "/home/jiangx6/data/10.GutFun/01.GutFunFind/GutFunFind/data/"
    #fun_name = "Mucic_and_Saccharic_Acid"
    detect_fun = retrieve_function_pipeline(
        database=args.database, fun_name=args.fun_name)

    #genome_prefix = "/home/jiangx6/data/10.GutFun/01.GutFunFind/GutFunFind/test/MGYG-HGUT-03390"
    #outprefix = "/home/jiangx6/data/10.GutFun/01.GutFunFind/GutFunFind/test/output"
    (system_dict, status, genomeObj) = detect_fun(
        genome_prefix=args.genome_prefix, outprefix=args.outprefix)

    if status:
        print(
            "Detect function:{fun_name} in genome {genome_prefix}\n{mp} out of {m} mandatory components present\n{ap} out of {a} accessory components present ".format(
                genome_prefix=args.genome_prefix,
                fun_name=system_dict["name"],
                mp=system_dict["completeness"]["mandatory_presence"],
                m=system_dict["completeness"]["mandatory"],
                ap=system_dict["completeness"]["accessory_presence"],
                a=system_dict["completeness"]["accessory"]))
    else:
        print(
            "Failed to detect function:{fun_name} in genome {genome_prefix}\n{mp} out of {m} mandatory components present\n{ap} out of {a} accessory components present ".format(
                genome_prefix=args.genome_prefix,
                fun_name=system_dict["name"],
                mp=system_dict["completeness"]["mandatory_presence"],
                m=system_dict["completeness"]["mandatory"],
                ap=system_dict["completeness"]["accessory_presence"],
                a=system_dict["completeness"]["accessory"]))
