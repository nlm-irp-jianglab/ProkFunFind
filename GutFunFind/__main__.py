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
from typing import Callable
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
    path_to_fun = database.rstrip("/")+"/" + fun_name

    # 1.1 check the existance of the function
    if not os.path.exists(path_to_fun):
        print("Function {} doesn't exist in GutFun".format(fun_name))
        exit()

    # 1.2 check the existance of function configuration
    config = ConfigParser()
    config_path = path_to_fun + "/config.ini"

    if not os.path.exists(config_path):
        sys.exit('Can not find {}. Please check the path of config.ini'.format(config_path))
    config.read(config_path)

    import importlib
    # 2. Run the correct detect tool
    detect_tool = config.get('main', 'detect.tool')
    detect_cf_path = config.get('main', 'detect.config')

    if not os.path.exists(detect_cf_path):
        sys.exit('Can not find {}. Please check the detect.inc file'.format(detect_cf_path))

    detect_module = importlib.import_module(
        "GutFunFind.detect" + "." + module_name(detect_tool), package=None)

    # 3. Run the cluster method
    cluster_tool = config.get('main', 'cluster.tool')
    cluster_cf_path = config.get('main', 'cluster.config')

    if not os.path.exists(cluster_cf_path):
        sys.exit('Can not find {}. Please check the cluster.inc file'.format(cluster_cf_path))

    cluster_module = importlib.import_module(
        "GutFunFind.cluster" + "." + module_name(cluster_tool), package=None)

    system_file = config.get('main', 'system.file')

    if not os.path.exists(system_file):
        sys.exit('Can not find {}. Please check the system_file'.format(system_file))


    def function_analysis(genome_prefix: str, outprefix: str) -> Genome:

        # varaible for the path for genome(fna),annotations(gff),proteins(faa)
        gff_file = genome_prefix + ".gff"
        if not os.path.exists(gff_file):
            sys.exit('Can not find {}. Please check the gff file'.format(gff_file))
        
        fna_file = genome_prefix + ".fna"
        if not os.path.exists(fna_file):
            sys.exit('Can not find {}. Please check the fna file'.format(fna_file))

        # create genome object from genome and annotations file and sort the
        # genes
        genomeObj = GetGenomeFromGFF(gff3file=gff_file, fnafile=fna_file)
        [ct.sort() for ct in genomeObj.contigs.values()]

        # detect function related gene
        if detect_tool == "interproscan":
            xml_file = genome_prefix + ".xml"
            if not os.path.isfile(xml_file):
                tsv_file = genome_prefix + "_InterProScan.tsv"
                if not os.path.isfile(tsv_file):
                    sys.exit('Neither file {} nor file {} exists'.format(xml_file,tsv_file))
                else:
                    detect_list = detect_module.pipeline(
                            config_file=detect_cf_path,
                            in_file= tsv_file,
                            fmt = "tsv")
            else:
                detect_list = detect_module.pipeline(
                        config_file=detect_cf_path,
                        in_file=xml_file,
                        fmt = "xml")
        else:
            faa_file = genome_prefix + ".faa"
            if not os.path.exists(faa_file):
                sys.exit('Can not find {}. Please check the faa file'.format(faa_file))
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
            "Detect function:{fun_name} in genome {genome_prefix}\n{mp} out of {m} essential components present\n{ap} out of {a} nonessential components present ".format(
                genome_prefix=args.genome_prefix,
                fun_name=system_dict["name"],
                mp=system_dict["completeness"]["essential_presence"],
                m=system_dict["completeness"]["essential"],
                ap=system_dict["completeness"]["nonessential_presence"],
                a=system_dict["completeness"]["nonessential"]))
    else:
        print(
            "Failed to detect function:{fun_name} in genome {genome_prefix}\n{mp} out of {m} essential components present\n{ap} out of {a} nonessential components present ".format(
                genome_prefix=args.genome_prefix,
                fun_name=system_dict["name"],
                mp=system_dict["completeness"]["essential_presence"],
                m=system_dict["completeness"]["essential"],
                ap=system_dict["completeness"]["nonessential_presence"],
                a=system_dict["completeness"]["nonessential"]))
