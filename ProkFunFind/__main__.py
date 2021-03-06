#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import importlib
import logging

from typing import Callable, Tuple, Dict
from configparser import ConfigParser
from argparse import ArgumentParser

from ProkFunFind import examine
from ProkFunFind import report
from ProkFunFind.read import Genome, GetGenomeFromGFF
from ProkFunFind.toolkit.utility import (find_file_in_folder,
                                         check_path_existence,
                                         read2orthoDict)
from ProkFunFind.annotate.genomes import parse_gtab

logging.basicConfig(level=logging.DEBUG)

switcher = {
    'blast': 'blast_search',
    'interproscan': 'ipr_search',
    'hmmer': 'hmmer_search',
    'kofamscan': 'kofam_search',
    'emapper': 'emap_search',
    'DBSCAN': 'DBSCAN'
}


def module_name(arg: str) -> str:
    """TODO: return the module name for name

    :arg: TODO
    :returns: TODO

    """
    if not switcher.get(arg):
        sys.exit("can not find module {}".format(arg))
    else:
        return switcher.get(arg)


# Write a funtion pipeline for function of interest for a individual genome
def retrieve_function_pipeline(database: str, fun_name: str, args) -> Callable:

    # 1. Parse configuration files and search files
    # 1.1 Obtain the configuration and check the exec as well as the database
    # required
    logging.info("Checking configuration files")
    path_to_fun = database.rstrip("/") + "/" + fun_name + "/"

    # 1.2 check the existance of the function
    if not os.path.exists(path_to_fun):
        sys.exit("Function {} doesn't exist in GutFun".format(fun_name))

    # 1.3 check the existance of function configuration
    config = ConfigParser()
    config_path = path_to_fun + "config.ini"
    config_path = check_path_existence(config_path)

    config.read(config_path)

    # 1.4 Parse genome search table
    search_list = []
    gids = parse_gtab(args.gtab)
    for genome, p in gids.items():
        # fna_path = None
        prefix = p+"/"+genome
        search_list.append(prefix)

    # 1.5 Parse ortholog search table
    ortho_file = check_path_existence(path_to_fun +
                                      config['main']['search_terms'])
    OrthScore_dict, search_approaches \
        = read2orthoDict(ortho_pair_file=ortho_file)

    # 2. Check for annotation file existence for all requested searches
    for detect_tool in search_approaches:
        if detect_tool in ['interproscan', 'kofamscan', 'emapper']:
            for genome, p in gids.items():
                if detect_tool == "interproscan":
                    # Need to handle both tsv and xml outputs.
                    check_path_existence(
                        p + "/" + genome +
                        config['interproscan']['annot_suffix'])
                elif detect_tool == "kofamscan":
                    check_path_existence(
                        p + "/" + genome + config['kofamscan']['annot_suffix'])
                elif detect_tool == "emapper":
                    check_path_existence(
                        p + '/' + genome +
                        config['emapper']['annot_suffix'])

    # 3. Set up clustering and parse system file
    cluster_tool = config.get('main', 'cluster.tool')

    cluster_module = importlib.import_module(
        "ProkFunFind.cluster" + "." + module_name(cluster_tool), package=None)

    system_file = path_to_fun + config.get('main', 'system.file')
    system_file = check_path_existence(system_file)

    # Function to run analysis on a given genome
    def function_analysis(genome_prefix: str,
                          outprefix: str) -> Tuple[Dict, int, Genome]:

        # get the paths to the the GFF and fna files
        gff_file = genome_prefix + config['main']['gff_suffix']
        gff_file = check_path_existence(gff_file)
        fna_file = genome_prefix + config['main']['fna_suffix']
        fna_file = check_path_existence(fna_file)

        outprefix = os.path.abspath(outprefix)
        out_base = os.path.basename(outprefix)
        out_dir = os.path.dirname(outprefix)
        if not out_base and out_dir:
            sys.exit(
                     "Please provide prefix of output file, Not directory \
                     name {}".format(out_dir))
        else:
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

        # create genome object from genome and gff file and sort the genes
        genomeObj = GetGenomeFromGFF(gff3file=gff_file, fnafile=fna_file)
        [ct.sort() for ct in genomeObj.contigs.values()]
        detect_list = []

        logging.info('Searching for function')
        # Run detect methods
        if "interproscan" in search_approaches:
            detect_module = importlib.import_module(
                "ProkFunFind.detect" + "." +
                module_name('interproscan'), package=None)
            xml_file = genome_prefix + ".xml"
            if not os.path.isfile(xml_file):
                tsv_file = genome_prefix + \
                    config['interproscan']['annot_suffix']
                if not os.path.isfile(tsv_file):
                    sys.exit(
                        "Neither file {} nor file {} exists".format(
                            xml_file, tsv_file))
                else:
                    detect_list = detect_module.pipeline(
                        config=config,
                        in_file=tsv_file,
                        fmt="tsv",
                        basedir=path_to_fun,
                        OrthScore_dict=OrthScore_dict['interproscan'],
                        q_list=detect_list)
            else:
                detect_list = detect_module.pipeline(
                    config=config,
                    in_file=xml_file,
                    fmt="xml",
                    basedir=path_to_fun,
                    OrthScore_dict=OrthScore_dict['interproscan'],
                    q_list=detect_list)
        if "kofamscan" in search_approaches:
            detect_module = importlib.import_module(
                "ProkFunFind.detect" + "." +
                module_name('kofamscan'), package=None)
            kofam_file = genome_prefix + config['kofamscan']['annot_suffix']
            kofam_file = check_path_existence(kofam_file)
            detect_list = detect_module.pipeline(
                config=config,
                in_file=kofam_file,
                basedir=path_to_fun,
                OrthScore_dict=OrthScore_dict['kofamscan'],
                q_list=detect_list)
        if "emapper" in search_approaches:
            detect_module = importlib.import_module(
                "ProkFunFind.detect" + "." +
                module_name('emapper'), package=None)
            emap_file = genome_prefix + config['emapper']['annot_suffix']
            emap_file = check_path_existence(emap_file)
            detect_list = detect_module.pipeline(
                config=config,
                in_file=emap_file,
                basedir=path_to_fun,
                OrthScore_dict=OrthScore_dict['emapper'],
                q_list=detect_list)

        if 'blast' in search_approaches:
            detect_module = importlib.import_module(
                "ProkFunFind.detect" + "." +
                module_name('blast'), package=None)
            faa_file = genome_prefix + config['main']['faa_suffix']
            faa_file = check_path_existence(faa_file)
            detect_list = detect_module.pipeline(
                config=config,
                protein_file=faa_file,
                outprefix=outprefix,
                basedir=path_to_fun,
                OrthScore_dict=OrthScore_dict['blast'],
                q_list=detect_list)
        if 'hmmer' in search_approaches:
            detect_module = importlib.import_module(
                "ProkFunFind.detect" + "." +
                module_name('hmmer'), package=None)
            faa_file = genome_prefix + config['main']['faa_suffix']
            faa_file = check_path_existence(faa_file)
            detect_list = detect_module.pipeline(
                config=config,
                protein_file=faa_file,
                outprefix=outprefix,
                basedir=path_to_fun,
                OrthScore_dict=OrthScore_dict['hmmer'],
                q_list=detect_list)

        # Attach search results to genome object
        for query in detect_list:
            setattr(genomeObj.genes[query.id], query.detect_tool, query)

        logging.info("Identifying gene clusters")
        # identify gene cluster at genome object
        genomeObj = cluster_module.pipeline(
            config=config,
            genome_object=genomeObj,
            detect_tools=search_approaches)

        # examine the genome to see if it contain the function and annotate
        # related genes
        logging.info("Summarizing function presence and genes")
        (system_dict, status, genomeObj) = examine.pipeline(
            system_file=system_file, genome_object=genomeObj,
            detect_tools=search_approaches)
        genome_name = genome_prefix.split("/")[len(genome_prefix.split("/"))-1]
        report.report_all(
            system_dict=system_dict,
            status=status,
            genomeObj=genomeObj,
            outprefix=outprefix+'.'+genome_name,
            detect_tools=search_approaches,
            cluster_tool=cluster_tool)

        if status:
            print(
                "Detected function:{fun_name} in genome {genome_prefix}\n{mp}" \
                " out of {m} essential components present\n{ap} out of {a}" \
                " nonessential components present ".format(
                    genome_prefix=genome_prefix,
                    fun_name=system_dict['name'],
                    mp=system_dict['completeness']['essential_presence'],
                    m=system_dict['completeness']['essential'],
                    ap=system_dict['completeness']['nonessential_presence'],
                    a=system_dict['completeness']['nonessential']))
        else:
            print(
                "Failed to detect function:{fun_name} in" \
                " genome {genome_prefix}\n{mp} out of {m} essential" \
                " components present\n{ap} out of {a} nonessential" \
                " components present".format(
                    genome_prefix=genome_prefix,
                    fun_name=system_dict['name'],
                    mp=system_dict['completeness']['essential_presence'],
                    m=system_dict['completeness']['essential'],
                    ap=system_dict['completeness']['nonessential_presence'],
                    a=system_dict['completeness']['nonessential']))

    return function_analysis, search_list


def main_individual(args):
    detect_fun, search_list = retrieve_function_pipeline(
        database=args.database, fun_name=args.fun_name, args=args)

    for prefix in search_list:
        detect_fun(genome_prefix=prefix, outprefix=args.outprefix)


def main():
    parser = ArgumentParser(
        description="Identify genes related function of interest")
    # subparsers = parser.add_subparsers(dest="command")

    # parser_rep = subparsers.add_parser(
        # "rep", help="Analyze an individual genome")

    parser.add_argument(
        "-b",
        "--databasedir",
        help="The base dir of function",
        required=True,
        dest="database",
        metavar="")
    parser.add_argument(
        "-f",
        "--function",
        help="Name of the function",
        required=True,
        dest="fun_name",
        metavar="")
    parser.add_argument(
        "-o",
        "--outputprefix",
        help="The output file name prefix",
        required=True,
        dest="outprefix",
        metavar="")
    parser.add_argument(
        "--gtab",
        help="Table of genomes to search",
        required=True,
        dest="gtab",
        metavar="")
    parser.set_defaults(func=main_individual)

    options = parser.parse_args()
    options.func(options)
