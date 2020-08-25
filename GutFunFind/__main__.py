#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import glob
import os
import sys

from GutFunFind import examine
from GutFunFind import cluster
from GutFunFind import detect
from GutFunFind import report
from GutFunFind.read import GetGenomeFromGFF, Genome, Roarycsv2pangenome, GetGenomeFromGzipGFF

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from typing import Callable, List, Tuple, Dict, AnyStr
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
    if not switcher.get(arg):
        sys.exit("can not find module {}".format(arg))
    else:
        return switcher.get(arg)


def find_file_in_folder(folder: AnyStr, pattern: AnyStr) -> List:
    import fnmatch
    fileList = []
    for dName, sdName, fList in os.walk(folder):
        for fileName in fList:
            if fnmatch.fnmatch(fileName, pattern):
                fileList.append(os.path.join(dName, fileName))
    return fileList

def check_path_existence(path):
    abspath = os.path.abspath(path)
    if not os.path.exists(abspath):
        sys.exit("Can not find {}! Please check!".format(abspath))
    return(abspath)

# Write a funtion pipeline for function of interest for a individual genome
def retrieve_function_pipeline(database: str, fun_name: str) -> Callable:

    # 1. Obtain the configuration and check the exec as well as the database
    # required
    path_to_fun = database.rstrip("/") + "/" + fun_name

    # 1.1 check the existance of the function
    if not os.path.exists(path_to_fun):
        sys.exist("Function {} doesn't exist in GutFun".format(fun_name))

    # 1.2 check the existance of function configuration
    config = ConfigParser()
    config_path = path_to_fun + "/config.ini"
    config_path = check_path_existence(config_path)

    config.read(config_path)

    import importlib
    # 2. Run the correct detect tool
    detect_tool = config.get('main', 'detect.tool')
    detect_cf_path = config.get('main', 'detect.config')
    detect_cf_path = check_path_existence(detect_cf_path)

    detect_module = importlib.import_module(
        "GutFunFind.detect" + "." + module_name(detect_tool), package=None)

    # 3. Run the cluster method
    cluster_tool = config.get('main', 'cluster.tool')
    cluster_cf_path = config.get('main', 'cluster.config')
    cluster_cf_path = check_path_existence(cluster_cf_path)

    cluster_module = importlib.import_module(
        "GutFunFind.cluster" + "." + module_name(cluster_tool), package=None)

    system_file = config.get('main', 'system.file')
    system_file = check_path_existence(system_file)

    def function_analysis(genome_prefix: str,
                          outprefix: str) -> Tuple[Dict, int, Genome]:

        # varaible for the path for genome(fna),annotations(gff),proteins(faa)
        gff_file = genome_prefix + ".gff"
        gff_file = check_path_existence(gff_file)
        fna_file = genome_prefix + ".fna"
        fna_file = check_path_existence(fna_file)


        outprefix = os.path.abspath(outprefix)
        out_base = os.path.basename(outprefix)
        out_dir = os.path.dirname(outprefix)
        if not out_base and out_dir:
            sys.exit(
                "Please provide prefix of output file, Not directory name {}".format(out_dir))
        else:
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

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
                    sys.exit(
                        'Neither file {} nor file {} exists'.format(
                            xml_file, tsv_file))
                else:
                    detect_list = detect_module.pipeline(
                        config_file=detect_cf_path,
                        in_file=tsv_file,
                        fmt="tsv")
            else:
                detect_list = detect_module.pipeline(
                    config_file=detect_cf_path,
                    in_file=xml_file,
                    fmt="xml")
        else:
            faa_file = genome_prefix + ".faa"
            faa_file = check_path_existence(faa_file)
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


def main_individual(args):

    detect_fun = retrieve_function_pipeline(
        database=args.database, fun_name=args.fun_name)

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


# Write a funtion pipeline for function of interest for pangenome
def retrieve_function_pipeline_pan(database: str, fun_name: str) -> Callable:

    # 1. Obtain the configuration and check the exec as well as the database
    # required
    path_to_fun = database.rstrip("/") + "/" + fun_name

    # 1.1 check the existance of the function
    if not os.path.exists(path_to_fun):
        sys.exist("Function {} doesn't exist in GutFun".format(fun_name))

    # 1.2 check the existance of function configuration
    config = ConfigParser()
    config_path = path_to_fun + "/config.ini"
    config_path = check_path_existence(config_path)

    config.read(config_path)

    import importlib
    # 2. Run the correct detect tool
    detect_tool = config.get('main', 'detect.tool')
    detect_cf_path = config.get('main', 'detect.config')
    detect_cf_path = check_path_existence(detect_cf_path)

    detect_module = importlib.import_module(
        "GutFunFind.detect" + "." + module_name(detect_tool), package=None)

    # 3. Run the cluster method
    cluster_tool = config.get('main', 'cluster.tool')
    cluster_cf_path = config.get('main', 'cluster.config')
    cluster_cf_path = check_path_existence(cluster_cf_path)

    cluster_module = importlib.import_module(
        "GutFunFind.cluster" + "." + module_name(cluster_tool), package=None)

    system_file = config.get('main', 'system.file')
    system_file = check_path_existence(system_file)

    def function_analysis(pangenome_path: str, outprefix: str, folder: str):

        roaryfile = pangenome_path + "/genes_presence-absence_locus.csv"
        roaryfile = check_path_existence(roaryfile)
        genome_prefix = pangenome_path + "/pan-genome"

        outprefix = os.path.abspath(outprefix).rstrip("/") + "/"
        if os.path.exists(outprefix):
            sys.exit(
                "Directory name {} exists! Please delete the directory first!".format(outprefix))
        else:
            os.makedirs(outprefix)

        # detect function related gene
        if detect_tool == "interproscan":
            xml_file = genome_prefix + ".xml"
            if not os.path.isfile(xml_file):
                tsv_file = genome_prefix + "_InterProScan.tsv"
                if not os.path.isfile(tsv_file):
                    sys.exit(
                        'Neither file {} nor file {} exists'.format(
                            xml_file, tsv_file))
                else:
                    detect_list = detect_module.pipeline(
                        config_file=detect_cf_path,
                        in_file=tsv_file,
                        fmt="tsv")
            else:
                detect_list = detect_module.pipeline(
                    config_file=detect_cf_path,
                    in_file=xml_file,
                    fmt="xml")
        else:
            faa_file = genome_prefix + ".faa"
            faa_file = check_path_existence(faa_file)
            detect_list = detect_module.pipeline(
                config_file=detect_cf_path,
                protein_file=faa_file,
                outprefix=outprefix + "pan-genome")

        # obtain pan genome information and transfer detect result to the
        # genomes in the pan-genome set

        pangenome = Roarycsv2pangenome(roaryfile)
        result_dict = OrderedDict()

        for genome_name in pangenome.genomes:
            genome_paths = find_file_in_folder(
                folder=folder, pattern=genome_name + ".gff.gz")

            if len(genome_paths) == 1:
                genome_path = genome_paths[0]
            else:
                sys.exit(
                    'Can not find {f} or more than one {f} is found'.format(
                        f=genome_name + ".gff.gz"))

            genomeObj = GetGenomeFromGzipGFF(genome_path)
            [ct.sort() for ct in genomeObj.contigs.values()]

            for query in detect_list:
                gene_group = pangenome.get_genegroup(query.id)
                gene_list_to_annotate = gene_group[genome_name]
                for i in gene_list_to_annotate:
                    if i: # if i == "" will not setattr
                        setattr(genomeObj.genes[i], detect_tool, query)
                        setattr(
                            genomeObj.genes[i],
                            "pangenome_group",
                            gene_group)

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
                outprefix=outprefix + genome_name,
                detect_tool=detect_tool,
                cluster_tool=cluster_tool)

            result_dict.update({genome_name: [system_dict, status, genomeObj]})

        return result_dict

    return function_analysis


def main_pan(args):
    detect_fun = retrieve_function_pipeline_pan(
        database=args.database, fun_name=args.fun_name)
    res_dict = detect_fun(
        pangenome_path=args.pangenome_path,
        outprefix=args.outprefix,
        folder=args.folder)

    for res_key in res_dict:
        system_dict, status, genomeObj = res_dict[res_key]

        if status:
            print(
                "Detect function:{fun_name} in genome {genome_prefix}\n{mp} out of {m} essential components present\n{ap} out of {a} nonessential components present ".format(
                    genome_prefix=res_key,
                    fun_name=system_dict["name"],
                    mp=system_dict["completeness"]["essential_presence"],
                    m=system_dict["completeness"]["essential"],
                    ap=system_dict["completeness"]["nonessential_presence"],
                    a=system_dict["completeness"]["nonessential"]))
        else:
            print(
                "Failed to detect function:{fun_name} in genome {genome_prefix}\n{mp} out of {m} essential components present\n{ap} out of {a} nonessential components present ".format(
                    genome_prefix=res_key,
                    fun_name=system_dict["name"],
                    mp=system_dict["completeness"]["essential_presence"],
                    m=system_dict["completeness"]["essential"],
                    ap=system_dict["completeness"]["nonessential_presence"],
                    a=system_dict["completeness"]["nonessential"]))


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Identify genes related function of interest')
    subparsers = parser.add_subparsers(
        dest='command', title="Usage: {} <command> [options]".format(
            os.path.basename(__file__)), metavar="Commands:")

    parser_rep = subparsers.add_parser(
        "rep", help='Analyze an individual genome')

    parser_rep.add_argument(
        '-b',
        '--databasedir',
        help='The base dir of function',
        required=True,
        dest='database',
        metavar='')
    parser_rep.add_argument(
        '-f',
        '--function',
        help='Name of the function',
        required=True,
        dest='fun_name',
        metavar='')

    parser_rep.add_argument(
        '-g',
        '--genomeprefix',
        help='The prefix of genome',
        required=True,
        dest='genome_prefix',
        metavar='')
    parser_rep.add_argument(
        '-o',
        '--outputprefix',
        help='The output prefix',
        required=True,
        dest='outprefix',
        metavar='')
    parser_rep.set_defaults(func=main_individual)

    parser_pan = subparsers.add_parser(
        "pan", help='Analyze genomes of the same pangenome')

    parser_pan.add_argument(
        '-b',
        '--databasedir',
        help='The base dir of function',
        required=True,
        dest='database',
        metavar='')
    parser_pan.add_argument(
        '-f',
        '--function',
        help='Name of the function',
        required=True,
        dest='fun_name',
        metavar='')

    parser_pan.add_argument(
        '-p',
        '--pangenome',
        help='The path to the pan-genome info dir',
        required=True,
        dest='pangenome_path',
        metavar='')
    parser_pan.add_argument(
        '-g',
        '--genomeset',
        help='the dir of gzipped gff file',
        required=True,
        dest='folder',
        metavar='')
    parser_pan.add_argument(
        '-d',
        '--outputdir',
        help='The output directory',
        required=True,
        dest='outprefix',
        metavar='')
    parser_pan.set_defaults(func=main_pan)

    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    options.func(options)
