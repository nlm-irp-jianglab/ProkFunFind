from configparser import ConfigParser
import logging
import os
from collections import defaultdict
import csv
import operator
from typing import Dict, IO, List, Union, AnyStr


def read_config(config_file: str) -> ConfigParser:
    config = ConfigParser()
    config.read(config_file)
    return(config)


def find_file_in_folder(folder: AnyStr, pattern: AnyStr) -> List:
    """ Find files with a given pattern in a given file path

        Arguments:
            folder: the directory to search
            pattern: the pattern to search for

        Returns:
            A list of path for the files with a given pattern in the file path
    """
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
        # logging.warning
        raise OSError("Can not find {}! Please check!".format(abspath))
    return(abspath)


def read2orthoDict(ortho_pair_file: Union[str, IO]) -> Dict:
    OrthScore_dict = defaultdict(dict)
    with open(ortho_pair_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        header = next(csv_reader)
        col_num = len(header)
        if col_num == 3:
            OrthScore_dict[header[2]][header[1]] = \
                {'queryID': header[0], 'precision': 1}
            for row in csv_reader:
                OrthScore_dict[row[2]][row[1]] = \
                    {'queryID': row[0], 'precision': 1}
        else:
            OrthScore_dict[header[2]][header[1]] = {
                'queryID': header[0], 'precision': float(header[3])}
            for row in csv_reader:
                OrthScore_dict[row[2]][row[1]] = {
                    'queryID': row[0], 'precision': float(row[3])}
    search_set = set()
    for i in OrthScore_dict.keys():
        if i not in ['kofamscan', 'interproscan', 'emapper', 'blast', 'hmmer']:
            logging.error('Search approach {} is not supported. \
                           Please check ortho table file'.format(i))
        else:
            search_set.add(i)
    return OrthScore_dict, search_set

def parse_system_yaml(system_dict):
    ops = {
        '<=': operator.le,
        '>=': operator.ge,
        '>': operator.gt,
        '<': operator.lt,
        '==': operator.eq,
        '!=': operator.ne
    }
    search_set = set()
    OrthScore_dict  = defaultdict(dict)
    filter_dict = defaultdict(lambda: defaultdict(list))
    curr_quer = ""
    def recursedict(d, curr_quer):
        for k,v in d.items():
            if k == "queryID":
                curr_quer = v
            if isinstance(v, dict):
                recursedict(v, curr_quer)
            elif k == 'terms':
                    for gene in v:
                        search_set.add(gene['method'])
                        method = gene['method']
                        OrthScore_dict[method][gene['id']] = {'queryID': curr_quer, 'precision': gene.get('precision', 1)}
                        if 'ident_pct' in gene:
                            filter_dict[method][gene['id']].append({'attr': 'ident_pct', 'cpfun': ops['>='], 'value': float(gene['ident_pct'])})
                        if 'evalue' in gene:
                            filter_dict[method][gene['id']].append({'attr': 'evalue', 'cpfun': ops['<='], 'value':float(gene['evalue'])})
                        if 'threshold' in gene:
                            filter_dict[method][gene['id']].append({'attr': 'threshold', 'cpfun': ops['=='], 'value':float(gene['threshold'])})
            elif isinstance(v, list):
                for w in v:
                    recursedict(w, curr_quer)

    recursedict(system_dict, curr_quer)
    return OrthScore_dict, search_set, filter_dict
