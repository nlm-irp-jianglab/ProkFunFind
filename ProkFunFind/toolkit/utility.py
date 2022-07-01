from configparser import ConfigParser
import logging
import warnings
import os
from collections import defaultdict
import csv
import sys
from typing import Dict, Any, Callable, Iterable, IO, List, Optional, Union, AnyStr


# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # for import by others without wrapping, pylint: disable=unused-import
    from Bio import SearchIO


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
            A list of path for the files with a given pattern in a given file path
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

    #############################################################
    #  User can change the last column to indicate precision    #
    #############################################################

    with open(ortho_pair_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        header = next(csv_reader)
        col_num = len(header)
        if col_num == 3:
            OrthScore_dict[header[2]][header[1]] = {'orthoID': header[0], 'precision': 1}
            for row in csv_reader:
                OrthScore_dict[row[2]][row[1]] = {'orthoID': row[0], 'precision': 1}
        else:
            OrthScore_dict[header[2]][header[1]] = {
                'orthoID': header[0], 'precision': float(header[3])}
            for row in csv_reader:
                OrthScore_dict[row[2]][row[1]] = {
                    'orthoID': row[0], 'precision': float(row[3])}
    search_set = set()
    for i in OrthScore_dict.keys():
        if i not in ['kofamscan', 'interproscan', 'emapper', 'blast', 'hmmer']:
            logging.error('Search approach {} is not supported. Please check ortho table file'.format(i))
        else:
            search_set.add(i)
    return OrthScore_dict, search_set
