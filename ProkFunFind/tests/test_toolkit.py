import pytest
from os import path
import configparser
from ProkFunFind.toolkit.utility import (check_path_existence, read_config,
                                        find_file_in_folder, read2orthoDict)


test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test-data",
                          "utility")


class TestFileIdentification:

    def test_check_path_existence_abspath(self):
        f = check_path_existence(test_dir_path+'/mockin.tsv')
        assert f == test_dir_path+'/mockin.tsv'

    def test_check_path_existence_relpath(self):
        f = check_path_existence('./test-data/utility/mockin.tsv')
        assert f == test_dir_path+'/mockin.tsv'

    def test_check_path_existence_nonexist(self):
        with pytest.raises(OSError):
            check_path_existence('./test-data/utility/notreal.tsv')

    def test_find_file_in_folder_exists(self):
        fl = find_file_in_folder(folder=test_dir_path, pattern='mockin.tsv')
        assert fl == [test_dir_path+'/mockin.tsv']

    def test_find_file_in_folder_nonexist(self):
        fl = find_file_in_folder(folder=test_dir_path, pattern='notreal.tsv')
        assert fl == []


class TestFileParsing:

    def test_read_config(self):
        d = read_config(test_dir_path+'/config.ini')
        config = configparser.ConfigParser()
        config['main'] = {'cluster.tool': 'DBSCAN'}
        config['DBSCAN'] = {'cluster.eps': 4}
        assert d == config

    def test_read2orthoDict_justortho(self):
        d = read2orthoDict(ortho_pair_file=test_dir_path+'/orthogenes.tsv')
        assert d == {'gene1': {'orthoID': 'ortho1', 'precision': 1},
                     'gene2': {'orthoID': 'ortho2', 'precision': 1},
                     'gene3': {'orthoID': 'ortho3', 'precision': 1}}

    def test_read2orthoDict_prec(self):
        d = read2orthoDict(ortho_pair_file=test_dir_path+'/orthoprec.tsv')
        assert d == {'gene1': {'orthoID': 'ortho1', 'precision': 1},
                     'gene2': {'orthoID': 'ortho2', 'precision': 2},
                     'gene3': {'orthoID': 'ortho3', 'precision': 3}}

    # def test_read2orthoDict_mix(self):
    #     d = read2orthoDict(ortho_pair_file=test_dir_path+'/orthomix.tsv')
    #     assert d == {'gene1':{'orthoID':'ortho1', 'precision':1},
    #                 'gene2':{'orthoID':'ortho2', 'precision':2},
    #                 'gene3':{'orthoID':'ortho3', 'precision':1}}
