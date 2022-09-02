from ProkFunFind.annotate.genomes import *
import pytest

test_dir_path = path.join(path.dirname(path.abspath(__file__)), "test-data",
                          "annotate")

class TestGtab():
    def test_parse_gtab(self)
