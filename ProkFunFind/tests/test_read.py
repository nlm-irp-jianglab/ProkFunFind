from ProkFunFind.read import Gene
import pytest
from Bio.SeqFeature import FeatureLocation
# from testfixtures import TempDirectory
#
#
# d = TempDirectory()
# with open(d.path+'/test.txt') as f:
#     f.write('eeee')
#
# print(d.path)
#
# for line in d.read('test.txt'):
#     print(line)
# d.cleanup()

class TestGene():
    g = Gene(contig='contig1', id='gene1', type='CDS', qualifiers='', location=FeatureLocation(1,150))


    def test_Gene_id(self):
        assert self.g.id == 'gene1'


    def test_Gene_contigid(self):
        assert self.g.contig == 'contig1'


    def test_Gene_type(self):
        assert self.g.type == 'CDS'


    def test_Gene_start(self):
        assert self.g.location.start == 1


    def test_Gene_end(self):
        assert self.g.location.end == 150


class TestContig():
    contig = {'contig1': Contig(id='contig1', seq='ATTAGAGAACGGTCGTAACATTATCGGTGGTTCTCTAACTACTATCAGTACCCATGACTCGACTCTGCCGCAGCTACCTATCGCCTGAAAGCCAGTTGGTGTTAAGGAGTGCTCTGTCCAGGACAACACGCGTAG',
        genes={'gene1':Gene(contig='contig1', id='gene1', type='CDS', qualifiers='', location=FeatureLocation(1,150))})}

    def test_contig_seq(self):
        assert self.contig.seq == 'ATTAGAGAACGGTCGTAACATTATCGGTGGTTCTCTAACTACTATCAGTACCCATGACTCGACTCTGCCGCAGCTACCTATCGCCTGAAAGCCAGTTGGTGTTAAGGAGTGCTCTGTCCAGGACAACACGCGTAG'


    def test_contig_id(self):
        assert self.contig.id == 'contig1'


    def test_contig_gene(self):
        
