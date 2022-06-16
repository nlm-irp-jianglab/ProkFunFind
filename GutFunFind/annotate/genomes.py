import csv
import subprocess

from Bio.SeqFeature import SeqFeature, FeatureLocation

from GutFunFind.read import GetGenomeFromGFF, GetGenomeFromGzipGFF


def extract_gene(contig, location, genType, translate=False):
    """ Extract gene sequence given a contig and location

        Arguments:
            contig:
            location:
            genType:
            translate:

        Returns:
            None
    """
    f = SeqFeature(FeatureLocation(location.start, location.end,
                   strand=location.strand), type=genType)
    seq = f.extract(contig.seq)
    if translate:
        return seq.translate(table=11, cds=False, to_stop=True)
    else:
        return seq


def create_genome_object(gff_file=None, fna_file=None, gff_gz_file=None):
    """

        Arguments:
            gff_file:
            fna_file:
            gff_gz_file:

        Returns:
            None
    """
    if gff_gz_file:
        genomeObj = GetGenomeFromGzipGFF(gzipfile=gff_gz_file)
    else:
        genomeObj = GetGenomeFromGFF(gff3file=gff_file, fnafile=fna_file)
    [ct.sort() for ct in genomeObj.contigs.values()]
    return genomeObj


def parse_gtab(genome_tab: str):
    """

        Arguments:
            genome_tab

        Returns:
            None
    """
    # Add in check to make sure that file exists.
    gtab = {}
    for row in csv.reader(open(genome_tab), delimiter="\t"):
        gtab[row[0]] = row[1]

    return gtab
