import csv

from Bio.SeqFeature import SeqFeature, FeatureLocation

from ProkFunFind.read import GetGenomeFromGFF, GetGenomeFromGzipGFF


def extract_gene(contig, location, genType, translate=False):
    """ Extract gene sequence given a contig and location

        Arguments:
            contig: A ProkFunFind.annotate.contig object
            location: A Bio.SeqFeature.FeatureLocation object
            genType: A string specifying the gene type property
            translate: Bool, if true then translate the gene seq

        Returns:
            seq: A nucleotide or amino acid Bio.Seq.Seq object
    """
    f = SeqFeature(FeatureLocation(location.start, location.end,
                   strand=location.strand), type=genType)
    seq = f.extract(contig.seq)
    if translate:
        return seq.translate(table=11, cds=False, to_stop=True)
    else:
        return seq


def create_genome_object(gff_file=None, fna_file=None, gff_gz_file=None):
    """Generates a genome object given a GFF and Fasta file.

        Arguments:
            gff_file: Path to the GFF annotation file
            fna_file: Path to the genome fasta file
            gff_gz_file: Path to combined gzipped GFF+Fasta file

        Returns:
            genomeObj: ProkFunFind.annotate.Genome object
    """
    if gff_gz_file:
        genomeObj = GetGenomeFromGzipGFF(gzipfile=gff_gz_file)
    else:
        genomeObj = GetGenomeFromGFF(gff3file=gff_file, fnafile=fna_file)
    [ct.sort() for ct in genomeObj.contigs.values()]
    return genomeObj


def parse_gtab(genome_tab: str):
    """Parses a genome input table file

        Arguments:
            genome_tab: Path to the genome input table

        Returns:
            gtab: A dictionary linking genome IDs to names
    """
    # Add in check to make sure that file exists.
    gtab = {}
    for row in csv.reader(open(genome_tab), delimiter="\t"):
        gtab[row[0]] = row[1]

    return gtab
