import csv
import subprocess

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from GutFunFind.read import GetGenomeFromGFF, GetGenomeFromGzipGFF

def run_prokka(config: dict, genome_dir: str, genome_tab: str):
    """

        Arguments:
            genome_tab

        Returns:
            None
    """
    gids = parse_gtab(genome_tab)
    for genome in gids.keys():
        cmd = ['prokka', '--outdir', genome_dir, '--locustag', gids[genome]+'_', genome_dir+'/'+genome]
        res = subprocess.run(cmd)
        if res.returncode != 0:
            raise RuntimeError("Failed to run: {}".format(" ".join(cmd)))

def export_proteins(config: dict, genome_dir: str, genome_tab: str, zipped: False):
    """ Generates protein and nucletotide fasta files from GFF and Genome file

        Arguments:
            config: configuration dictionary object
            genome_dir: Path to input directory of genomes to search
            genome_tab: Path to table of input file prefixes and IDs
            zipped: True or False value indicating if the input is gzipped

        Returns:
            A list of input input paths with genome prefixes the input genomes
    """
    # gff_path = './test-genomein/GTDB15451.comb.gff3'
    # fna_path = None
    # zipped
    prefixes = []
    gids = parse_gtab(genome_tab)
    for genome in gids.keys():
        # fna_path = None
        prefix = genome_dir+'/'+genome
        gff_path = prefix+'.gff'
        # if fna_path:
        #     g = create_genome_object(gff_file=gff_path, fna_file=fna_path)
        # else:
        #     if zipped:
        #         g = create_genome_object(gff_gz_file=gff_path)
        #     else:
        extract_fna(gff_path, prefix)
        g = create_genome_object(gff_file=prefix+'.pff.gff3', fna_file=prefix+'.pff.fna')

        g1 = open(prefix+".pff.ffn", "w")
        g2 = open(prefix+".pff.faa", "w")
        prefixes.append(prefix)


        for contig in g.contigs.values():
            for gene in contig.genes.values():
                if gene.type != "CDS":
                    continue
                nt_seq = extract_gene(contig, gene.location, gene.type, translate=False)
                aa_seq = extract_gene(contig, gene.location, gene.type, translate=True)
                g1.write(">"+gene.id+"\n")
                g1.write(str(nt_seq)+"\n")
                g2.write(">"+gene.id+"\n")
                g2.write(str(aa_seq)+"\n")

    return prefixes

def extract_fna(comb_gff_path: str, prefix):
    """ Generates separate genome fasta file and gff file from combined file

        Arguments:
            comb_gff_path: path to combined fasta+gff file
            prefix: Genome prefix for output file names

        Returns:
            None
    """
    with open(comb_gff_path) as f:
        txt = f.read()
        try:
            gff_str, fna_str = txt.split("##FASTA\n")
        except ValueError as e:
            print("Error: No fasta tag '##FASTA' found in gff file")
    with open(prefix+'.pff.gff3', 'w') as g:
        g.write(gff_str)
    with open(prefix+'.pff.fna', 'w') as g:
        g.write(fna_str)

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
    f = SeqFeature(FeatureLocation(location.start, location.end, strand=location.strand), type=genType)
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
    for row in csv.reader(open(genome_tab), delimiter='\t'):
        gtab[row[0]] = row[1]

    return gtab
