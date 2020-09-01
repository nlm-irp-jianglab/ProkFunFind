import gzip
import csv
from collections import OrderedDict
from typing import MutableMapping, List

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq

from GutFunFind.read.GFFParser import parse

####################################################################################################

class Gene(SeqFeature):
    def __init__(self, contig: str, id: str, **kwds) -> None:
        super().__init__(id=id, **kwds)
        self.contig = contig
        #self.id = id
        #self.type = type
        #self.nt = nt
        #self.aa = aa
        #self.location = location

    def __repr__(self):
        answer = "%s(%s" % (self.__class__.__name__, repr(self.location))
        if self.id and self.id != "<unknown id>":
            answer += ", id=%s" % repr(self.id)
            answer += ")"
        return answer


class Contig:
    def __init__(self,
                 id: str,
                 seq: Seq,
                 genes: MutableMapping[str, Gene],
                 sorted=False) -> None:
        self.id = id
        self.seq = seq
        self.genes = genes
        self.sorted = sorted

    def __repr__(self):
        return "{0}(seq={1!r}, id={2!r}, sorted={3}, #genes:{4})".format(
            self.__class__.__name__,
            self.seq, self.id, str(self.sorted), str(len(self.genes)))

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key: str):
        return self.genes[key]

    def __setitem__(self, key: str, gene: Gene) -> None:
        self.genes.update({key: gene})

    def __delitem__(self, key: str) -> None:
        try:
            del self.genes[key]
        except KeyError:
            print("Gene:{} not found".format(KeyError))

    def __iter__(self):
        return iter(self.genes)

    def __contains__(self, key: str):
        return key in self.genes.keys()

    def sort(self):
        self.genes = OrderedDict(
            sorted(
                self.genes.items(),
                key=lambda x: (
                    x[1].location.start.real,
                    x[1].location.end.real)))
        self.sorted = True
        return self


class Genome:
    """Docstring for Genome. """

    def __init__(self, contigs: MutableMapping[str, Contig]) -> None:
        self.contigs = contigs

    def __repr__(self):
        return "{0}(#contigs:{1})".format(
            self.__class__.__name__, str(len(self.contigs)))

    @property
    def genes(self):
        return dict(j for i in self.contigs.keys()
                    for j in self.contigs[i].genes.items())

    def __getitem__(self, key: str):
        return self.contigs[key]

    def __setitem__(self, key: str, contig: Contig):
        self.contigs.update({key: contig})

    def __delitem__(self, key: str) -> None:
        try:
            del self.contigs[key]
        except KeyError:
            print("Contig:{} not found".format(KeyError))

    def __iter__(self):
        return iter(self.contigs)

    def __contains__(self, key: str):
        return key in self.contigs.keys()

####################################################################################################


def GetGenomeFromGFF(gff3file: str, fnafile: str) -> Genome:

    in_seq_handle = open(fnafile)
    g1 = Genome(
        contigs={
            i.id: Contig(
                id=i.id,
                seq=i.seq,
                genes={}) for i in SeqIO.parse(in_seq_handle, "fasta")})
    in_seq_handle.close()

    try:
        with open(gff3file) as handle:
            gff_records = list(parse(handle))
            #gff_records = list(GFF.parse(handle))
        for rec in gff_records:
            ct = g1[rec.id]
            for f in rec.features:
                if f.type == "CDS":
                    gene = Gene(contig=ct.id,
                                id=f.id,
                                type=f.type,
                                qualifiers=f.qualifiers,
                                location=f.location)
                    ct[gene.id] = gene
                elif f.type in ["assembly_gap", "source"]:
                    pass
                else:
                    gene = Gene(contig=ct.id,
                                id=f.id,
                                type=f.type,
                                qualifiers=f.qualifiers,
                                location=f.location)
                    ct[gene.id] = gene
    except Exception as err:
        print("Error:{}".format(err))
    return g1


def GetGenomeFromGB(gbfile: str) -> Genome:
    """

    :gbfile: the path to genome file (formated in genbank)
    :returns: A Genome object

    """
    g1 = Genome(contigs={})
    for rec in SeqIO.parse(gbfile, "genbank"):
        ct = Contig(id=rec.id, seq=rec.seq, genes={})
        for f in rec.features:
            if f.type == "CDS":
                fid = f.qualifiers["locus_tag"][0]
                gene = Gene(contig=rec.id,
                            id=fid,
                            type=f.type,
                            qualifiers=f.qualifiers,
                            location=f.location)
                ct[gene.id] = gene
            elif f.type in ["assembly_gap", "source"]:
                pass
            else:
                fid = f.qualifiers["locus_tag"][0]
                gene = Gene(contig=rec.id,
                            id=fid,
                            qualifiers=f.qualifiers,
                            type=f.type,
                            location=f.location)
                ct[gene.id] = gene
        g1[ct.id] = ct
    return g1


def GetGenomeFromGzipGFF(gzipfile: str) -> Genome:

    in_seq_handle = gzip.open(gzipfile, "rt")
    g1 = Genome(
        contigs={
            i.id: Contig(
                id=i.id,
                seq=i.seq,
                genes={}) for i in SeqIO.parse(in_seq_handle, "fasta")})
    in_seq_handle.close()

    try:
        in_seq_handle = gzip.open(gzipfile, "rt")
        gff_records = list(parse(in_seq_handle))
        for rec in gff_records:
            ct = g1[rec.id]
            for f in rec.features:
                if f.type == "CDS":
                    gene = Gene(contig=ct.id,
                                id=f.id,
                                type=f.type,
                                qualifiers=f.qualifiers,
                                location=f.location)
                    ct[gene.id] = gene
                elif f.type in ["assembly_gap", "source"]:
                    pass
                else:
                    gene = Gene(contig=ct.id,
                                id=f.id,
                                type=f.type,
                                qualifiers=f.qualifiers,
                                location=f.location)
                    ct[gene.id] = gene
    except Exception as err:
        print("Error:{}".format(err))
    in_seq_handle.close()

    return g1


####################################################################################################

class PanGene:
    def __init__(self, id: str, GenomeID: str, GeneGroup: str) -> None:
        self.id = id
        self.GenomeID = GenomeID
        self.GeneGroup = GeneGroup

    def __repr__(self):
        return "{0}(id={1!r}, GenomeID:{2!r}, GeneGroup:{3!r})".format(
            self.__class__.__name__,
            self.id,
            self.GenomeID,
            self.GeneGroup
        )


class GeneGroup:
    def __init__(self,
                 id: str,
                 pangenes: MutableMapping[str, List[str]]
                 ) -> None:
        self.id = id
        self.pangenes = pangenes

    def __repr__(self):
        return "{0}(id={1!r}, #isolates={2!r}, #sequences={3!r})".format(
            self.__class__.__name__,
            self.id,
            len(self.pangenes),
            sum([len(self.pangenes[i]) for i in self.pangenes])
        )

    def __getitem__(self, key: str) -> List[str]:
        return self.pangenes[key]

    def __setitem__(self, key: str, genelist: List) -> None:
        self.pangenes.update({key: genelist})

    def __delitem__(self, key: str) -> None:
        try:
            del self.pangenes[key]
        except KeyError:
            print("Gene:{} not found".format(KeyError))

    def __iter__(self):
        return iter(self.pangenes)

    def __contains__(self, key: str):
        return key in self.pangenes.keys()


class PanGenome:
    """Docstring for Genome. """

    def __init__(self, gene_groups: MutableMapping[str, GeneGroup] = OrderedDict(), genomes: List = []) -> None:
        self.gene_groups = gene_groups
        self.genomes = genomes

    def __repr__(self):
        return "{0}(#gene_groups={1} #genomes={2})".format(
            self.__class__.__name__,
            len(self.gene_groups),
            len(self.genomes)
        )

    @property
    def pangenes(self):
        return {gene: PanGene(id=gene, GenomeID=key, GeneGroup=gene_group)
                for gene_group in self.gene_groups.values()
                for key in gene_group
                for gene in gene_group[key]
                }

    def __getitem__(self, key: str):
        return self.gene_groups[key]

    def __setitem__(self, key: str, value: GeneGroup):
        self.gene_groups.update({key: value})

    def __delitem__(self, key: str) -> None:
        try:
            del self.gene_groups[key]
        except KeyError:
            print("GeneGroup:{} not found".format(KeyError))

    def __contains__(self, key: str):
        return key in self.gene_groups.keys()

    def __iter__(self):
        return iter(self.gene_groups)

    def get_genegroup(self, genename):
        genegroup = self.pangenes[genename].GeneGroup
        return genegroup

####################################################################################################


def Roarycsv2pangenome(in_file):
    with open(in_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        headers = next(reader, None)
        pg = PanGenome(genomes=headers[14:])
        ncol = len(headers)
        for row in reader:
            geneG = GeneGroup(
                id=row[0],
                pangenes={headers[idx]: row[idx].split(
                    "\t") for idx in range(14, ncol, 1)}
            )
            # set Roary-specific annotation attributes
            setattr(geneG, "method", "Roary")
            setattr(geneG, "common_name", row[1])
            setattr(geneG, "annotation", row[2])
            setattr(geneG, "number", {"isolates": row[3], "sequences": row[4], "Avg sequences": row[5]})
            setattr(geneG, "overall", {"Fragment": row[6], "Order": row[7]})
            if row[8]:
                setattr(geneG, "accessory", {"Fragment": row[8], "Order": row[9]})
            if row[10]:
                setattr(geneG, "QC", row[10])
            setattr(geneG, "group_size", {"min": row[11], "max": row[12], "avg": row[13]})

            pg[row[0]] = geneG
    return pg
