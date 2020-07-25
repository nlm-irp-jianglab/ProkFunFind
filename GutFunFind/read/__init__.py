from typing import Dict, IO, List, Set, OrderedDict
from BCBio import GFF
from Bio import SeqIO

from Bio.SeqFeature import FeatureLocation,SeqFeature
from Bio.Seq import Seq

class Gene(SeqFeature):
    #def __init__(self,contig:str,id:str,type:str,location:FeatureLocation, qualifiers:OrderedDict=None) -> None:
    #    super().__init__(id=id, type=type, location=location, qualifiers=qualifiers)
    def __init__(self,contig:str,id:str, **kwds) -> None:
        super().__init__(id=id, **kwds)
        self.contig= contig
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
    def __init__(self,id:str,seq:Seq,genes:OrderedDict[str,Gene],sorted = False) -> None:
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

    def __getitem__(self, key:str): 
        return self.genes[key]
    
    def __setitem__(self, key:str, gene:Gene) -> None:
        self.genes.update({key:gene})

    def __delitem__(self, key:str) -> None:
        try:
            del self.genes[key]
        except KeyError:
            print("Gene:{} not found".format(KeyError))

    def __contains__(self, item):
        pass

    def __contains__(self, key:str):
        return key in self.genes.keys()

    def sort(self):
        self.genes = OrderedDict(sorted(self.genes.items(), key=lambda x: (x[1].location.start.real,x[1].location.end.real)))
        #self.genes = {k: v for k, v in sorted(self.genes.items(), key=lambda x: (x[1].location.start.real,x[1].location.end.real))}
        self.sorted = True
        return self

class Genome:
    """Docstring for Genome. """
    def __init__(self, contigs:OrderedDict[str,Contig]) -> None:
        self.contigs = contigs

    def __repr__(self):
        return "{0}(#contigs:{1})".format( 
                 self.__class__.__name__, str(len(self.contigs)))

    @property
    def genes(self):
        return dict(j for i in self.contigs.keys() for j in self.contigs[i].genes.items())

    def __getitem__(self, key:str): 
        return self.contigs[key]

    def __setitem__(self, key:str, contig:Contig):
        self.contigs.update({key:contig})

    def __delitem__(self, key:str) -> None:
        try:
            del self.contigs[key]
        except KeyError:
            print("Contig:{} not found".format(KeyError))

    def __contains__(self, key:str):
        return key in self.contigs.keys()


def GetGenomeFromGFF(gff3file:str,fnafile:str) -> Genome:

    in_seq_handle = open(fnafile)
    g1 = Genome(contigs={i.id:Contig(id=i.id,seq=i.seq,genes={}) for i in SeqIO.parse(in_seq_handle,"fasta")})
    in_seq_handle.close()

    try:
        with open(gff3file) as handle:
            gff_records = list(GFF.parse(handle))
        for rec in gff_records:
            ct = g1[rec.id]
            for f in rec.features:
                if f.type == "CDS":
                    gene = Gene(contig=ct.id,
                        id = f.id,
                        type = f.type,
                        qualifiers=f.qualifiers,
                        location = f.location)
                    ct[gene.id] = gene
                elif f.type in ["assembly_gap", "source"]:
                    pass
                else:
                    gene = Gene(contig=ct.id,
                        id = f.id,
                        type = f.type,
                        qualifiers=f.qualifiers, 
                        location = f.location)
                    ct[gene.id] = gene
    except Exception as err:
        print("Error:{}".format(err))
    return g1


def GetGenomeFromGB(gbfile:str) -> Genome:
    """

    :gbfile: the path to genome file (formated in genbank)
    :returns: A Genome object

    """
    g1 = Genome(contigs={})
    for rec in SeqIO.parse(gbfile,"genbank"):
        ct = Contig(id = rec.id, seq=rec.seq,genes={})
        for f in rec.features:
            if f.type == "CDS":
                fid = f.qualifiers["locus_tag"][0]
                gene = Gene(contig=rec.id,
                    id = fid,
                    type = f.type,
                    qualifiers=f.qualifiers,
                    location = f.location)
                ct[gene.id] = gene
            elif f.type in ["assembly_gap", "source"]:
                pass
            else:
                fid = f.qualifiers["locus_tag"][0]
                gene = Gene(contig=rec.id,
                    id = fid,
                    qualifiers=f.qualifiers,
                    type = f.type,
                    location = f.location)
                ct[gene.id] = gene
        g1[ct.id] = ct
    return g1
