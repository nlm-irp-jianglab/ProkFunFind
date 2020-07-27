import json
from GutFunFind.read import *
from typing import Any, Dict


def export_pickle(genomeObj: Genome, outprefix: str) -> None:
    import pickle
    handle = open(outprefix + ".pkl", "wb")
    pickle.dump(genomeObj, file=handle)
    handle.close()


def export_json(variable: Dict[str, Any], outprefix: str) -> None:
    """TODO: Docstring for export_json.
    :variable: TODO
    :returns: TODO
    """
    with open(outprefix + ".json", "w") as outfile:
        json.dump(variable, outfile, indent=4)


def export_gene_tab(
        genomeObj: Genome,
        detect_tool: str,
        cluster_tool: str,
        outprefix: str) -> None:
    """TODO: Docstring for export_gene_tab.

    """
    f = open(outprefix + ".tsv", "a")
    print("Gene_Name\tCluster_ID\tFunctions", file=f)
    for gene in genomeObj.genes.values():
        if hasattr(gene, detect_tool):
            cluster_annot = gene.contig + ":Cl_" + \
                str(getattr(gene, cluster_tool)) if hasattr(gene, cluster_tool) else "NA"
            function_annot = ";".join(
                gene.Functions) if hasattr(
                gene, "Functions") else "Unassigned Function"
            print(
                gene.id +
                "\t" +
                cluster_annot +
                "\t" +
                function_annot,
                file=f)
    f.close()


def export_gene_gff(
        genomeObj: Genome,
        detect_tool: str,
        cluster_tool: str,
        outprefix: str) -> None:
    """TODO: Docstring for export_gene_tab.

    """
    f = open(outprefix + ".annot.gff", "a")
    for gene in genomeObj.genes.values():
        if hasattr(gene, detect_tool) and detect_tool == "blast":
            qry = getattr(gene, detect_tool)
            hsp = qry.hsps[0]
            cluster_annot = "Cl_" + \
                str(getattr(gene, cluster_tool)) if hasattr(gene, cluster_tool) else "NA"
            #function_annot = ";".join(gene.Functions) if hasattr(gene,"Functions") else "Unassigned Function"
            print("{ct}\tGuFunFind\t{tp}\t{start}\t{end}\t.\t{strand}\t.\tID={id};Name={orthoID};Parent={cluster_ID};Target={Target};pct_identity={pct_identity};evalue={evalue}".format(
                ct=gene.contig,
                tp=gene.type,
                start=gene.location.start,
                end=gene.location.end,
                strand="+" if gene.strand == 1 else "-",
                id=gene.id,
                orthoID=qry.orthoID,
                cluster_ID=cluster_annot,
                Target=hsp.hit_id + " " + str(hsp.hit_start) + " " + str(hsp.hit_end),
                pct_identity=hsp.ident_pct,
                evalue=hsp.evalue
            ), file=f)
    f.close()


def report_all(system_dict: Dict[str,
                                 Any],
               status: int,
               genomeObj: Genome,
               outprefix: str,
               detect_tool: str,
               cluster_tool: str) -> None:

    export_pickle(genomeObj, outprefix)
    export_json(system_dict, outprefix)
    export_gene_tab(genomeObj, detect_tool, cluster_tool, outprefix)
    export_gene_gff(genomeObj, detect_tool, cluster_tool, outprefix)
