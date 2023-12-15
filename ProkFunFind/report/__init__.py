import yaml
from typing import Any, Dict

from ProkFunFind.read import Genome


def export_pickle(genomeObj: Genome, outprefix: str) -> None:
    import pickle
    handle = open(outprefix + ".pkl", "wb")
    pickle.dump(genomeObj, file=handle)
    handle.close()


def export_yaml(variable: Dict[str, Any], outprefix: str) -> None:
    """TODO: Docstring for export_yaml.
    :variable: TODO
    :returns: TODO
    """
    with open(outprefix + ".yaml", "w") as outfile:
        yaml.dump(variable, outfile, sort_keys=False)


def export_gene_tab(
        genomeObj: Genome,
        detect_tools: list,
        cluster_tool: str,
        outprefix: str) -> None:
    """TODO: Docstring for export_gene_tab.

    """
    f = open(outprefix + ".tsv", "w")
    f.write("Gene_Name\tCluster_ID\tFunctions\n")
    f_set = set()
    for gene in genomeObj.genes.values():
        for detect_tool in detect_tools:
            if hasattr(gene, detect_tool):
                cluster_annot = gene.contig + ":Cl_" + \
                    str(getattr(gene, cluster_tool)) \
                    if hasattr(gene, cluster_tool) else "NA"
                function_annot = ";".join(set(gene.Functions)) \
                    if hasattr(gene, "Functions") else "Unassigned Function"
                s = gene.id + "\t" + cluster_annot + "\t" + function_annot+"\n"
                f_set.add(s)
    for i in sorted(f_set):
        f.write(i)
    f.close()


def export_gene_gff(
        genomeObj: Genome,
        detect_tools: list,
        cluster_tool: str,
        outprefix: str) -> None:
    """TODO: Docstring for export_gene_tab.

    """
    f = open(outprefix + ".annot.gff", "w")
    for gene in genomeObj.genes.values():
        for detect_tool in detect_tools:
            if hasattr(gene, detect_tool):
                qry = getattr(gene, detect_tool)
                hsp = qry.hsps[0]
                for hsp_t in qry.hsps:
                    if hsp_t.evalue < hsp.evalue:
                        hsp = hsp_t
                    else:
                        continue
                hit = qry.hits[0]
                cluster_annot = "Cl_" + \
                    str(getattr(gene, cluster_tool)) if hasattr(
                        gene, cluster_tool) else "NA"

                if detect_tool == "blast":
                    f.write("{ct}\ProkFunFind\t{tp}\t{start}\t{end}" \
                        "\t.\t{strand}\t.\tID={id};Name={geneID};" \
                        "ClusterID={cluster_ID};Target={Target};" \
                        "pct_identity={pct_identity};evalue={evalue}".format(
                        ct=gene.contig,
                        tp=gene.type,
                        start=gene.location.start+1,
                        end=gene.location.end,
                        strand="+" if gene.strand == 1 else "-",
                        id=gene.id,
                        geneID=qry.geneID,
                        cluster_ID=cluster_annot,
                        Target=hsp.hit_id + " " +
                        str(hsp.hit_start) + " " + str(hsp.hit_end),
                        pct_identity=hsp.ident_pct,
                        evalue=hsp.evalue
                        ))

                    if hasattr(gene, "pangenome_group"):
                        gene_group = gene.pangenome_group
                        if gene_group.method == "Roary":
                            f.write(";pan_annot=" + gene_group.annotation)
                            f.write(";pan_isolates=" +
                                    gene_group.number['isolates'])
                            f.write(";pan_avg_seq=" +
                                    gene_group.number['Avg sequences'])
                            if hasattr(gene_group, "accessory"):
                                f.write(";pan_accessory=True")
                            else:
                                f.write(";pan_core=True")

                    f.write("\n")

                elif detect_tool == "interproscan" \
                        or detect_tool == "hmmer" \
                        or detect_tool == "kofamscan" \
                        or detect_tool == "emapper":
                    if detect_tool == 'emapper' or detect_tool == 'interproscan':
                        if hasattr(qry, 'evalue'):
                            eval = ';evalue='+str(qry.evalue)
                        else:
                            eval = ''
                    else:
                        if hasattr(hsp, 'evalue'):
                            eval = ';evalue='+str(hit.evalue)
                        else:
                            eval = ''
                    f.write("{ct}\tProkFunFind\t{tp}\t{start}\t{end}\t." \
                            "\t{strand}\t.\tID={id};Name={geneID};" \
                            "ClusterID={cluster_ID};" \
                            "Target={Target}{evalue}".format(
                        ct=gene.contig,
                        tp=gene.type,
                        start=gene.location.start+1,
                        end=gene.location.end,
                        strand="+" if gene.strand == 1 else "-",
                        id=gene.id,
                        geneID=qry.geneID,
                        cluster_ID=cluster_annot,
                        Target=hsp.hit_id,
                        evalue=eval))

                    if hasattr(gene, "pangenome_group"):
                        gene_group = gene.pangenome_group
                        if gene_group.method == "Roary":
                            f.write(";pan_annot=" + gene_group.annotation)
                            f.write(";pan_isolates=" +
                                    gene_group.number['isolates'])
                            f.write(";pan_avg_seq=" +
                                    gene_group.number['Avg sequences'])
                            f.write(";pan_accessory=True") \
                                if hasattr(gene_group, "accessory") \
                                else f.write("pan_core=True")

                    f.write("\n")

    f.close()


def report_all(system_dict: Dict[str, Any],
               status: int,
               genomeObj: Genome,
               outprefix: str,
               detect_tools: list,
               cluster_tool: str) -> None:

    export_pickle(genomeObj, outprefix)
    export_yaml(system_dict, outprefix)
    export_gene_tab(genomeObj, detect_tools, cluster_tool, outprefix)
    export_gene_gff(genomeObj, detect_tools, cluster_tool, outprefix)
