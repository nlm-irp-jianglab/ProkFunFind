from typing import List
import numpy as np
from sklearn.cluster import DBSCAN

from ProkFunFind.read import Genome
# from ProkFunFind.toolkit.utility import *


def DBSCANCluster(
        pos: List[int],
        weight: List[float] = None,
        eps: int = 2,
        min_samples: int = 1) -> List[int]:
    """Wrapper function to run DBSCAN clustering

       Arguments:
           pos:
           weight:
           eps:
           min_samples:

       Returns:
           list of cluster labels for genes
    """
    pos = np.array(pos).reshape(-1, 1)
    cluster = DBSCAN(
        eps=eps, min_samples=min_samples).fit(pos, sample_weight=weight)
    return list(cluster.labels_)


def pipeline(config: dict, genome_object, detect_tools: set) -> Genome:
    """Main pipeline function for DBSCAN clustering analysis

       Arguments:
           config: A config dictionary
           genome_object: A ProkFunFind.annotate.Genome object
           detect_tools: list of gene hit detection tools used
       Returns:
           genome_object: ProkFunFind.annotate.Genome object with cluster info
    """
    for ct in genome_object.contigs.values():

        # obtain weight for genes with blast attributes and the postion in the
        # contig(0-based)
        poslist, weightlist, genelist = [], [], []
        for detect_tool in detect_tools:
            for i, gene in enumerate(ct.genes.values()):
                if hasattr(gene, detect_tool):
                    poslist.append(i)
                    genelist.append(gene)
                    weightlist.append(
                        getattr(
                            getattr(
                                gene,
                                detect_tool),
                            "queryID_weight"))

        # if there is a hit in the contig, run DBSCANCluster
        if len(poslist) >= 1:
            cl = DBSCANCluster(pos=poslist,
                               weight=weightlist,
                               eps=float(config['DBSCAN']['cluster_eps']),
                               min_samples=int(config['DBSCAN'][
                                                 'cluster_min_samples']))

            for gene, cl_label in zip(genelist, cl):
                if cl_label == -1:
                    cl_label = "NA"
                setattr(gene, "DBSCAN", cl_label)

    return genome_object
