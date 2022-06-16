from typing import List
import numpy as np
from sklearn.cluster import DBSCAN

from GutFunFind.read import Genome
# from GutFunFind.toolkit.utility import *


def DBSCANCluster(
        pos: List[int],
        weight: List[float] = None,
        eps: float = 2,
        min_samples: float = 0.9) -> List[int]:
    pos = np.array(pos).reshape(-1, 1)
    cluster = DBSCAN(
        eps=eps, min_samples=min_samples).fit(pos, sample_weight=weight)
    return list(cluster.labels_)


def pipeline(config: dict, genome_object, detect_tool: str) -> Genome:

    # cf = read_config(config_file)["DBSCAN"]

    for ct in genome_object.contigs.values():

        # obtain weight for genes with blast attributes and the postion in the
        # contig(0-based)
        poslist, weightlist, genelist = [], [], []
        for i, gene in enumerate(ct.genes.values()):
            if hasattr(gene, detect_tool):
                poslist.append(i)
                genelist.append(gene)
                weightlist.append(
                    getattr(
                        getattr(
                            gene,
                            detect_tool),
                        "orthoID_weight"))
                # weightlist.append(gene.blast.orthoID_weight)

        # if there is a hit in the contig, run DBSCANCluster
        if len(poslist) >= 1:
            cl = DBSCANCluster(pos=poslist,
                               weight=weightlist,
                               eps=float(config['DBSCAN']['cluster.eps']),
                               min_samples=float(config['DBSCAN'][
                                                 'cluster.min_samples']))

            for gene, cl_label in zip(genelist, cl):
                if cl_label == -1:
                    cl_label = "NA"
                setattr(gene, "DBSCAN", cl_label)

    return genome_object
