import networkx as nx
import numpy as np
from collections import defaultdict

class Pangenome(object):

    def __init__(self):

       clusters = defaultdict(nx.Graph)

    def add_to_cluster(self, cluster: str, gene1: GeneNode, gene2: GeneNode,
                       distance: float):

        self.clusters[cluster].add_edge(gene1.gene, gene2.gene,
                                        weight=distance)

    def join_clusters(self, cluster1: str, cluster2: str,
                      gene1: Seq, gene2: Seq, middle_gene: str,
                      distance1: float, distance2: float):

        G = self.clusters.pop(cluster2)
        self.clusters[cluster1].add_edges_from(self.clusters[cluster2].edges())
        self.add_to_cluster(cluster1, gene1, middle_gene, distance1)
        self.add_to_cluster(cluster1, middle_gene, gene2, distance2)

    def add_founder(self, gene):
        """Creates the founding node of a new cluster graph.

        The graph is named after the founding gene, and is *not*
        updated when centre changes.
        """

        self.clusters[gene].add_node(gene.id, sequence=gene.seq)

class GeneNode(object):

    def __init__(self, genome, gene, sequence):

        genome = genome
        gene = gene
        sequence = sequence
