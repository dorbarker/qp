from collections import defaultdict
from qp import calculate_distance as d
import networkx as nx
import numpy as np

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

    def find_closest(self, cluster, gene, entry_point, done=None):

        done = done or []
        possible_nexts = {}
        d = partial(calculate_distance, gene1=gene.sequence)

        neighbours = set(nx.all_neighbors(cluster, entry_point)) - set(done)

        if len(neighbours) is not 0:
            for node in neighbours:
                if node not in done:

                    dist = d(gene.sequence, node.sequence)

                    gene.compared[node] = dist
                    node.compared[gene] = dist

                    possible_nexts[node] = dist

                best = min(possible_nexts,
                           key=lambda x: possible_nexts[x].compared)
                find_closest(cluster, gene, best, done)

        else:
            self.add_to_cluster(cluster, gene, node, d(node.sequence))

class GeneNode(object):

    def __init__(self, genome, gene, sequence):

        genome = genome
        gene = gene
        sequence = sequence
        compared = {}
