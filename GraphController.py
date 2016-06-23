from collections import defaultdict
from functools import partial
import networkx as nx
import numpy as np
import utilities

class Pangenome(object):

    def __init__(self):

       self.clusters = defaultdict(nx.Graph)

    def add_to_cluster(self, cluster: str,
                       gene1: 'GeneNode', gene2: 'GeneNode'):

        self.clusters[cluster].add_edge(gene1, gene2,
                                        weight=gene1.compared[gene2])

    def join_clusters(self, cluster1: str, cluster2: str,
                      gene1: 'Seq', gene2: 'Seq', middle_gene: str,
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

        self.clusters[gene].add_node(gene)

    def find_closest(self, cluster, gene, entry_point, done=None):

        def best_next(nexts):

            return min(nexts, key=lambda x: nexts[x])

        done = done or [entry_point]
        possible_nexts = {}

        neighbours = set(nx.all_neighbors(cluster, entry_point)) - set(done)

        for n in (node for node in neighbours if node not in done):

            gene.update_compared(n)

            possible_nexts[n] = gene.compared[n]

            done.append(n)

        best = best_next(possible_nexts) if possible_nexts else entry_point

        if not len(neighbours) or entry_point.compared[gene] < possible_nexts[best]:

            closest = entry_point

        else:

            closest = self.find_closest(cluster, gene, best, done)

        return closest

class GeneNode(object):

    def __init__(self, genome, gene, sequence):

        self.genome = genome
        self.gene = gene
        self.sequence = sequence
        self.compared = {}

    def update_compared(self, other):

        if other not in self.compared:

            dist = utilities.calculate_distance(self, other)

            self.compared[other] = dist
            other.compared[self] = dist

            return dist

        else:

            return self.compared[other]

    def __repr__(self):
        return self.gene
