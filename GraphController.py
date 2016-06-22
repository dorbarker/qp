from collections import defaultdict
from functools import partial
import networkx as nx
import numpy as np
import utilities

class Pangenome(object):

    def __init__(self):

       self.clusters = defaultdict(nx.Graph)

    def add_to_cluster(self, cluster: str, gene1: 'GeneNode', gene2: 'GeneNode',
                       distance: float):

        self.clusters[cluster].add_edge(gene1, gene2, weight=distance)

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
        d = partial(utilities.calculate_distance, gene1=gene.sequence)

        neighbours = set(nx.all_neighbors(cluster, entry_point)) - set(done)

        if len(neighbours):
            for node in neighbours:
                if node in done:
                    continue
                #dist = d(gene.sequence, node.sequence)
                #dist = d(gene2=node.sequence)
                #gene.compared[node] = dist
                #node.compared[gene] = dist

                gene.update_compared(node)

                possible_nexts[node] = gene.compared[node]

                done.append(node)

            best = best_next(possible_nexts)

            if entry_point.compared[gene] < possible_nexts[best]:

                self.add_to_cluster(cluster, gene, entry_point,
                                    entry_point.compared[gene])

            else:
                self.find_closest(cluster, gene, best, done)

        else:
            self.add_to_cluster(cluster, gene, entry_point,
                                entry_point.compared[gene])

class GeneNode(object):

    def __init__(self, genome, gene, sequence):

        self.genome = genome
        self.gene = gene
        self.sequence = sequence
        self.compared = {}

    def update_compared(self, other):

        dist = utilities.calculate_distance(self, other)

        self.compared[other] = dist
        other.compared[self] = dist

    def __repr__(self):
        return self.gene
