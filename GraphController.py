from collections import defaultdict
from functools import partial
from multiprocessing import cpu_count
import networkx as nx
import numpy as np
import utilities
from concurrent.futures import ThreadPoolExecutor

class Pangenome(object):

    def __init__(self):

       self.clusters = defaultdict(nx.Graph)

    def add_to_cluster(self, cluster: str,
                       gene1: 'GeneNode', gene2: 'GeneNode'):

        self.clusters[cluster].add_edge(gene1, gene2,
                                        weight=gene1.compared[gene2])

    def import_nodes(self, cluster1, cluster2):

        if len(self.clusters[cluster2]) is 1:

            self.clusters[cluster1].add_nodes_from(self.clusters[cluster2])

        else:

            self.clusters[cluster1].add_edges_from(self.clusters[cluster2])

        del self.clusters[cluster2]

    def join_clusters(self, cluster1: str, cluster2: str, gene: 'GeneNode'):


        c1_closest = self.find_closest(cluster1.cluster, gene, cluster1.entry)
        c2_closest = self.find_closest(cluster2.cluster, gene, cluster2.entry)

        self.import_nodes(cluster1.cluster, cluster2.cluster)

        self.add_to_cluster(cluster1.cluster, gene, c1_closest)
        self.add_to_cluster(cluster1.cluster, gene, c2_closest)

    def add_founder(self, gene):
        """Creates the founding node of a new cluster graph.

        The graph is named after the founding gene, and is *not*
        updated when centre changes.
        """

        self.clusters[gene].add_node(gene)

    def find_closest(self, cluster, gene, entry_point, cpus, done=None):

        def best_next(nexts):
            try:
                return min(nexts, key=lambda x: nexts[x])[0]

            except ValueError:
                return entry_point

        def explore_neighbours(node, cpus):

            gene.update_compared(node, cpus)

            return node, gene.compared[node]


        done = done or set([entry_point])

        all_neighbours = nx.all_neighbors(self.clusters[cluster], entry_point)
        neighbours =  set(all_neighbours) - set(done)
        n_neighbours = len(neighbours)

        threads = utilities.calculate_search_threads(n_neighbours, cpus)

        with ThreadPoolExecutor(threads) as executor:

            f = partial(explore_neighbours, cpus=threads)
            search = executor.map(f, neighbours)

        done = done | neighbours
        best = best_next(search)

        if not len(neighbours) or \
           entry_point.compared[gene] < possible_nexts[best]:

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

    def update_compared(self, other, cpus=cpu_count()):

        if other not in self.compared:

            dist = utilities.calculate_distance(self, other, cpus)

            self.compared[other] = dist
            other.compared[self] = dist

            return dist

        else:

            return self.compared[other]

    def __repr__(self):
        return self.gene
