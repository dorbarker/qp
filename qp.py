#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import json
import os
import GraphController as GC
import networkx as nx
import utilities
import matplotlib.pyplot as plt

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', nargs='+')

    return parser.parse_args()


def annotations(ffn):

    def genome_basename(n):
        return os.path.splitext(os.path.basename(n))[0]

    with open(ffn) as f:
        for rec in SeqIO.parse(f, 'fasta'):

            yield GC.GeneNode(genome_basename(ffn),
                              rec.id, rec.seq)

def add_to_graph(gene, pangenome, threshold):


    def recursive_cluster_join(c1, clusters):

        if clusters:

            c2, *clusters  = clusters

            pangenome.join_clusters(c1, c2, gene)

            recursive_cluster_join(c1, clusters)

    if pangenome.clusters:

        matching_clusters = []

        for cluster in pangenome.clusters:

            centre = min(nx.center(pangenome.clusters[cluster]),
                         key=gene.update_compared)

            if gene.compared[centre] < threshold:

                cluster_entry = utilities.Cluster_Entry(cluster, centre)

                matching_clusters.append(cluster_entry)

        # found a new cluster
        if not matching_clusters:
            pangenome.add_founder(gene)

        # add to existing cluster
        elif len(matching_clusters) is 1:

            clust, entry = matching_clusters[0]
            closest = pangenome.find_closest(clust, gene, entry)
            pangenome.add_to_cluster(clust, gene, closest)

        # merge two or more clusters
        else:

            c1, *matching_clusters = matching_clusters

            recursive_cluster_join(c1, matching_clusters)

    else:
        pangenome.add_founder(gene)

def main():

    args = arguments()

    pangenome = GC.Pangenome()
    counter = 1
    for i in args.input:
        for a in annotations(i):

            add_to_graph(a, pangenome, 0.1)

            G = nx.Graph()
            for c in pangenome.clusters.values():

                if len(c) is 1:
                    G.add_nodes_from(c.nodes())
                else:

                    G.add_edges_from(c.edges())

            nx.draw(G, nx.spring_layout(G), with_labels=True)
            plt.savefig('{}.png'.format(counter))
            plt.close()

            counter += 1
if __name__ == '__main__':
    main()
