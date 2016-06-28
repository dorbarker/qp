#!/usr/bin/env python3

from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from multiprocessing import cpu_count
import GraphController as GC
import argparse
import json
import matplotlib.pyplot as plt
import networkx as nx
import os
import utilities

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', nargs='+', dest='infile',
                        help='Paths to one or more FASTA files')

    parser.add_argument('-t', '--threshold', type=float, default=0.15,
                        help='Maximum distance for cluster inclusion [0.15]')

    parser.add_argument('-p', '--cpus', type=int, default=cpu_count(),
                        help='Number of processor threads to use [all]')

    parser.add_argument('--aln-program', default='mafft',
                        choices=('mafft', 'clustalo'),
                        help='Alignment program [mafft]')

    return parser.parse_args()

def annotations(ffn):

    def genome_basename(n):
        return os.path.splitext(os.path.basename(n))[0]

    with open(ffn) as f:
        for rec in SeqIO.parse(f, 'fasta'):

            yield GC.GeneNode(genome_basename(ffn),
                              rec.id, rec.seq)

def add_to_graph(gene, pangenome, threshold, prog, cpus):

    def recursive_cluster_join(c1, clusters):

        if clusters:

            c2, *clusters  = clusters

            pangenome.join_clusters(c1, c2, gene, prog, cpus)

            recursive_cluster_join(c1, clusters)

    if pangenome.clusters:

        with ThreadPoolExecutor(cpus) as executor:

            f = partial(utilities.cluster_match,
                        gene=gene, threshold=threshold)

            matches = executor.map(f, pangenome.clusters.items())

        matching_clusters = list(filter(None, matches))
        # found a new cluster
        if not matching_clusters:
            pangenome.add_founder(gene)

        # add to existing cluster
        elif len(matching_clusters) is 1:

            clust, entry = matching_clusters[0]
            closest = pangenome.find_closest(clust, gene, entry, prog, cpus)
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

    for i in args.infile:
        for a in annotations(i):

            add_to_graph(a, pangenome, args.threshold,
                         args.aln_program, args.cpus)

    G = nx.Graph()
    for c in pangenome.clusters.values():

        if len(c) is 1:
            G.add_nodes_from(c.nodes())
        else:

            G.add_edges_from(c.edges())

    nx.draw(G, nx.spring_layout(G), with_labels=True)
    plt.savefig('test_data.png')
    plt.close()

if __name__ == '__main__':
    main()
