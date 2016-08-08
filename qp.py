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

def annotations(ffn, pangenome):

    def genome_basename(n):
        return os.path.splitext(os.path.basename(n))[0]

    with open(ffn) as f:
        for rec in SeqIO.parse(f, 'fasta'):
            if rec.id in pangenome.gene_lookup:
                continue
            yield rec.id, GC.GeneNode(genome_basename(ffn),
                                      rec.id, rec.seq)

def cluster_incoming_genome(pangenome, annotations, threshold, cpus):

    incoming_genome =  utilities.ClusterMatch(pangenome, annotations)

    return incoming_genome.cluster(threshold, cpus)

def add_to_graph(gene, pangenome, matches, threshold, prog, cpus):

    def recursive_cluster_join(c1, clusters):

        if clusters:

            c2, *clusters  = clusters

            if c1 != c2:
                pangenome.join_clusters(c1, c2, gene, prog, cpus)

            recursive_cluster_join(c1, clusters)

    if pangenome.clusters:

        # found a new cluster
        if not matches:
            pangenome.add_founder(gene)

        # add to existing cluster
        elif len(matches) is 1:

            clust, entry = matches[0]
            closest = pangenome.find_closest(clust, gene, entry, prog, cpus)
            pangenome.add_to_cluster(clust, gene, closest)

        # merge two or more clusters
        else:

            c1, *matching_clusters = matches

            recursive_cluster_join(c1, matches)

    else:
        pangenome.add_founder(gene)

def main():
    import time
    args = arguments()

    pangenome = GC.Pangenome()

    genomes = 1

    ref = time.time()
    for i in args.infile:
        annots = dict(annotations(i, pangenome))
        matching_clusters = cluster_incoming_genome(pangenome, annots,
                                                    args.threshold, args.cpus)

        for a in annots:

            gene = annots[a]
            add_to_graph(gene, pangenome, matching_clusters[a], args.threshold,
                         args.aln_program, args.cpus)

        now = time.time()
        print(len(pangenome.clusters))
        print(now-ref)
        ref = now
        genomes += 1

    G = nx.Graph()
    for c in pangenome.clusters.values():

        if len(c.nodes()) is 1:
            G.add_nodes_from(c.nodes())
        else:

            G.add_edges_from(c.edges())

    print("Edges:", len(G.edges()), "Nodes:", len(G.nodes()), "Clusters:", len(pangenome.clusters))
    #nx.draw(G, nx.spring_layout(G), with_labels=False, edge_color='k')
    #plt.savefig('test_data.png')
    #plt.close()

if __name__ == '__main__':
    main()
