#!/usr/bin/env python3

from Bio import SeqIO
from functools import partial
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

    if pangenome.clusters:
        for cluster in pangenome.clusters.values():

            closest = min(nx.center(cluster), key=gene.update_compared)

            if gene.compared[closest] < threshold:
                pangenome.find_closest(cluster, gene, closest)
                break # need to change later to make cluster joining work
        else:
            pangenome.add_founder(gene)

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

            nx.draw(G, nx.random_layout(G), with_labels=True)
            plt.savefig('{}.png'.format(counter))
            plt.close()

            counter += 1
if __name__ == '__main__':
    main()
