#!/usr/bin/env python3

from Bio import SeqIO
from functools import partial
import argparse
import json
import os
import GraphController as GC
import networkx as nx
import utilities

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', nargs='+')

    return parser.parse_args()


def annotations(ffn):

    def genome_basename(n):
        return os.path.splitext(os.path.basename(n))[0]

    with open(ffn) as f:
        for rec in SeqIO.parse(f, 'fasta'):
            print('here')
            yield GC.GeneNode(genome_basename(ffn),
                              rec.id, rec.seq)

def add_to_graph(gene, pangenome, threshold):

    d = partial(utilities.calculate_distance, gene2=gene)
    if pangenome.clusters:
        for cluster in pangenome.clusters.values():
            closest = min(nx.center(cluster), key=d)

            gene.update_compared(closest)

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

    for i in args.input:
        for a in annotations(i):
            print(a.gene)
            add_to_graph(a, pangenome, 0.5)
            print(len(pangenome.clusters))
            print(pangenome.clusters.keys())
if __name__ == '__main__':
    main()
