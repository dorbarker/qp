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

            yield GC.GeneNode(genome_basename(ffn),
                              rec.id, rec.seq)

def add_to_graph(gene, pangenome, threshold):
    print(gene)
    d = partial(utilities.calculate_distance, gene2=gene)
    if pangenome.clusters:
        for cluster in pangenome.clusters.values():
            try:
                closest = min(nx.center(cluster), key=d)
            except:
                print("EXCEPTION ON", cluster)
                quit()
            gene.update_compared(closest)

            if gene.compared[closest] < threshold:
                pangenome.find_closest(cluster, gene, closest)
                break # need to change later to make cluster joining work
        else:
            pangenome.add_founder(gene)

        print(gene, gene.compared)  # diag

    else:
        pangenome.add_founder(gene)

def main():

    args = arguments()

    pangenome = GC.Pangenome()

    for i in args.input:
        for a in annotations(i):


            try:
                add_to_graph(a, pangenome, 0.5)
            except AttributeError:

                print('\n', pangenome.clusters)
                print(pangenome.clusters['NCTC11168_00001'].nodes())

                print(a)
if __name__ == '__main__':
    main()
