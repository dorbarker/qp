#!/usr/bin/env python3

from Bio import SeqIO
from functools import partial
import argparse
import json
import os
import subprocess
import GraphController as GC
import networkx as nx

def arguments():

    parser = argparse.ArgumentParser()


    return parser.parse_args()

def calculate_distance(seq1, seq2):

    pass

def annotations(annot_dir):

    def genome_basename(n):
        return os.path.splitext(os.path.basename(n))[0]

    ffn_path = partial(os.path.join, annot_dir)
    ffns = (ffn_path(x) for x in os.listdir(annot_dir) if '.ffn' in x)

    for annotation in ffns:
        with open(annotation) as f:
            for rec in SeqIO.parse(f, 'fasta'):
                yield GC.GeneNode(genome_basename(annotation),
                                  rec.id, rec.seq)

def add_to_graph(gene, pangenome, threshold):

    d = partial(calculate_distance, seq1=gene)
    for cluster in pangenome.clusters:
        closest = min(nx.center(cluster), key=d)

def main():

    args = arguments()

    pangenome = GC.Pangenome()

if __name__ == '__main__':
    main()
