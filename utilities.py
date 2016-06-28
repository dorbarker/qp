import subprocess
import collections
from Bio import SeqIO
import os
import networkx as nx
import tempfile

Cluster_Entry = collections.namedtuple('Cluster_Entry', ('cluster', 'entry'))

def format_input(seq1, seq2):
    return '>seq1\n{}\n>seq2\n{}\n'.format(seq1, seq2)

def calculate_distance(gene1, gene2,prog, cpus):

    def hamming(a, b):
        return sum(int(i != j) for i, j in zip(a,b)) / len(a)

    def aln_parse(aln):

        lines = aln.strip().split()

        sec_head = lines.index('>seq2')

        return ''.join(lines[1:sec_head]), ''.join(lines[sec_head + 1:])

    if prog == 'mafft':

        run = ('mafft', '--quiet', '--retree', '1', '--thread', str(cpus), '-')
    else:

        run = ('clustalo', '--threads', str(cpus), '-i', '-')

    fasta = format_input(gene1.sequence, gene2.sequence)

    aln = subprocess.check_output(run, input=fasta, universal_newlines=True)

    return hamming(*aln_parse(aln))

def calculate_search_threads(n_neighbours, cpus):

    if n_neighbours:

        simul_searches = n_neighbours if n_neighbours < cpus else cpus

        per_search = cpus // simul_searches

        out = per_search if per_search > 1 else 1

    else:
        out = cpus

    return out

def cluster_match(clustername_cluster, gene, threshold):

    clustername, cluster = clustername_cluster
    centre = nx.center(cluster)[0]  # centre node of the cluster

    fasta = format_input(gene.sequence, centre.sequence)

    with tempfile.NamedTemporaryFile(mode='r+') as i:
        
        i.write(fasta)

        with tempfile.NamedTemporaryFile(mode='r+') as o:
            i.seek(0)
            cdhit = ('cdhit',
                     '-i', i.name,
                     '-o', o.name,
                     '-c', str(1.0 - threshold))

            subprocess.call(cdhit, stdout=subprocess.DEVNULL)

            data = str(o.read()).split()

    if len([x for x in data if '>' in x]) is 1:

        return Cluster_Entry(clustername, centre)
