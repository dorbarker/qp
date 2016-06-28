import subprocess
import collections
from Bio import SeqIO
import os
import networkx as nx

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

def cluster_match(pangenome, gene, threshold, cpus):

    def determine_cluster_entry(cluster_keys, genes, match):

        cluster = cluster_keys[genes.index(match)]
        return Cluster_Entry(cluster, pangenome.gene_lookup[match])

    if not os.access('/tmp/qptemp/', os.F_OK):
        os.mkdir('/tmp/qptemp/')

    temp_in = '/tmp/qptemp/{}.fasta'.format(gene.gene)
    temp_out = '/tmp/qptemp/{}_out.fasta'.format(gene.gene)

    cdhit = ('cdhit',
             '-i', temp_in,
             '-o', temp_out,
             '-c', str(1.0 - threshold),
             '-T', str(cpus))

    clusts = list(pangenome.clusters.keys())
    centres = [nx.center(c)[0] for c in pangenome.clusters.values()]
    genes = [g.gene for g in centres]
    fastas = ('>{}\n{}'.format(x.gene, x.sequence) for x in centres)

    fasta = '>{}\n{}\n'.format(gene.gene, gene.sequence) + '\n'.join(fastas)

    with open(temp_in, 'w') as o:
        o.write(fasta)

    subprocess.check_call(cdhit, stdout=subprocess.DEVNULL)

    with open(temp_out, 'r') as f:
        d = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

    similar = set(genes) - set(d.keys())

    matches = [determine_cluster_entry(clusts, genes, x) for x in similar]

    # tidy up temp directory
    for i in (temp_in, temp_out, temp_out + '.clstr'):
        os.remove(i)

    return matches
