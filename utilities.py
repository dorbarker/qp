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

class ClusterMatch(object):

    def __init__(self, clusters, genome):

        self.clusters = clusters
        self.genome = genome  # list of GeneNode objects 
        self.fasta = self.assemble_multifasta()

    def assemble_multifasta(self):
        
        def fasta(gene):
            return '>{}\n{}'.format(gene.gene, gene.sequence)

        incoming = '\n'.join(fasta(g) for g in self.genome)
        return '\n'.join((incoming, self.get_representatives()))

    def cluster(self, threshold, cpus):

        with tempfile.TemporaryDirectory() as d:

            fasta_path = os.path.join(d.name, 'clusters.fasta')
            out_path = os.path.join(d.name, 'clusters.out')

            cdhit = ('cd-hit', '-i', fasta_path, '-o', out_path,
                     '-c', str(1.0 - threshold), '-T', str(cpus),
                     '-M', '0', '-d', '0')

            with open(fasta_path, 'w') as f:
                f.write(self.fasta)

            subprocess.check_call(cdhit)

            return self.parse_clusters(out_path + '.clstr')

    def get_representatives(self):

       return ((k, nx.center(v)[0].sequence) for k, v in self.clusters.items())

    def parse_clusters(self, path):

        def get_clust_entry(gene):
            Cluster_Entry(gene, nx.center(self.clusters[gene])[0])

        pt = re.compile('>.*\.\.\.')

        matching_clusters = collections.defaultdict(list)

        with open(path) as f:
            data = f.read().split('>Cluster')

        # get membership for each cluster
        clusts = ([x.strip('>.') for x in re.findall(pt, y)] for y in data if y)

        for c in clusts:

            # genes from the incoming genome
            incoming = [gene for gene in c if gene not in self.clusters]

            # representative genes from the pangenome clusters
            reps = [get_clust_entry(g) for gene in c if g in self.clusters]

            # pair incoming genes with matching pangenome clusters
            for g in incoming:
                matching_clusters[g] = already_clustered

        return matching_clusters
