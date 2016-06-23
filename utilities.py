import subprocess
from multiprocessing import cpu_count
import collections

def calculate_distance(gene1, gene2):

    def format_input(seq1, seq2):
        return '>seq1\n{}\n>seq2\n{}\n'.format(seq1, seq2)

    def hamming(a, b):
        return sum(int(i != j) for i, j in zip(a,b)) / len(a)

    clustalo = ('clustalo',
                '-i', '-',
                '--wrap', str(int(1e9)),
                '--threads', str(cpu_count()))

    aln = subprocess.check_output(clustalo,
                                 input=format_input(gene1.sequence,
                                                    gene2.sequence),
                                 universal_newlines=True)

    return hamming(*filter(lambda x: '>' not in x, aln.split()))

Cluster_Entry = collections.namedtuple('Cluster_Entry', ('cluster', 'entry'))
