import subprocess
import collections

Cluster_Entry = collections.namedtuple('Cluster_Entry', ('cluster', 'entry'))

def calculate_distance(gene1, gene2, cpus):

    def format_input(seq1, seq2):
        return '>seq1\n{}\n>seq2\n{}\n'.format(seq1, seq2)

    def hamming(a, b):
        return sum(int(i != j) for i, j in zip(a,b)) / len(a)

    clustalo = ('clustalo',
                '-i', '-',
                '--wrap', str(int(1e9)),
                '--threads', str(cpus))

    aln = subprocess.check_output(clustalo,
                                 input=format_input(gene1.sequence,
                                                    gene2.sequence),
                                 universal_newlines=True)

    return hamming(*filter(lambda x: '>' not in x, aln.split()))

def calculate_search_threads(n_neighbours, cpus):

    if n_neighbours:

        simul_searches = n_neighbours if n_neighbours < cpus else cpus

        per_search = cpus // simul_searches

        out = per_search if per_search > 1 else 1

    else:
        out = cpus

    return out
