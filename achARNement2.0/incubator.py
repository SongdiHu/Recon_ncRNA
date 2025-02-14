import os
import sys
import argparse
from math import fsum, exp
from random import random, choice
from tempfile import NamedTemporaryFile as NTF
from subprocess import check_output
from multiprocessing import Pool

from Bio import Phylo

NUCLEOTIDES = set('ACGU')
NB_PROCS = 24
BOLTZMANN = 0.0019872041
T = 310.15

def random_weighted_sampling(l_samples):

    tot = fsum(x[1] for x in l_samples)
    scaled_weights = [x[1]/tot for x in l_samples]
    rand_nb = random()
    accumulation = 0
    for i,x in enumerate(scaled_weights):
        accumulation += x
        if accumulation > rand_nb:
            return l_samples[i][0]
    return l_samples[-1][0]


def generate_mutante(sequence, mutation_rate, sankoff=True):
    """Given a sequence oand a mutation rate, will generate one mutant"""
    #ACGU
    #0213 2021 1202 3120
    d_weights = {'A':[('C',1), ('G',2), ('U',1)],                                                                                                                                  
                  'C':[('A',1), ('G',1), ('U',2)],                                                                                                                                 
                  'G':[('A',2), ('C',1), ('G',1)],                                                                                                                                 
                  'U':[('A',1), ('C',2), ('G',1)]}                                                                                                                                 


    new_seq = ['']*len(sequence)
    for i,x in enumerate(sequence):
        if random() < mutation_rate:
            if sankoff:
                new_seq[i] = random_weighted_sampling(d_weights[sequence[i]])
            else:
                new_seq[i] = choice(list(NUCLEOTIDES - set(x)))
        else:
            new_seq[i] = x
    return ''.join(new_seq)


def rnafold_mfe_energy(seq):
    """Input is a sequence
    returns the mfe and its energy"""
    tmp_file = NTF(dir='.', prefix='.tmp_', delete=False, mode='w')
    tmp_file.file.write('>\n' + seq)
    tmp_file.file.close()

    try:
        output = check_output(['RNAfold', '--noPS'], stdin=open(tmp_file.name))
    except OSError, e:
        print e.child_traceback
        print e.strerror
        print "Error with RNAfold"
        return 
    finally:
        os.remove(tmp_file.name)


    data = output.split('\n')[2].strip()
    mfe, data = data.split(None, 1)
    dot_pos = data.find('.')
    start = min(i for i in range(dot_pos) if data[i].isdigit())
    end = max(i for i in range(dot_pos, len(data)) if data[i].isdigit())
    return mfe, float(data[start:end+1])


def rnadistance(ss1, ss2, mode='f'):
    """Input is a two structs and mode for RNAdistance
    returns distance"""
    tmp_file = NTF(dir='.', prefix='.tmp_', delete=False, mode='w')
    tmp_file.file.write(ss1 + '\n' + ss2)
    tmp_file.file.close()

    try:
        output = check_output(['RNAdistance', '-D' , mode], stdin=open(tmp_file.name))
    except OSError, e:
        print e.child_traceback
        print e.strerror
        print e
        print "Error with RNAdist"
        sys.exit(1)
    finally:
        os.remove(tmp_file.name)

    dist = float(output[output.find(':')+1:].strip())
    return dist


def slave_fitness(seq_struct):
    """For Pool"""
    seq, struct = seq_struct
    mfe, energy = rnafold_mfe_energy(seq)
    dist = rnadistance(mfe, struct)
    return dist


def dist_fitness(fit):
    return exp((-fit)/(BOLTZMANN*T))


def fitness(l_samples, struct):
    """Returns a list of fitnesses for the sequences given a structure, 
    depending of the dist_fitness function
    ordered as the input"""
    to_do = ((seq, struct) for seq in l_samples)
    #pool = Pool(processes=NB_PROCS)
    dists = [slave_fitness(x) for x in to_do]
    fit =  [dist_fitness(e) for e in dists]
    return fit


def generate_childs(seed, struct, nb_childs, mutation_rate, nb_samples):
    """Given one seed sequence and one struct
    Generates nb_samples random mutants of seed with the mutation rate mutation_rate
    returns nb_childs sequences samples given the fitness function"""
    l_samples = [generate_mutante(seed, mutation_rate) for _ in range(nb_samples)]
    l_fitness = fitness(l_samples, struct)
    tot_fitness = fsum(l_fitness)
    try:
        l_fitness = [f/tot_fitness for f in l_fitness]
    except ZeroDivisionError:
        return [choice(l_samples) for _ in range(nb_childs)]
    return_childs = []
    for _ in range(nb_childs):
        tot = 0
        r = random()
        for i, f in enumerate(l_fitness):
            tot += f
            if tot > r:
                return_childs.append(l_samples[i])
                break
    return return_childs


def rec_populate_tree(tree, l_seed, l_struct, mutation_rate, nb_samples):
    tree.name = '-'.join(l_seed)
    nb_clades = len(tree.clades)
    if not nb_clades:
        return
    new_seed = zip(*[generate_childs(l_seed[i], l_struct[i], nb_clades, mutation_rate, nb_samples)
                for i in range(len(l_seed))])
    for i, c in enumerate(tree.clades):
        rec_populate_tree(c, new_seed[i], l_struct, mutation_rate, nb_samples)


def populate_tree(newick_path, l_seed, l_struct, mutation_rate=0.01, nb_samples=1000):
    tree = Phylo.read(open(newick_path), "newick", rooted=True)
    root = tree.root
    rec_populate_tree(root, l_seed, l_struct, mutation_rate, nb_samples)
    return tree
   

def newick_empty_str_bin_tree(depth):
    if depth == 0:
        return ''
    return '(' + newick_empty_str_bin_tree(depth-1) + ',' +\
            newick_empty_str_bin_tree(depth-1) + ')'


def call_frnakenstein(l_structs):
    tmp_file = NTF(dir='.', prefix='.tmp_', delete=False, mode='w')
    tmp_file.file.write('\n\n'.join(l_structs))
    tmp_file.file.close()

    path_frnakenstein = os.path.join('frnakenstein', 'frnakenstein.py')

    try:
        out = check_output(['python', path_frnakenstein], 
                              stdin=open(tmp_file.name),
                              stderr=open(os.devnull, 'w'))
    except OSError, e:
        print e.child_traceback
        print e.strerror
        print "Error with Frankenstein"
        return 
    finally:
        os.remove(tmp_file.name)

    out = [x.strip() for x in out.split('\n')]

    l_seqs = []
    nb = 0
    while nb < len(out):
        line = out[nb]
        if line.startswith('Sequence'):
            seq = line.split()[-1]
            obj_fit = out[nb+1].split()
            obj, fit = float(obj_fit[1][:-1]), float(obj_fit[-1])
            l_seqs.append((seq, obj, fit))
            nb += 4
        else:
            nb += 1
    return l_seqs


def best_obj_frnakenstein(l_structs):
    return max(call_frnakenstein(l_structs), key=lambda x:x[1])[0]


def main(ss_path, output_path, depth, mutation_rate):
    l_structs = [x.strip() for x in open(ss_path) if x.strip()]
    if all(x in 'ACGU' for x in l_structs[-1]):
        l_seqs = [l_structs[-1]]*(len(l_structs)-1)
        l_structs = l_structs[:-1]
        print "Seqs retrieved:", l_seqs
    else:
        l_seqs = [best_obj_frnakenstein(l_structs)]*len(l_structs)

    try:
        tmp_tree_path = NTF(dir='.', prefix='.tmptree_', delete=False, mode='w')
        tmp_tree_path.file.write(newick_empty_str_bin_tree(depth) + ';')
        tmp_tree_path.file.close()
        tree = populate_tree(tmp_tree_path.name, l_seqs, l_structs, mutation_rate=mutation_rate)
    except:
        print sys.exc_info()[0]
        print "Error with populate_tree"
        sys.exit(1)
    finally:
        os.remove(tmp_tree_path.name)
    Phylo.write(tree, output_path, 'newick')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='jesuideshonore!') 
    parser.add_argument('-s', '--secondary_structure_path', type=str,
                        help='''path to file containing the secondary 
                        structures, one per line. 
                        If a sequence is given in the third line, 
                        it will be used as seed instead of generating 
                        one with Frnakenstein.''',
                        dest='ss_path', required=True)
    parser.add_argument('-d', '--depth', type=int,
                        help='depth of binary tree',
                        dest='depth', required=True)
    parser.add_argument('-m', '--mutation', type=float,
                        help='mutation rate',
                        dest='mutation_rate', required=True)
    parser.add_argument('-o', '--output', type=str,
                        help='''path to output file''',
                        dest='output_path', required=True)
    args = parser.parse_args()


    main(args.ss_path, args.output_path, args.depth, args.mutation_rate)
