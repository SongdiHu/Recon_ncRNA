import os
import sys
from tempfile import NamedTemporaryFile as NTF
from subprocess import call
from multiprocessing import Pool
import itertools

import incubator

NB_PROCS = 24
PATH_SS = 'ss.txt'
PATH_OUT = 'frnakenstein_seeds.txt'

REPETITION = 100

def slave(to_do):
    return (to_do[0], to_do[1], incubator.best_obj_frnakenstein(to_do))

def main():
    l_ss = [x.strip() for x in open(PATH_SS) if x]
    pool = Pool(processes=NB_PROCS)
    to_do = ((s1, s2) for s1, s2 in itertools.combinations(l_ss, 2) for _ in range(REPETITION))
    data = {}

    for s1, s2, seq in pool.imap_unordered(slave, to_do):
        print s1, s2, seq
        data.setdefault((s1, s2), []).append(seq)

    with open(PATH_OUT, 'w') as f:
        first = True
        for s1,s2 in data:
            if not first:
                f.write('\n')
            f.write('%s\n%s\n' % (s1, s2))
            f.write('\n'.join(data[s1,s2]))
            first = False


   
if __name__ == '__main__':
    main()
