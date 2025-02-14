import os
import sys
from tempfile import NamedTemporaryFile as NTF
from subprocess import call
from multiprocessing import Pool
from random import shuffle

PATH_DATA = 'frnakenstein_seeds.txt'
PATH_SS = 'ss.txt'
NB_PROCS = 8
#TO_DO_MUT = ['0.01', '0.05',  '0.1', '0.15', '0.2']
TO_DO_MUT = ['0.01', '0.05',  '0.1']

DEPTH = '6'

def slave(args):
    for MUT in TO_DO_MUT:
        i,j,k,seq,sss = args
        name_output = os.path.join('..', 'Data_sank', 
                        'incubated_%s_%s_%s_%s.newick' % (i, j, k, MUT))
        if os.path.isfile(name_output):
            print "PASSED"
            continue
        print "Structs: %s %s at Repetition %s and mut rate %s" % (i,j,k, MUT)
        tmp_file = NTF(delete=False, dir='.', prefix='.tmp', mode='w')
        tmp_file.file.write(sss[i] + '\n' + sss[j] + '\n' + seq)
        tmp_file.file.close()
        print tmp_file.name

        try:
            print name_output
            comd = ['python', 'incubator.py', '-s', tmp_file.name, 
                  '-d', DEPTH, '-o', name_output, '-m', MUT]
            call(comd) 
        except OSError, e:
            print e
            print i, j, k
            sys.exit(1)
        finally:
            os.remove(tmp_file.name)

def main():
    data = {}
    sss = tuple(x.strip() for x in open(PATH_SS))
    with open(PATH_DATA) as f:
        line = f.readline().strip()
        while line:
            if all(x in '(.)' for x in line):
                ss = (sss.index(line), sss.index(f.readline().strip()))
            elif all(x in 'ACGU' for x in line):
                data.setdefault(ss, []).append(line)
            line = f.readline().strip()

    pool = Pool(processes=NB_PROCS)
    to_do = [(i,j,k,seq,sss) for i,j in data for k, seq in enumerate(data[i,j])]
    shuffle(to_do)
    list(pool.imap_unordered(slave, to_do))


if __name__ == '__main__':
    main()
