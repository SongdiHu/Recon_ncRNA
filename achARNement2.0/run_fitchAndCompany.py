# Python3.4
import os
import multiprocessing
import time
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool
from tempfile import NamedTemporaryFile as NTF

# PATH_PYTHON3 = os.path.join('/', 'home', 'hus', '.virtualenvs', 'Recon_ncRNA', 'bin', 'python3')
# /home/hus/anaconda3/envs/recon_rna_310/bin
# PATH_PYTHON3 = os.path.join('/', 'mnt', 'd', 'Recon_ncRNA', 'venv', 'bin', 'python')
PATH_PYTHON3 = os.path.join('/', 'mnt', 'd', 'Projects', 'Recon_ncRNA', 'venv', 'bin', 'python')
# PATH_PYTHON3 = os.path.join('/', 'home', 'hus', 'anaconda3', 'envs', 'recon_rna_310', 'bin', 'python3')
NB_PROCS = multiprocessing.cpu_count()

PATH_SS = os.path.join('..', 'Resources', 'ss.txt')
PATH_DATA = os.path.join('..', 'Data_sank')
PATH_RESULTS = os.path.join('..', 'Data_sank_results_bf')
# ALGOS = list(map(str, range(1, 6)))
ALGOS = ["5"]


def worker(args):
    # print(os.environ['PYTHONPATH'])
    # exit()
    t, ss1, ss2, nb = args
    struct_pair = str(t[1][0])+str(t[1][1])
    t = t[0]
    # print("structure pair:", ss1,ss2)
    path_out = os.path.join(PATH_RESULTS, t[:t.rfind('.')] + '_%s_.txt' % nb)
    # print(path_out)

    if os.path.isfile(path_out):
        if os.stat(path_out).st_size > 0:
            return
    else:
        open(path_out, 'w').close()

    temp_file = NTF(prefix='.tmp', dir='.', delete=False, mode='w')
    temp_file.write('\n'.join([ss1, ss2]))
    temp_file.write('\n#\n')
    with open(os.path.join(PATH_DATA, t)) as f:
        temp_file.write(f.read())
    temp_file.close()

    try:
        check_call([PATH_PYTHON3, 'achARNement_V2.py', '-e', nb, '2', temp_file.name, '0.001',  struct_pair, '0'],
                   stdout=open(path_out, 'w'), timeout=172800)
        # check_call([PATH_PYTHON3, 'achARNement_V2.py', nb, '2', temp_file.name, '0.003', struct_pair],
        #            stdout=open(path_out, 'w+'),
        #            env={'PYTHONPATH': '/home/hus/.virtualenvs/Recon_ncRNA/lib/python3.10/site-packages:'
        #                               '/mnt/d/Projects/Recon_ncRNA/achARNement2.0/tree-diet/python:'
        #                               '/mnt/d/Projects/Recon_ncRNA/achARNement2.0/tree-diet/lib'})
    except CalledProcessError as e:
        print(e)
        print("Error with fitchAndCompany", t)
    finally:
        print(t, nb, "done.")
        os.remove(temp_file.name)
        pass


# 246,204,222
if __name__ == '__main__':
    list_ss = [x.strip() for x in open(PATH_SS, 'r')]
    trees = [(x.strip(), list(map(int, x[:x.find('.')].split('_')[1:3]))) for x in os.listdir(PATH_DATA)
             if x.endswith('newick') and os.stat(os.path.join(PATH_DATA, x)).st_size > 0 and
             x[x.rfind('_') + 1:x.rfind('.')] in ('0.01', '0.05', '0.1', '0.15', '0.2')]

    # print(trees)
    # exit()
    # # test
    # print("PATH_DATA:", os.path.exists(PATH_DATA))
    # print("PATH_RESULTS:", os.path.exists(PATH_RESULTS))
    # print("list_ss:", list_ss)
    # print("ALGOS:", ALGOS)
    # print("ALGO #:", ALGOS)
    to_do = []
    for nb in ALGOS:
        for t in trees:
            mut = t[0][t[0].rfind('_') + 1:t[0].rfind('.')]
            if (mut not in ('0.15', '0.2')) or (mut == '0.15' and nb in ('4', '5')) or (mut == '0.2' and nb == '5'):
                to_do.append((t, list_ss[t[1][0]], list_ss[t[1][1]], nb))

    # os.environ['PYTHONPATH'] += ":/home/hus/anaconda3/envs/recon_rna_310/lib/python3.10/site-packages:/home/hus/Documents/Projects/Recon_ncRNA/tree-diet/python:/home/hus/Documents/Projects/Recon_ncRNA/tree-diet/lib"
    # print("TODO list:\n")
    # for params in to_do:
    #     print(params)

    print("# of cores to use:", NB_PROCS)
    start_time = time.time()
    # pool = Pool(processes=NB_PROCS)
    # list(pool.imap_unordered(worker, to_do))
    with Pool(NB_PROCS) as p:
        p.map(worker, to_do)

    print("time elapsed:", time.time() - start_time)
