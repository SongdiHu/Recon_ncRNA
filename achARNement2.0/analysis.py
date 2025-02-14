import os
import sys
import re
from math import log, fsum

from numpy import std, average
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

MUT = ['0.01', '0.05', '0.1', '0.15', '0.2']
num_algs = 5

DATA = os.path.join('..', 'Data_sank_results')
re_tot = r'Total: (\d+) errors / (\d+) nucleotides \((\d+\.\d+) % error\)'  # For the best
re_tot_uns = r'Total for unstructured positions: (\d+) errors / (\d+) nucleotides \((\d+.\d+) % error\)'
re_tot_str = r'Total for structured positions: (\d+) errors / (\d+) nucleotides \((\d+.\d+) % error\)'
re_avg = r'Average \(over (\d+) optimal sequences\): (\d+.\d+) aver. errors per optimal sequence of length = (\d+.\d+) nucleotides \((\d+.\d+) % error\)'
re_Uns = r'Average for unstructured positions: (\d+.\d+) aver. errors per optimal sequence, considering (\d+.\d+) nucleotides per seq. \((\d+.\d+) % error\)'
re_Struct = r'Average for structured positions: (\d+.\d+) aver. errors per optimal sequence, considering (\d+.\d+) nucleotides per seq. \((\d+.\d+) % error\)'
re_tot_err = r'Total: (\d+) errors'  # Best root 4 next lines
re_opt = r'Average \(over (\d+) optimal sequences\): (\d+.\d+) aver. errors per optimal sequence of length = 100.0 nucleotides \((\d+.\d+) % error\)'  # All optimales
re_opt_uns = r'Average for unstructured positions: (\d+.\d+) aver. errors per optimal sequence, considering (\d+.\d+) nucleotides per seq. \((\d+.\d+) % error\)'
re_opt_str = r'Average for structured positions: (\d+.\d+) aver. errors per optimal sequence, considering (\d+.\d+) nucleotides per seq. \((\d+.\d+) % error\)'


# Average pair de structures
# Taux d'erreurs Y / mutation X
# % Total, Uns, Struct

def parser(lines):
    data = {x: [] for x in
            ('tot', 'tot_uns', 'tot_str', 'avg', 'Uns', 'Struct', 'tot_err', 'opt', 'opt_uns', 'opt_str')}
    algos = zip(('tot', 'tot_uns', 'tot_str', 'avg', 'Uns', 'Struct', 'tot_err', 'opt', 'opt_uns', 'opt_str'),
                (re_tot, re_tot_uns, re_tot_str, re_avg, re_Uns, re_Struct, re_tot_err, re_opt, re_opt_uns, re_opt_str))

    z = 0
    for x in lines:
        if "-For the root only:" in x:
            z += 1
        if z == 3: break
        for name, alg in algos:
            g = re.match(alg, x)
            if g:
                data[name].append(g.groups())
                continue
    return data


def root_parser(lines):
    found = False
    z = 0
    tot = None
    for x in lines:
        if "-For the root only:" in x:
            z += 1
        if z == 3:
            found = True
        if found:
            g = re.match(re_opt, x)
            if g:
                tot = g.groups()[0]
    return tot


def graph_by_alg(data):
    def plot_honor(data, color):
        tot = []
        tot_err = []
        str = []
        str_err = []
        uns = []
        uns_err = []
        ind_data = data.set_index('mut')

        all_plot = [('total', tot, tot_err),
                    ('structured', str, str_err),
                    ('unstructured', uns, uns_err)]

        for mut in MUT:
            try:
                tot.append(ind_data.loc[mut, 'tot'])
                # print tot[-1]
                tot_err.append(ind_data.loc[mut, 'tot_err'])
                str.append(ind_data.loc[mut, 'str'])
                str_err.append(ind_data.loc[mut, 'str_err'])
                uns.append(ind_data.loc[mut, 'uns'])
                uns_err.append(ind_data.loc[mut, 'uns_err'])
            except KeyError:
                print mut

        cols = {'total': 'k', 'structured': 'red', 'unstructured': 'blue'}
        for n, p, e in all_plot:
            if len(p) != 3:
                continue
            plt.errorbar([0.01, 0.05, 0.1, 0.15, 0.2], p, yerr=e, label=n, color=cols[n])
        plt.xlabel([0.01, 0.05, 0.1, 0.15, 0.2])
        plt.legend()

    g = {}
    for name in data:
        _, s1, s2, rep, mut, alg, _ = name.split('_')
        g.setdefault((s1, s2, alg, mut), {})
        g[s1, s2, alg, mut].setdefault('tot', []).append(
            float(data[name]['tot'][-1][-1]))
        g[s1, s2, alg, mut].setdefault('Uns', []).append(
            float(data[name]['Uns'][-1][-1]))
        g[s1, s2, alg, mut].setdefault('Struct', []).append(
            float(data[name]['Struct'][-1][-1]))

    d = {'ss': [], 'alg': [], 'tot': [],
         'tot_err': [], 'uns': [], 'uns_err': [],
         'str': [], 'str_err': [],
         'mut': []}
    for s1, s2, alg, mut in g:
        d['mut'].append(mut)
        d['ss'].append(s1 + s2)
        d['alg'].append(alg)
        d['tot'].append(average(g[s1, s2, alg, mut]['tot']))
        d['tot_err'].append(std(g[s1, s2, alg, mut]['tot']))
        d['str'].append(average(g[s1, s2, alg, mut]['Struct']))
        d['str_err'].append(std(g[s1, s2, alg, mut]['Struct']))
        d['uns'].append(average(g[s1, s2, alg, mut]['Uns']))
        d['uns_err'].append(std(g[s1, s2, alg, mut]['Uns']))
    data = pd.DataFrame(d)

    sns.set(font_scale=1.2)
    g = sns.FacetGrid(data, col="ss", row='alg', row_order=['1', '2', '3', '4', '5'],
                      col_order=['01', '02', '12'], aspect=1.236)
    g.map_dataframe(plot_honor).set_axis_labels("Mutation rate", "Percentage Error").add_legend().set(
        xticks=[0.01, 0.05, 0.1, 0.15, 0.2])

    plt.tight_layout()
    plt.savefig('test.pdf')
    plt.close()
    # plt.show()


def _graph_by_Region(data):
    def plot_honor(data, color):
        d_algs_val = {str(x): [] for x in range(1, num_algs + 1)}
        d_algs_err = {str(x): [] for x in range(1, num_algs + 1)}
        d_muts = {str(x): [] for x in range(1, num_algs + 1)}
        ind_data = data.set_index('alg')

        for x in range(1, num_algs + 1):
            x = str(x)
            m = ind_data.loc[x, 'mut']
            m = {z: i for i, z in enumerate(m)}
            vals = ind_data.loc[x, 'vals']
            errs = ind_data.loc[x, 'errs']
            for mut in MUT:
                if mut in m:
                    d_algs_val[x].append(vals[m[mut]])
                    d_algs_err[x].append(errs[m[mut]])
                    d_muts[x].append(mut)

        cols = dict(zip('12345', 'bgkry'))
        for x in range(1, num_algs + 1):
            x = str(x)
            print d_algs_val[x]
            eb = plt.errorbar(d_muts[x], d_algs_val[x], yerr=d_algs_err[x], label=x, color=cols[x],
                              alpha=0.5 if x != '1' else 1, linewidth=3.0)
            eb[-1][0].set_linestyle('--')
            eb[-1][0].set_linewidth(1)
        # plt.xlabel([0.01, 0.05, 0.1])
        plt.legend()

    g = {}
    for name in data:
        _, s1, s2, rep, mut, alg, _ = name.split('_')
        g.setdefault((s1, s2, alg, mut), {})
        g[s1, s2, alg, mut].setdefault('tot', []).append(
            float(data[name]['tot'][-1][-1]))
        g[s1, s2, alg, mut].setdefault('Uns', []).append(
            float(data[name]['Uns'][-1][-1]))
        g[s1, s2, alg, mut].setdefault('Struct', []).append(
            float(data[name]['Struct'][-1][-1]))

    d = {'ss': [], 'alg': [], 'errs': [],
         'mut': [], 'Region': [], 'vals': []}
    for s1, s2, alg, mut in g:
        for x in ('tot', 'Struct', 'Uns'):
            d['mut'].append(mut)
            d['ss'].append(s1 + s2)
            d['alg'].append(alg)
            d['Region'].append(x)
            d['vals'].append(average(g[s1, s2, alg, mut][x]))
            d['errs'].append(std(g[s1, s2, alg, mut][x]))
    data = pd.DataFrame(d)

    sns.set(font_scale=1.2)
    g = sns.FacetGrid(data, col="ss", row='Region', row_order=['tot', 'Struct', 'Uns'],
                      col_order=['01', '02', '12'], aspect=1.236)
    g.map_dataframe(plot_honor).set_axis_labels("Mutation rate", "Percentage Error").add_legend().set(
        xticks=[0.01, 0.05, 0.1, 0.15, 0.2], yticks=[0, 2, 4, 6, 8])
    plt.tight_layout()
    plt.savefig('error.png')
    plt.close()

    # plt.show()


def graph_by_Region(data):
    def plot_honor(data, color):
        d_algs_val = {str(x): [] for x in range(1, num_algs + 1)}
        d_algs_err = {str(x): [] for x in range(1, num_algs + 1)}
        d_muts = {str(x): [] for x in range(1, num_algs + 1)}
        ind_data = data.set_index('alg')

        for x in range(1, num_algs + 1):
            x = str(x)
            m = ind_data.loc[x, 'mut']
            m = {z: i for i, z in enumerate(m)}
            vals = ind_data.loc[x, 'vals']
            errs = ind_data.loc[x, 'errs']
            for mut in MUT:
                if mut in m:
                    d_algs_val[x].append(vals[m[mut]])
                    d_algs_err[x].append(errs[m[mut]])
                    d_muts[x].append(mut)

        cols = dict(zip('12345', 'bgkry'))
        for x in range(1, num_algs + 1):
            if x == 3:
                continue
            x = str(x)
            # print d_algs_val[x]
            eb = plt.errorbar(d_muts[x], d_algs_val[x], yerr=d_algs_err[x], label=x, color=cols[x],
                              alpha=0.5 if x != '1' else 1, linewidth=3.0)
            eb[-1][0].set_linestyle('--')
            eb[-1][0].set_linewidth(1)
        # plt.xlabel([0.01, 0.05, 0.1])
        plt.legend()

    g = {}
    for name in data:
        # print name, data[name]
        _, s1, s2, rep, mut, alg, _ = name.split('_')
        g.setdefault((s1, s2, alg, mut), {})
        g[s1, s2, alg, mut].setdefault('tot', []).append(
            float(data[name]['tot'][-1][-1]))
        g[s1, s2, alg, mut].setdefault('Uns', []).append(
            float(data[name]['Uns'][-1][-1]))
        g[s1, s2, alg, mut].setdefault('Struct', []).append(
            float(data[name]['Struct'][-1][-1]))

    d = {'ss': [], 'alg': [], 'errs': [],
         'mut': [], 'Region': [], 'vals': []}
    for s1, s2, alg, mut in g:
        for x in ('tot', 'Struct', 'Uns'):
            if x == 'tot':
                continue
            d['mut'].append(mut)
            d['ss'].append(s1 + s2)
            d['alg'].append(alg)
            d['Region'].append(x)
            d['vals'].append(average(g[s1, s2, alg, mut][x]))
            d['errs'].append(std(g[s1, s2, alg, mut][x]))
    data = pd.DataFrame(d)
    # print(type(data.loc[0].at['mut']))  # test
    # print(data)

    sns.set(font_scale=1.2)
    g = sns.FacetGrid(data, col="ss", row='Region', row_order=['Struct', 'Uns'], col_order=['01', '02', '12'],
                      aspect=1.236)
    g.map_dataframe(plot_honor)
    g.set_axis_labels("Mutation rate", "Percentage Error")
    g.set(xticks=["0.01", "0.05", "0.1", "0.15", "0.2"], yticks=[0, 4, 8, 12, 16, 20])
    g.fig.subplots_adjust(left=0.05, top=0.95)
    hand = [mpatches.Patch(color=c, label=alg)
            for c, alg in zip('bgkry', ['Fitch', 'Sankoff', 'achARNement 1.0', 'achARNement 2.0'])]
    plt.legend(handles=hand, loc="upper left", bbox_to_anchor=(0, 2.25), framealpha=0.0)
    plt.savefig('error.png')
    plt.close()

    # plt.show()


def graph_tot_sol(data):
    def plot_honor(data, color):
        d_algs_val = {str(x): [] for x in range(1, num_algs + 1)}
        d_algs_err = {str(x): [] for x in range(1, num_algs + 1)}
        d_muts = {str(x): [] for x in range(1, num_algs + 1)}
        ind_data = data.set_index('alg')

        for x in range(1, num_algs + 1):
            x = str(x)
            m = ind_data.loc[x, 'mut']
            m = {z: i for i, z in enumerate(m)}
            vals = ind_data.loc[x, 'vals']
            errs = ind_data.loc[x, 'errs']
            for mut in MUT:
                if mut in m:
                    d_algs_val[x].append(vals[m[mut]])
                    d_algs_err[x].append((errs[m[mut]]))
                    d_muts[x].append(mut)

        cols = dict(zip('12345', 'bgkry'))
        print "TOTAL:", d_algs_val
        for x in range(1, num_algs + 1):
            if x == 3:
                continue
            x = str(x)
            plt.yscale('log', nonposy='clip')
            eb = plt.errorbar(d_muts[x], d_algs_val[x], yerr=d_algs_err[x], label=x, color=cols[x],
                              alpha=0.5 if x != '1' else 1, linewidth=3.0)
            eb[-1][0].set_linestyle('--')
            eb[-1][0].set_linewidth(1)
        # plt.xlabel([0.01, 0.05, 0.1])
        plt.legend()

    g = {}
    for name in data:
        _, s1, s2, rep, mut, alg, _ = name.split('_')
        g.setdefault((s1, s2, alg, mut), []).append(
            float(data[name]['avg'][-1][0]))

    d = {'ss': [], 'alg': [], 'errs': [],
         'mut': [], 'vals': []}
    for s1, s2, alg, mut in g:
        d['mut'].append(mut)
        d['ss'].append(s1 + s2)
        d['alg'].append(alg)
        d['vals'].append(average(g[s1, s2, alg, mut]))
        d['errs'].append(std(g[s1, s2, alg, mut]))
    data = pd.DataFrame(d)

    sns.set(font_scale=1.2)
    g = sns.FacetGrid(data, col="ss", row_order=['tot', 'Struct', 'Uns'],
                      col_order=['01', '02', '12'], aspect=1.236)
    g.map_dataframe(plot_honor)
    g.set_axis_labels("Mutation rate", "Number Sequences")
    g.set(xticks=["0.01", "0.05", "0.1", "0.15", "0.2"])

    hand = [mpatches.Patch(color=c, label=alg)
            for c, alg in zip('bgkry', ['Fitch', 'Sankoff', 'achARNement 1.0', 'achARNement 2.0'])]
    plt.legend(handles=hand, loc="upper left", bbox_to_anchor=(0, 1.08), framealpha=0.0)

    plt.tight_layout()
    plt.savefig('nb_all.png')
    plt.close()

    # plt.show()


def graph_nb_sol(data):
    def plot_honor(data, color):
        d_algs_val = {str(x): [] for x in range(1, num_algs + 1)}
        d_algs_err = {str(x): [] for x in range(1, num_algs + 1)}
        d_muts = {str(x): [] for x in range(1, num_algs + 1)}
        ind_data = data.set_index('alg')

        for x in range(1, num_algs + 1):
            x = str(x)
            m = ind_data.loc[x, 'mut']
            m = {z: i for i, z in enumerate(m)}
            vals = ind_data.loc[x, 'vals']
            errs = ind_data.loc[x, 'errs']
            for mut in MUT:
                if mut in m:
                    d_algs_val[x].append(vals[m[mut]])
                    d_algs_err[x].append(errs[m[mut]])
                    d_muts[x].append(mut)
        print "ROOT", d_algs_val
        cols = dict(zip('12345', 'bgkry'))
        for x in range(1, num_algs + 1):
            if x == 3:
                continue
            x = str(x)
            eb = plt.errorbar(d_muts[x], d_algs_val[x], yerr=d_algs_err[x], label=x, color=cols[x],
                              alpha=0.5 if x != '1' else 1, linewidth=3.0)
            eb[-1][0].set_linestyle('--')
            eb[-1][0].set_linewidth(1)
            plt.yscale('log', nonposy='clip')
        # plt.xlabel([0.01, 0.05, 0.1])
        plt.legend()

    g = {}
    for name in data:
        _, s1, s2, rep, mut, alg, _ = name.split('_')
        # print(s1, s2, alg, mut)
        # print(data[name])
        g.setdefault((s1, s2, alg, mut), []).append(float(data[name]))

    d = {'ss': [], 'alg': [], 'errs': [],
         'mut': [], 'vals': []}
    for s1, s2, alg, mut in g:
        d['mut'].append(mut)
        d['ss'].append(s1 + s2)
        d['alg'].append(alg)
        d['vals'].append(average(g[s1, s2, alg, mut]))
        d['errs'].append(std(g[s1, s2, alg, mut]))
    data = pd.DataFrame(d)

    sns.set(font_scale=1.2)
    g = sns.FacetGrid(data, col="ss", col_order=['01', '02', '12'], aspect=1.236)
    g.map_dataframe(plot_honor)
    g.set_axis_labels("Mutation rate", "Number Sequences")
    g.set(xticks=["0.01", "0.05", "0.1", "0.15", "0.2"])
    hand = [mpatches.Patch(color=c, label=alg)
            for c, alg in zip('bgkry', ['Fitch', 'Sankoff', 'achARNement 1.0', 'achARNement 2.0'])]
    plt.legend(handles=hand, loc="upper left", bbox_to_anchor=(0, 1.08), framealpha=0.0)
    plt.tight_layout()
    plt.savefig('root_nb.png')
    plt.close()

    # plt.show()


def graph_root(path_data):
    def retrieve_root(path):
        seqs = set()
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if all(x in 'ACGU' for x in line):
                    seqs.add(line)
                if "Family 1:" in line:
                    print seqs, path, len(seqs)
                    return seqs

    to_do = []
    for x in os.listdir(path_data):
        alg = x[-6]
        if alg != '1':
            continue
        to_do.append(x)
    d = {}
    for x in to_do:
        _, s1, s2, rep, mut, alg, _ = x.split('_')
        sank = os.path.join(path_data, x[:-6] + '2' + x[-5:])
        four = os.path.join(path_data, x[:-6] + '4' + x[-5:])
        if not os.path.isfile(sank) or not os.path.isfile(four):
            continue
        s_sank = retrieve_root(sank)
        s_four = retrieve_root(four)
        d.setdefault((s1, s2), {}).setdefault(mut, []).append(fsum(1 for z in s_four if z in s_sank) / len(s_four))
    print d


def main():
    # graph_root(DATA)
    # sys.exit()

    data = {}
    total = {}
    for file_name in os.listdir(DATA):
        if not file_name.endswith('txt') or os.stat(os.path.join(DATA, file_name)).st_size == 0:
            continue
        data[file_name] = parser(open(os.path.join(DATA, file_name)))
        total[file_name] = root_parser(open(os.path.join(DATA, file_name)))
    # print data
    # exit()
    # with open('data_sank.cPickle', 'w') as f:
    #     data = cPickle.dump((data, total), f, -1)

    # with open('data_sank.cPickle') as f:
    #     data, total = cPickle.load(f)

    graph_by_Region(data)
    
    graph_tot_sol(data)

    graph_nb_sol(total)


if __name__ == '__main__':
    main()
