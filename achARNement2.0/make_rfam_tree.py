import os
import sys
import numpy as np
from pprint import pprint

PATH_DATA = os.path.join('..')
PATH_FAM_1 = os.path.join(PATH_DATA, 'RF00083_full.txt')
PATH_FAM_2 = os.path.join(PATH_DATA, 'RF00128_full.txt')
PATH_FAM = os.path.join(PATH_DATA, 'RF00083_RF00128_merged.fasta')
PATH_TREE = os.path.join(PATH_DATA, 'theTree.out')

PATH_OUT = os.path.join(PATH_DATA, 'theTree.newick')

FAM_1_SS = '<<...<<<<<<<<<<.......<<<<<...<<.....................................................>>...>>>>>......>>>>>>>.>>.>...>>...........<<..<<<<<.....<<<<..........<<<<...............................................................................>>>>..........>>>>...<<..........................>>.>>>>>.>>............<<<.......................................>>>........................<<<<<<<<<<<.................................................>>>>>>>>>>>.....'

FAM_2_SS = '....<<<<<<<<.........<<<<<<..<<<.................................................>>>>>>>>>.......>>.>>>>.>>..................<<<<<<<..<<<<<........<<<<<<.......................................................>..>>>>>............>>>>>.<<<..........>>>>>>>>>>.....'


def dict_rf(lines):
    d = {}
    for x in lines:
        id, seq = x.split()
        d[id.split("/")[0]] = seq
    return d

def dict_fam(lines):
    d = {}
    name = ''
    seq = ''
    for x in lines:
        x = x.strip()
        if not x:
            continue
        if x.startswith('>'):
            new_name = x[1:]
            #if name == 'AFVR01000019.1':
            if seq:
                d.setdefault(name, []).append(seq.replace('-', '.'))
            name = new_name
            seq = ''
        else:
            seq += x

    d.setdefault(name, []).append(seq.replace('-', '.'))
    return d

def col(l_seqs):
    to_remove = set()
    for seq in l_seqs:
        for i, x in enumerate(seq):
            if x not in 'ACGU':
                to_remove.add(i)
    return to_remove


def rm_pos(seq, forb):
    return ''.join(x for i,x in enumerate(seq) if i not in forb)

def ss2bps(ss):
    l = []
    bps = []
    for i,x in enumerate(ss):
        if x in '<(':
            l.append(i)
        elif x in ')>':
            bps.append((l.pop(), i))
    return bps

def align_map(small, large):
    if ''.join(x for x in small if x != '.') != ''.join(x for x in large if x != '.'):
        print "aligning diff seqs"
        sys.exit(1)
    z = 0
    map = {}
    for i, x in enumerate(small):
        if x == '.':
            continue
        while large[z] == '.':
            z += 1
        map[i] = z
        z += 1
    return map


def fix_ss(fam, fam_ali):
    global FAM_1_SS
    global FAM_2_SS

    to_remove_1 = col([z[0] for z in fam.values()])
    to_remove_2 = col([z[1] for z in fam.values()])
    for l,r in ss2bps(FAM_1_SS):
        if l in to_remove_1 or r in to_remove_1:
            FAM_1_SS = FAM_1_SS[:l] + '.'+ FAM_1_SS[l+1:]
            FAM_1_SS = FAM_1_SS[:r] + '.'+ FAM_1_SS[r+1:]
    for l,r in ss2bps(FAM_2_SS):
        if l in to_remove_2 or r in to_remove_2:
            FAM_2_SS = FAM_2_SS[:l] + '.'+ FAM_2_SS[l+1:]
            FAM_2_SS = FAM_2_SS[:r] + '.'+ FAM_2_SS[r+1:]

    fam1 = ((k, s[0]) for k, s in fam.iteritems())
    fam2 = ((k, s[1]) for k, s in fam.iteritems())
    while 1:
        s1_pos = []
        k1, s = fam1.next()
        z = -1
        for i, x in enumerate(s):
            if x != '.':
                z += 1
            if FAM_1_SS[i] != '.':
                if x== '.':
                    break
                s1_pos.append(z)
        if len(s1_pos) == sum(1 for x in FAM_1_SS if x != '.'):
            break

    while 1:
        s2_pos = []
        k2, s = next(fam2)
        z = 0
        for i, x in enumerate(s):
            if FAM_2_SS[i] != '.':
                if x== '.':
                    break
                s2_pos.append(z)
            if x != '.':
                z ++ 1
        if len(s2_pos) == sum(1 for x in FAM_2_SS if x != '.'):
            break

    m1 = align_map(fam[k1][0], fam_ali[k1][0])
    m2 = align_map(fam[k2][1], fam_ali[k2][1])
    new_bps_1 = [(m1[bp[0]], m1[bp[1]]) for bp in ss2bps(FAM_1_SS)]
    new_bps_2 = [(m2[bp[0]], m2[bp[1]]) for bp in ss2bps(FAM_2_SS)]
    newss1 = ['.']*len(fam_ali.values()[0][0])
    for l,r in new_bps_1:
        newss1[l] = '('
        newss1[r] = ')'
    FAM_1_SS = ''.join(newss1)
    newss2 = ['.']*len(fam_ali.values()[0][0])
    for l,r in new_bps_2:
        newss2[l] = '('
        newss2[r] = ')'
    FAM_2_SS = ''.join(newss2)

def pprint_table(maps):
    table = []
    for x in maps:
        table.append(r'%s & %s\\' % (x, maps[x]))
        table[-1] = table[-1].replace('_', r'\_')
    print '\n'.join(table)

def main():
    tree = [x.strip() for x in open(PATH_TREE) if x.strip()]
    tree, maps = tree[0], dict([x.rsplit(':', 1) for x in tree[2:]])
    maps = {x:maps[x] for x in maps if x in tree }

    pprint_table(maps)
    sys.exit()

    fam1 = {k:v for k,v in dict_rf(open(PATH_FAM_1)).iteritems() if k in maps.values()}
    fam2 = {k:v for k,v in dict_rf(open(PATH_FAM_2)).iteritems() if k in maps.values()}

    """
    with open('fam1.fasta', 'w') as f:
        f.write('\n'.join(">%s\n%s" % (k,'\n'.join(v[i:i+80].replace('.', '-') for i in range(0,len(v),80))) for k,v in fam1.iteritems()))
    with open('fam2.fasta', 'w') as f:
        f.write('\n'.join(">%s\n%s" % (k,'\n'.join(v[i:i+80].replace('.', '-') for i in range(0,len(v),80))) for k,v in fam2.iteritems()))
    sys.exit(0)
    """

    fam_ali = dict_fam(open(PATH_FAM))
    fam_ali = {x:fam_ali[x] for x in maps.values()}



    fam = {x:(fam1[x], fam2[x]) for x in maps.values()}


    fix_ss(fam, fam_ali)

    """
    with open("outfam.txt", 'w') as f:
        out = '\n'.join("%s\t%s\t%s\n" % (x, fam[x][0], fam[x][1]) for x in fam)
        f.write(out)
    sys.exit(0)
    """


    to_remove = col([x for z in fam_ali.values() for x in z])
    for x in fam_ali:
        data = fam_ali[x]
        fam_ali[x] = rm_pos(data[0], to_remove) + '-' + rm_pos(data[1], to_remove)

    global FAM_1_SS
    global FAM_2_SS


    for l,r in ss2bps(FAM_1_SS):
        if l in to_remove or r in to_remove:
            FAM_1_SS = FAM_1_SS[:l] + '.'+ FAM_1_SS[l+1:]
            FAM_1_SS = FAM_1_SS[:r] + '.'+ FAM_1_SS[r+1:]
    for l,r in ss2bps(FAM_2_SS):
        if l in to_remove or r in to_remove:
            FAM_2_SS = FAM_2_SS[:l] + '.'+ FAM_2_SS[l+1:]
            FAM_2_SS = FAM_2_SS[:r] + '.'+ FAM_2_SS[r+1:]

    new_seq = ''
    for i, x in enumerate(FAM_1_SS):
        if i in to_remove:
            continue
        new_seq += x
    FAM_1_SS = new_seq
    new_seq = ''
    for i, x in enumerate(FAM_2_SS):
        if i in to_remove:
            continue
        new_seq += x
    FAM_2_SS = new_seq
    print FAM_1_SS
    print FAM_2_SS
    for k, v in maps.iteritems():
        tree = tree.replace(k, fam_ali[v])

    with open(PATH_OUT, 'w') as f:
        f.write('%s\n%s\n#\n%s' % (FAM_1_SS, FAM_2_SS, tree))

if __name__ == '__main__':
    main()
