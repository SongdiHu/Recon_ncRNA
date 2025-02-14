#This script reads 2 RFAM alignments and outputs the list of species that are found in both files
import os
import sys
from Bio import SeqIO
import itertools
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

def name(id):
    if '.' in id:
        id = id[:id.find('.')]
    url = 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % id
    try:
        f = urlopen(url)
    except HTTPError as e:
        return "HTTP Error:" + str(e.code) + "--> failed to find name for id: " + id
        #print("PDB url %s: download failed" % url)
        #sys.exit(1)
    except URLError as e:
        return "URL Error:" + e.reason + "--> failed  to find name for id: " + id
        #print("PDB url %s: download failed" %  url)
        #sys.exit(1)
    for x in f:
        x = x.decode().strip()
        if '<title>' in x:
            return x[x.find('>')+1:-len(r' - Nucleotide - NCBI<\title>')]


#Reads an RFAM alignment and puts it in a dictionary {"name": alignment. Returns the dictionary.}
def readAlignment(filename):

    rfamFile = open(filename, 'r')

    nameToAlignDict = {}
    read_seq = False
    sequence_tmp = ""
    spec_name = ""
    for line in rfamFile:
        if '/' in line:
            if read_seq:
                nameToAlignDict[spec_name] = sequence_tmp
            spec_name = line.split("/")[0][1:]
            read_seq = True
            sequence_tmp = ""
        elif read_seq:
            sequence_tmp += line.strip()

    return nameToAlignDict

#Finds and returns a list of keys that are present in both dicts.
def intersectDictKeys(dict1, dict2):

    listNamesIntersected = []

    for key in dict1:
        if key in dict2:
            listNamesIntersected.append(key)

    return listNamesIntersected

if __name__ == '__main__':
    #print(name('AFGV00000000'))

    if len(sys.argv) < 2:
        print("USAGE: python rfamAlignFile1 rfamALignFile2")
        sys.exit()

    files = [f for f in os.listdir(path=sys.argv[1]) if f.endswith(".fa")]
    # print(files)
    combs = list(itertools.combinations(files, 2))
    intersectionOfNames = []
    family_pair = []
    for comb in combs:
        # print(comb[0], comb[1])
        dictFam1 = readAlignment(sys.argv[1] + comb[0])
        dictFam2 = readAlignment(sys.argv[1] + comb[1])
    # # dictFam1 = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
    # # dictFam2 = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
    #
    # print(dictFam1)
    # print(dictFam2)
        intersectionOfNames_tmp = intersectDictKeys(dictFam1, dictFam2)
        print(comb[0], comb[1])
        print(len(intersectionOfNames_tmp))
        if len(intersectionOfNames_tmp) > len(intersectionOfNames):
            # print(comb[0], comb[1], str(len(intersectionOfNames_tmp)))
            family_pair = [comb[0], comb[1]]
            intersectionOfNames = intersectionOfNames_tmp

    print(family_pair, len(intersectionOfNames))
    for speciesID in intersectionOfNames:
        print(speciesID + ": " + name(speciesID))
