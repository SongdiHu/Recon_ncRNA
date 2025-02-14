import sys
from urllib.request import urlopen
from urllib.error import HTTPError, URLError


def name(id):
    if '.' in id:
        id = id[:id.find('.')]
    url = 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % id
    try:
        f = urlopen(url)
    except HTTPError as e:
        print("HTTP Error:", e.code)
        print("PDB url %s: download failed" % url)
        sys.exit(1)
    except URLError as e:
        print("URL Error:", e.reason)
        print("PDB url %s: download failed" % url)
        sys.exit(1)
    for x in f:
        x = x.decode().strip()
        if '<title>' in x:
            return x[x.find('>')+1:-len(r' - Nucleotide - NCBI<\title>')]


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("USAGE : python ncby_name rfam_names.txt")
        sys.exit()

    rfam_names = open(sys.argv[1], 'r')
    ncby_write = open("ncby_names.txt", "w")
    for line in rfam_names:
        line = line.split(":")
        ncby_write.write(name(line[0])+'\n')
        # print(name(line[0]))
    ncby_write.close()
    # print(name('252393.30'))
