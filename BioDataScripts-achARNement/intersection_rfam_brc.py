import sys
import pandas as pd

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("USAGE: python intersection_rfam_brc rfam_names brc_genomes")
        sys.exit()

    setSeqsFile1 = set()
    setSeqsFile2 = set()

    rfam_file = open(sys.argv[1], 'r')
    brc_file = pd.read_csv(sys.argv[2])

    rfam_ids = []
    for line in rfam_file:
        rfam_ids.append(line[:line.find('.')])
    print(len(rfam_ids), rfam_ids)

    genbank_brc = brc_file['GenBank Accessions'].to_list()
    print(len(genbank_brc), genbank_brc)
    for ind in range(len(genbank_brc)):
        line = genbank_brc[ind]
        if line == 'CP003410.1':
            print("check")
        tmp_list = [line]
        if ',' in tmp_list[0]:
            tmp_list = tmp_list[0].split(',')
            # print(tmp_list)
        for ind_tmp in range(len(tmp_list)):
            item = tmp_list[ind_tmp]
            if '.' in item:
                tmp_list[ind_tmp] = item[:item.find('.')]
        genbank_brc[ind] = tmp_list[0]
        if len(tmp_list) > 1:
            # for item in tmp_list[1:]:
            #     genbank_brc.append(item)
            genbank_brc += tmp_list[1:]
    print(len(genbank_brc), genbank_brc)

    intersection = set(rfam_ids).intersection(set(genbank_brc))
    # # # for genBank_rfam in rfam_ids:
    # # #     if genBank_rfam in genbank_brc:
    # # #         intersection.append(genBank_rfam)
    # #
    print(len(intersection), intersection)
    print('JTHL01000001' in genbank_brc)

