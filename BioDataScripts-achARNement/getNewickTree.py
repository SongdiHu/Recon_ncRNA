# This script takes a large newick tree and a list of species. The script outputs the subtree with that list of species only.

import sys
from Bio import Phylo
from io import StringIO
import re
import csv
from ncby_name import name

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("USAGE : python getNewickTree newickTree listOfSpeciesToKeep [remove plasmid, contig from name]")
        sys.exit()

    removePlasmid = False
    if len(sys.argv) == 4:
        removePlasmid = True

    with open('../CL00012/BVBRC_genome.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            print(row['Genome ID'], row['GenBank Accessions'])

    newickFile = open(sys.argv[1], 'r')
    for line in newickFile:
        # print("LINE", line)
        # line = re.sub("\(.*?:", "(" + name() + ":", line)
        # print(line)
        theTree = line  # there should be only one line

    phyloTree = Phylo.read(StringIO(theTree), "newick")

    print("Tree before pruning:")
    print(str(phyloTree))
    leaves = phyloTree.get_terminals(order='preorder')
    print("number of leaves:", len(leaves))
    print("number of leaves:", phyloTree.count_terminals())
    print("leaves:", leaves)

    names = lookup_by_names(phyloTree)
    for bv_id in names:
        # cld.name = name(bv_id)
        print(bv_id)
    ################################################ New Module ########################################################

    # phyloTree
    #
    # ####################################################################################################################
    #
    # listToKeep = []
    # listNodesToDel = []
    #
    # listToKeepFile = open(sys.argv[2], 'r')
    # dictNameToID = {}
    # for line in listToKeepFile:
    #
    #     if line == "\n":
    #         # print("empty line")
    #         continue
    #
    #     if len(line) < 6:  # the first line might contain the size of the intersection
    #         continue
    #
    #     colonSplit = line.split(": ")
    #
    #     name = colonSplit[0].split(",")[0]
    #
    #     if removePlasmid:  # we keep everything before the " cont" or " plasmid", etc.
    #         if " cont" in name:
    #             name = name.split(" cont")[0]
    #         if " plasm" in name:
    #             name = name.split(" plasm")[0]
    #         if ".cont" in name:
    #             name = name.split(".cont")[0]
    #         if " ctg" in name:
    #             name = name.split(" ctg")[0]
    #         if ".assembly" in name:
    #             name = name.split(".assembly")[0]
    #         if " gcont" in name:
    #             name = name.split(" gcont")[0]
    #         if " Cont" in name:
    #             name = name.split(" Cont")[0]
    #         if "_Cont" in name:
    #             name = name.split("_Cont")[0]
    #
    #     name = name.replace(" ", "_")
    #     name = name.strip()
    #     print("The name: " + name)
    #     dictNameToID[name] = colonSplit[0]
    #
    #     listToKeep.append(name)
    #
    # # listToKeep = ["Escherichia_coli_TA124", "Escherichia_coli_DEC2B"]
    #
    # print("List to keep size = " + str(len(listToKeep)))
    #
    # for element in phyloTree.find_elements():
    #     if element.name is not None and element.name not in listToKeep:
    #         # print("Adding node to listNodesToDel")
    #         listNodesToDel.append(element.name)
    #
    # for node in listNodesToDel:
    #     # print("Deleting node: " + node)
    #     try:
    #         phyloTree.prune(target=node)
    #     except:
    #         print("Could not delete node: " + node)
    #
    # print("Tree after pruning:")
    # print(str(phyloTree))
    #
    # outputTreeFilename = sys.argv[2].split(".")[0] + ".tree"
    #
    # outputNewick = open(outputTreeFilename, 'w')
    # Phylo.write(phyloTree, outputNewick, format="newick")
    #
    # outputNewick.write("\n\nMapping:\n")
    #
    # for key, value in dictNameToID.items():
    #     outputNewick.write(key + ":" + value + "\n")
