#This script takes a fitchAndCompany output file (with the initial tree printed at the beginning) and outputs the tree in newick format, with names changed to something short
#and the sequences for each of the 2 families in separated fasta files (header will do the mapping)

######### Work in progress... unfinished, I did it by hand it was much easier :P

import sys
from Bio import Phylo
from io import StringIO



if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("USAGE : python readOutputFile-getTreeThatWasUsed  fitchAndCompany.out")
        sys.exit()

    nucleos = ['A','C','G','U']

    nodeIdToChildrenPair = {}  #dict that maps a node id to a pair of ids: (leftChild, rightChild)
    nodeIdToSeqPair = {}  #dict that maps a node id to a pair of sequences (only for leaves)
    nodeIdToNewName = {}  #dict that maps a leaf node to its new name
    nodeIdToPhyloClade = {}

    phyloTree = BaseTree.Tree()

    seqFam0Coming = False
    seqFam1Coming = False

    lastNodeId = 0
    seqFam0 = ""
    seqFam1 = ""
    familyId = -1

    #output
    prefix = sys.argv[1].split(".")[0]
    seqsFam0 = open(prefix + "_fam0.fa", 'w')
    seqsFam1 = open(prefix + "_fam1.fa", 'w')
    treeFile = open(prefix + ".theFinalTree", 'w')
        
    outputFile = open(sys.argv[1], 'r')
    for line in outputFile:
        if line[:4] == "--- ":  #we have an ID line
            twoSpacesSplit = line.split("  ")
            lastNodeId = twoSpacesSplit[0].split("Id = ")[1]
            lc = twoSpacesSplit[2].split[" = "][1]
            rc = twoSpacesSplit[3].split[" = "][1]

            if currentId in nodeIdToChildrenPair:  #We have read all the tree
                break

            nodeIdToChildrenPair[lastNodeId] = (lc, rc)  #left child id, followed by right child id

        elif line[0] == ".":
            if line[9] == '0':
                seqFam0 = ""
                seqFam1 = ""  #cleaning old sequences just to be sure
                familyId = 0
            elif line[9] == '1':
                familyId = 1

        elif line[0] == in nucleos:  #we are reading a sequence
            if familyId == -1:
                print("BIG PROBLEM: I didn't expect to have familyId == -1 here...")
                sys.exit(0)
            elif familyId == 0:
                seqFam0 = line.strip()
            elif familyId == 1:
                seqFam1 = line.strip()
                #We should have both seqs in memory here
                nodeIdToSeqPair[lastNodeId] = (seqFam0, seqFam1)
                #We know we have a leaf node, so we rename it
                nodeIdToNewName[lastNodeId] = "Species-" + str(lastNodeId)
        
    
    phyloTree = Phylo.read(StringIO(theTree), "newick")

    listToKeep = []
    listNodesToDel = []

    listToKeepFile = open(sys.argv[2], 'r')
    dictNameToID = {}
    for line in listToKeepFile:

        if line == "\n":
            #print("empty line")
            continue

        if len(line) < 6:   #the first line might contain the size of the intersection
            continue
        
        colonSplit = line.split(": ")

        name = colonSplit[1].split(",")[0]

        if removePlasmid:   #we keep everything before the " cont" or " plasmid", etc.
            if " cont" in name:
                name = name.split(" cont")[0]
            if " plasm" in name:
                name = name.split(" plasm")[0]
            if ".cont" in name:
                name = name.split(".cont")[0]
            if " ctg" in name:
                name = name.split(" ctg")[0]
            if ".assembly" in name:
                name = name.split(".assembly")[0]
            if " gcont" in name:
                name = name.split(" gcont")[0]
            if " Cont" in name:
                name = name.split(" Cont")[0]
            if "_Cont" in name:
                name = name.split("_Cont")[0]
        
        name = name.replace(" ", "_")
        name = name.strip()
        print("The name: " + name)
        dictNameToID[name] = colonSplit[0]
        
        listToKeep.append(name)

    #listToKeep = ["Escherichia_coli_TA124", "Escherichia_coli_DEC2B"]

    print("List to keep size = " + str(len(listToKeep)))
    
    for element in phyloTree.find_elements():
        if element.name is not None and element.name not in listToKeep:
            #print("Adding node to listNodesToDel")
            listNodesToDel.append(element.name)

    
    for node in listNodesToDel:
        #print("Deleting node: " + node)
        try:
            phyloTree.prune(target=node)
        except:
            print("Could not delete node: " + node)
            

    print("Tree after pruning:")
    print(str(phyloTree))

    outputTreeFilename = sys.argv[2].split(".")[0] + ".tree"
    
    outputNewick = open(outputTreeFilename, 'w')
    Phylo.write(phyloTree, outputNewick, format="newick")

    outputNewick.write("\n\nMapping:\n")

    for key, value in dictNameToID.items():
        
        outputNewick.write(key + ":" + value + "\n")
