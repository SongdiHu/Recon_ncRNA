import math
import os
import sys
import pickle
import time
import itertools
import RNA
from tqdm import tqdm
from io import StringIO
from collections import defaultdict
from collections import deque
from enum import Enum
from optparse import OptionParser
from Bio import Phylo
from structure_preparation_module import StructureMerge


#####################################################################################
# Ancestral RNA sequence reconstruction using sequence and 2D structure information
#
# Contains also the regular Fitch and Sankoff algorithms.
#####################################################################################


# A context manager for doing a "deep suppression" of stdout and stderr in
# Python, i.e. will suppress all print, even if the print originates in a
# compiled C/Fortran sub-function.
#    This will not suppress raised exceptions, since exceptions are printed
# to stderr just before a script exits, and after the context manager has
# exited (at least, I think that is why it lets exceptions through).
class suppress_stdout_stderr(object):

    def __init__(self):
        # Open a pair of null files
        # self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        self.null_fds = [os.open(os.devnull, os.O_RDWR)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        # self.save_fds = [os.dup(1), os.dup(2)]
        self.save_fds = [os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        # os.dup2(self.null_fds[0],1)
        # os.dup2(self.null_fds[1],2)
        os.dup2(self.null_fds[0], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        # os.dup2(self.save_fds[0],1)
        # os.dup2(self.save_fds[1],2)
        os.dup2(self.save_fds[0], 2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


###
### Class SpeciesTree: information on the whole species tree, methods to reconstruct ancestral sequences.
###
class SpeciesTree:
    ### Class variables
    # an Enum representing the different scoring functions for the ancestral sequence reconstruction
    ALGO = Enum('ALGO', ['Fitch', 'Sankoff', 'SeqAndStruct', 'TwoStructsApprox', 'StructsTreeDecomp'])
    NB_MUTATIONS = 1  # the maximum number of mutations to consider for the enumerate algorithm

    ###
    ### Static methods
    ###

    # To compare 2 trees: "simulTree" is the simulated tree and "inferredTree" is the tree containing the inferred
    # ancestral sequences. Obviously, the 2 trees are required to have the exact same topology and node ids. This
    # returns the total number of errors followed by the total number of nucleotides in the simulated tree (for all
    # sequences in the tree) in a list.
    def compareTrees(self, simulTree, inferredTree, verbose=False):

        if verbose:
            print("\nCOMPARING TREES")
        queue = deque()
        queue.appendleft(simulTree.root)
        queue.appendleft(inferredTree.root)  # Put the nodes from each tree one after the other in the same queue.

        ### THE STATS ###

        # All variables are lists of ints, representing, in this order: [family1, family2, all]

        # For the complete tree:
        totalNumberOfNucleos = [0, 0, 0]  # Will count the total number of nucleotides in the simulated tree (for all sequences of all families)
        totalNumberOfErrors = [0, 0, 0]  # Will count the total number of errors relative to the simulated sequences in all the tree

        totalNumberOfOptimalSeqs = [0, 0, 0]  # Will count the total number of optimal sequences for all the nodes in the species tree; this will allow us to compute the average nb of errors
        averageNumberOfErrors = [0, 0, 0]  # Will count the average number of errors over all the optimal solutions; will be divided by totalNumberOfOptimalSeqs

        # For the positions with/without structure
        totalNumberOfNucleos_withStruct = [0, 0, 0]  # Will count the total number of nucleotides that are part of a structure
        totalNumberOfErrors_withStruct = [0, 0, 0]  # Will count the total number of errors relative to the simulated sequences in all the tree at structured positions

        totalNumberOfNucleos_withoutStruct = [0, 0, 0]  # Will count the total number of nucleotides that are not part of a structure
        totalNumberOfErrors_withoutStruct = [0, 0, 0]  # Will count the total number of errors relative to the simulated sequences in all the tree at unstructured positions

        averageNumberOfErrors_withStruct = [0, 0, 0]  # Will count the average number of errors over all the optimal solutions in all structured positions; divided by totalNumberOfOptimalSeqs
        averageNumberOfErrors_withoutStruct = [0, 0, 0]  # Will count the average number of errors over all the optimal solutions in all unstructured positions; divided by totalNumberOfOptimalSeqs
        numberOfNodesInTheTree = [0, 0, 0]  # Will count the number of nodes for family 1, family 2, and the 2 together (result should be [x, x, 2x])

        # For the root only
        totalNumberOfErrors_root = [0, 0, 0]
        totalNumberOfOptimalSeqs_root = [0, 0, 0]  # With the new findSuperAncestorAndUpdateMatrix method, we should have the same number for fam 0 and fam 1
        averageNumberOfErrors_root = [0, 0, 0]
        averageNumberOfErrors_root_withoutStruct = [0, 0, 0]
        averageNumberOfErrors_root_withStruct = [0, 0, 0]

        #################

        while queue:  # while the queue is not empty
            theNodeSimul = queue.pop()
            theNodeInfer = queue.pop()

            # Making sure the ids are identical
            if theNodeSimul.id != theNodeInfer.id:
                print("WARNING: Node ids are different. I did not expect that...")

            if verbose:
                print("- Node: " + str(theNodeSimul.id))

            # Comparing the sequences
            for familyIndex in range(simulTree.getNbStructs()):
                # simulSeq = theNodeSimul.getSequence(familyIndex, 0)  #there is only one sequence for each family in the simulated tree
                simulSeq = theNodeSimul.getOnlySequence(familyIndex)  # NEW: updated to deal with sets
                # test
                # print("The only sequence:", simulSeq)

                minNbErrors = float("inf")
                minNbErrorsWithStruct = float("inf")
                minNbErrorsWithoutStruct = float("inf")

                totalNumberOfOptimalSeqs[familyIndex] += theNodeInfer.getNbSeqs(
                    familyIndex)  # adding the number of seqs for the current familyIndex
                totalNumberOfOptimalSeqs[2] += theNodeInfer.getNbSeqs(
                    familyIndex)  # adding the number of seqs for all (family1 and family2)
                numberOfNodesInTheTree[familyIndex] += 1
                numberOfNodesInTheTree[2] += 1

                # Counting the number of positions with structure and without structure (only once for each family for each node)
                for pos in range(len(simulSeq)):
                    pairedPos = theNodeInfer.getPairedPos(familyIndex, pos)  # theNodeInfer and theNodeSimul should have the same structs and same seq length
                    if pairedPos == -1:  # unpaired
                        totalNumberOfNucleos_withoutStruct[
                            familyIndex] += 1  # Those numbers will have to be divided by the totalNumberOfOptimalSeqs
                        totalNumberOfNucleos_withoutStruct[2] += 1

                    else:  # paired
                        totalNumberOfNucleos_withStruct[familyIndex] += 1
                        totalNumberOfNucleos_withStruct[2] += 1

                for seq in theNodeInfer.getSequencesListsSet(familyIndex):
                    if len(simulSeq) != len(seq):
                        print("Big problem: sequences are not the same length...")

                    else:
                        nbErrors = 0
                        nbErrorsWithStruct = 0
                        nbErrorsWithoutStruct = 0
                        for pos in range(len(seq)):
                            pairedPos = theNodeInfer.getPairedPos(familyIndex, pos)

                            if SpeciesNode.nucToNucDict[simulSeq[pos]] != SpeciesNode.nucToNucDict[seq[pos]]:
                                # Use the translation, just in case.
                                nbErrors += 1
                                averageNumberOfErrors[
                                    familyIndex] += 1  # Add all the errors, we'll divide it by totalNumberOfOptimalSeqs at the end.
                                averageNumberOfErrors[2] += 1  # for all

                                if pairedPos == -1:  # unpaired
                                    nbErrorsWithoutStruct += 1
                                    averageNumberOfErrors_withoutStruct[familyIndex] += 1
                                    averageNumberOfErrors_withoutStruct[2] += 1

                                else:  # paired
                                    nbErrorsWithStruct += 1
                                    averageNumberOfErrors_withStruct[familyIndex] += 1
                                    averageNumberOfErrors_withStruct[2] += 1

                        if nbErrors < minNbErrors:
                            minNbErrors = nbErrors
                            minNbErrorsWithStruct = nbErrorsWithStruct
                            minNbErrorsWithoutStruct = nbErrorsWithoutStruct

                if verbose:
                    # Printing a report
                    print("Family : " + str(familyIndex) + " --> minNbErrors = " + str(minNbErrors) + " (" + str(
                        minNbErrors / len(simulSeq) * 100) + " %)")

                # updating stats for the whole tree
                totalNumberOfErrors[familyIndex] += minNbErrors
                totalNumberOfErrors[2] += minNbErrors  # updating all in the list
                totalNumberOfNucleos[familyIndex] += len(simulSeq)
                totalNumberOfNucleos[2] += len(simulSeq)

                totalNumberOfErrors_withoutStruct[familyIndex] += minNbErrorsWithoutStruct
                totalNumberOfErrors_withoutStruct[2] += minNbErrorsWithoutStruct

                totalNumberOfErrors_withStruct[familyIndex] += minNbErrorsWithStruct
                totalNumberOfErrors_withStruct[2] += minNbErrorsWithStruct

            if not theNodeSimul.isLeaf():  # only checking one, both of them should be leaves at the same time anyway
                queue.appendleft(theNodeSimul.leftChild)
                queue.appendleft(theNodeInfer.leftChild)
                queue.appendleft(theNodeSimul.rightChild)
                queue.appendleft(theNodeInfer.rightChild)

                if theNodeInfer.isRoot():  # We save the stats for the root only
                    for index in range(len(totalNumberOfErrors)):
                        totalNumberOfErrors_root[index] = totalNumberOfErrors[index]
                        totalNumberOfOptimalSeqs_root[index] = totalNumberOfOptimalSeqs[index]
                        averageNumberOfErrors_root[index] = averageNumberOfErrors[index] / \
                                                            totalNumberOfOptimalSeqs_root[
                                                                index]  # We calculate the averages immediately
                        averageNumberOfErrors_root_withoutStruct[index] = averageNumberOfErrors_withoutStruct[index] / \
                                                                          totalNumberOfOptimalSeqs_root[index]
                        averageNumberOfErrors_root_withStruct[index] = averageNumberOfErrors_withStruct[index] / \
                                                                       totalNumberOfOptimalSeqs_root[index]

        for x in range(3):  # Calculating the averages -> numbers will be for one sequence only. (average nb errors for 1 optimal sequence in the tree)
            averageNumberOfErrors[x] /= totalNumberOfOptimalSeqs[x]
            averageNumberOfErrors_withStruct[x] /= totalNumberOfOptimalSeqs[x]
            averageNumberOfErrors_withoutStruct[x] /= totalNumberOfOptimalSeqs[x]

        print("totalNumberOfNucleos_withStruct:", totalNumberOfNucleos_withStruct)
        return totalNumberOfErrors, totalNumberOfNucleos, totalNumberOfOptimalSeqs, averageNumberOfErrors, totalNumberOfNucleos_withStruct, totalNumberOfErrors_withStruct, totalNumberOfNucleos_withoutStruct, totalNumberOfErrors_withoutStruct, averageNumberOfErrors_withStruct, averageNumberOfErrors_withoutStruct, numberOfNodesInTheTree, totalNumberOfErrors_root, totalNumberOfOptimalSeqs_root, averageNumberOfErrors_root, averageNumberOfErrors_root_withoutStruct, averageNumberOfErrors_root_withStruct

    # To translate a parenthesis representation of a structure into a list of paired positions
    def parenthesisToPairedPos(self, parenthesisString):
        pairedPosList = [-1] * len(parenthesisString)  # initializing the list with -1 everywhere

        openParenthesisStack = deque()  # stack that will be used to put the positions of open parentheses, will be poped when we read a closing parenthesis.
        for pos in range(len(parenthesisString)):

            char = parenthesisString[pos]

            if char == '(':
                openParenthesisStack.append(pos)

            elif char == ')':
                posPair = openParenthesisStack.pop()
                pairedPosList[pos] = posPair
                pairedPosList[posPair] = pos

        return pairedPosList

    # Constructor
    def __init__(self, structuresList_dot_bracket, structure_w, bags_in_order, root=None):
        ### Data attributes
        self.structuresList_dot_bracket = structuresList_dot_bracket
        self.nbStructs = len(structuresList_dot_bracket)  # the number of different RNA families
        self.structuresList = []  # A list a all the structures (one per family). A structure is represented by a list(int) of size sequenceLength.
        self.sequenceLength = len(structuresList_dot_bracket[0])  # The length of every sequence. (right now we have the same length for family 1 and family 2)
        self.structure_weight = structure_w
        self.bags_in_order = bags_in_order
        self.root = root  # The root of the tree (SpeciesNode)
        self.maxDepth = -1  # the maximum depth of the tree, will be updated after the call to setDepths

        for structure in structuresList_dot_bracket:
            # print(structure)
            self.structuresList.append(self.parenthesisToPairedPos(structure))
        # print(self.structuresList)

    ###
    ### Accessors
    ###

    def getStructs_dot_bracket(self):
        return self.structuresList_dot_bracket

    def getNbStructs(self):
        return self.nbStructs

    def getSeqLength(self):
        return self.sequenceLength

    def getRoot(self):
        return self.root

    def isRoot(self, node):
        if node == self.root:
            return True

        return False

    def isLeaf(self, node):
        if node.leftChild is None and node.rightChild is None:
            return True

        return False

    # Returns the list of paired positions for familyIndex
    def getStruct(self, familyIndex):
        return self.structuresList[familyIndex]

    # To print the whole tree in level order.
    # printCosts is a boolean indicating if we want to print the cost of every nucleotide for every position.
    def printTree(self, printCosts=False):
        nbOptSeqs = 0
        # test
        self.root.printNode(printCosts)
        # queue = deque()
        # queue.append(self.root)
        #
        #
        # while queue:  # while the queue is not empty
        #
        #     theNode = queue.pop()
        #     nbOptSeqs += theNode.printNode(printCosts)
        #     if not theNode.isLeaf():
        #         queue.appendleft(theNode.leftChild)
        #         queue.appendleft(theNode.rightChild)

        return nbOptSeqs

    # To print everything below the given node (including the node)
    def printSubtree(self, node):
        if node is None:
            return

        node.printNode()
        self.printSubtree(node.leftChild)
        self.printSubtree(node.rightChild)

    # To create a clone of the tree for which we remove the sequences at the ancestral nodes (to create a copy of the
    # simulated tree, without the answers (ie anc. sequences)).
    def cloneWithoutAncSeqs(self):
        treeClone = SpeciesTree(self.structuresList_dot_bracket, self.structure_weight, self.bags_in_order)  # Creating the new tree with the same structures and root -> not using a copy of structuresList, but shouldn't be a problem.
        treeClone.root = self.root.clone(treeClone)  # setting the root

        self.cloneWithoutAncSeqsRec(self.root, treeClone.root, treeClone)  # calling the recursive method

        return treeClone

    # Recursive method to clone the tree
    # internalNode and internalNodeClone correspond to the same node
    def cloneWithoutAncSeqsRec(self, internalNode, internalNodeClone, treeClone):
        if internalNode.isLeaf():
            return

        leftChildClone = internalNode.leftChild.clone(treeClone)
        rightChildClone = internalNode.rightChild.clone(treeClone)

        internalNodeClone.leftChild = leftChildClone
        internalNodeClone.rightChild = rightChildClone

        leftChildClone.parent = internalNodeClone
        rightChildClone.parent = internalNodeClone

        # Recursive calls
        self.cloneWithoutAncSeqsRec(internalNode.leftChild, leftChildClone, treeClone)
        self.cloneWithoutAncSeqsRec(internalNode.rightChild, rightChildClone, treeClone)

    ###
    ### Mutators
    ###

    # Method to set the depth of every node in the tree (root has depth 0)
    def setDepths(self):
        depth = 0
        self.setDepthsRec(self.root, depth)  # call to the recursive method

    # Recursive method to set the depths of the nodes.
    # Parameter 'node' is the node to update and 'depth' is the depth it has.
    def setDepthsRec(self, node, depth):
        node.depth = depth
        if depth > self.maxDepth:
            self.maxDepth = depth
        if not node.isLeaf():
            self.setDepthsRec(node.leftChild, depth + 1)
            self.setDepthsRec(node.rightChild, depth + 1)

    # General method to do the reconstruction; calls the bottom-up step and then the top-down step
    def reconstructAncestralSeqs(self, algo, p_bar=None):
        # print("Calculating costs...")
        # st_time = time.time()
        self.calculateCosts(algo, p_bar)
        new_st = time.time()
        # print("Costs Calculated in", new_st - st_time, ".")
        self.selectNucleos(algo, p_bar)
        print("Nucleos selected in", time.time() - new_st, ".")

    # The top-down step for the Enumerate algorithm. optimalMutSeq is the MutatedSequence object that is on the
    # optimal path from the root to the 'node' parameter.
    def selectMutatedSequences(self, familyIndex, node, optimalMutSeq):
        if not node.isLeaf():  # otherwise, there is nothing to do

            node.addSequence(familyIndex,
                             optimalMutSeq.sequence)  # adding the sequence as an optimal ancestral sequence

            self.selectMutatedSequences(familyIndex, node.leftChild, optimalMutSeq.leftChild)
            self.selectMutatedSequences(familyIndex, node.rightChild, optimalMutSeq.rightChild)  # recursive calls

    # The bottom-up step, calculating the costs using the algo specified by the parameter 'algo'.
    def calculateCosts(self, algo, c_bar):
        self.root.calculateCosts(algo, c_bar)

    # The top-down step, called after calculateCosts; this method selects the nucleotide of minimum cost for each
    # position and builds the ancestral sequences.
    def selectNucleos(self, algo, s_bar):
        self.findSuperAncestorAndUpdateMatrix(algo)

        # # test
        # print("Super anc found.")
        numberOfOptimalSeqsAtRoot = 1  # Counts the total number of opt. seqs. at the root. (for one family only, since it's the same number now because of Fitch at the root)
        listOfListOfOptNucleos = []  # For each position in the ancestral seq, a list of all optimal nucleos.
        maxPossibleNucleo = 0  # Stores the largest number of optimal nucleos that we can find at one position in the ancestral sequence.
        keepTrying = True  # Indicates if we have to stop trying to enumerate all optimal sequences.

        if algo != SpeciesTree.ALGO(5):
            # We start by choosing the nucleotides of minimum cost in the root. WARNING: For now, we only build 1 of all
            # the possible optimal sequences
            for familyIndex in range(self.getNbStructs()):
                # One optimal sequence
                sequenceList = [[" " for x in range(
                    self.getSeqLength())]]  # A list of optimal sequences (each sequence is a list of chars here)

                for pos in range(self.getSeqLength()):
                    # print(pos)
                    # print(len(sequenceList))
                    minCost = float("inf")
                    nucleoOfMinCost = []  # a list of all optimal nucleos for the position

                    for key in self.root.posToNucleoCosts[familyIndex][pos]:  # Finding all nucleotides of minimum cost

                        cost = self.root.getCostOf(familyIndex, pos, key)
                        if cost == minCost:
                            # print("appending " + key)
                            nucleoOfMinCost.append(key)  # appending to the list
                        if cost < minCost:
                            # print("Updating minCost: minCost = " + str(cost) + " and nucleo = " + nucleo)
                            minCost = cost
                            nucleoOfMinCost = [key]  # starting a new list

                    if len(nucleoOfMinCost[0]) == 1:
                        if keepTrying:
                            keepTrying = self.root.addNucleoToSequencesList(sequenceList, pos, nucleoOfMinCost)

                        # Only need to do this once (the 2 families have the same ancestral sequences at the root now
                        # because of Fitch step.
                        if familyIndex == 0:
                            listOfListOfOptNucleos.append(nucleoOfMinCost)
                            numberOfOptimalSeqsAtRoot *= len(nucleoOfMinCost)
                            if len(nucleoOfMinCost) > maxPossibleNucleo:
                                maxPossibleNucleo = len(nucleoOfMinCost)

                    # May 27, 2016: I believe this next condition will never happen anymore because of the Fitch that
                    # we do at the root (findSuperAncestralSequences)
                    elif len(nucleoOfMinCost[0]) == 2:  # base pair position
                        posPair = self.root.getPairedPos(familyIndex, pos)
                        if pos < posPair:  # otherwise, it has been updated already
                            self.root.addNucleoToSequencesList(sequenceList, pos, [nuc[0] for nuc in nucleoOfMinCost])
                            self.root.addNucleoToSequencesList(sequenceList, posPair, [nuc[1] for nuc in nucleoOfMinCost])

                if keepTrying:
                    for seq in sequenceList:  # adding all the sequences
                        # print(seq)
                        self.root.addSequence(familyIndex, "".join(seq))

        else:
            for familyIndex in range(self.getNbStructs()):
                sequence_template = ['_'] * self.sequenceLength
                optim_seqs_dict = [math.inf, []]
                # # test
                # st_time = time.time()
                self.root.selectNucleos4Bag_rec(familyIndex, 0, sequence_template, 0, optim_seqs_dict, order=0)
                self.root.selectNucleos4Bag_rec(familyIndex, (len(self.bags_in_order)-1), sequence_template, 0, optim_seqs_dict, order=1)
                sequenceList = optim_seqs_dict[1]
                # # test
                # print("Number of potential sequences:", len(sequenceList))
                # print("Optimal sequences done. (one family)", time.time() - st_time)
                numberOfOptimalSeqsAtRoot = len(sequenceList)
                if numberOfOptimalSeqsAtRoot > 25000000:
                    print("Too many sequences to enumerate.(root)")
                    sys.exit()
                if keepTrying:
                    for seq in sequenceList:  # adding all the sequences
                        self.root.addSequence(familyIndex, "".join(seq))

        if s_bar is not None:
            s_bar.update(1)

        # To print the possible nucleotides at each position and the size of the optimal ancestral sequences at the
        # root.
        ancRepresentationString = "\nNumber of optimal sequences at the root = " + str(
            int(numberOfOptimalSeqsAtRoot)) + "\n"

        if algo != SpeciesTree.ALGO(5):
            for x in range(maxPossibleNucleo):
                for y in range(len(listOfListOfOptNucleos)):
                    if len(listOfListOfOptNucleos[y]) > x:
                        ancRepresentationString += listOfListOfOptNucleos[y][x]
                    else:
                        ancRepresentationString += " "
                ancRepresentationString += "\n"

        print(ancRepresentationString)

        if not options.rootOnly:
            # Now that the root sequences have been found, we call the SpeciesNode selectNucleos method on the root
            # to select the nucleos for its two children nodes.
            # # test
            # print("Start selecting children.")
            # st_time = time.time()
            self.root.selectNucleos(algo, s_bar)
            # print("Finished selecting children.", time.time() - st_time)


    # This is non-recursive.
    # Nucleotides that gives the minimum cost for each bag will be selected.
    def selectNucleos4Bag(self, fam, sequence_template_orig):
        sequence_template = sequence_template_orig
        sequence_template = list(sequence_template)
        sequence_template = [[x] for x in sequence_template]
        optim_bags = []

        for bag_index in range(len(self.bags_in_order)):

            cost_nucleos = self.root.bagsCostMatrix[fam][bag_index]
            bag = self.bags_in_order[bag_index]
            args_bag = [['A', 'C', 'G', 'U'] for x in range(len(bag))]
            min_cost = float('inf')
            min_cost_list = []  # A list containing all combinations for a bag that gives the least cost.

            # Find positions in the sequence where values are assigned in previous bag(s).
            for ind in range(len(bag)):
                nucleotides = sequence_template[bag[ind]]
                if nucleotides[0] != '_':
                    args_bag[ind] = nucleotides

            for combination in itertools.product(*args_bag):
                cost_tmp = cost_nucleos[combination]
                if cost_tmp < min_cost:
                    min_cost_list = [combination]
                    min_cost = cost_tmp
                elif cost_tmp == min_cost:
                    min_cost_list.append(combination)

            # TODO: 'min_cost_list' # of branches?

            for ind in range(len(bag)):
                for comb in min_cost_list:
                    if comb[ind] in sequence_template[bag[ind]]:
                        pass
                    elif sequence_template[bag[ind]][0] != '_':
                        sequence_template[bag[ind]].append(comb[ind])
                    else:
                        sequence_template[bag[ind]] = [comb[ind]]

            optim_bags.append(min_cost_list)

        # print("SEQ_TEMPLATE:", sequence_template)
        # print("OPTIM_BAGS:", optim_bags)
        return optim_bags

    # To calculate the matrix costs, using Fitch for the super ancestor (ancestor of family 1 and family 2 of the
    # root). Then, this method will update the costMatrix of the root with this super ancestor matrix for both
    # families (same matrix copied for both)
    def findSuperAncestorAndUpdateMatrix(self, algo):
        # First, we need to create a fake tree of 3 nodes, so that we can call computeFitch normally. We don't need
        # it. We need only one list, so that it will be like there is only one family. The length of the
        # structureList is important, as it defines the length of the sequences in the tree.
        fakeSpeciesTree = SpeciesTree(self.getStructs_dot_bracket(), self.structure_weight, self.bags_in_order)
        superAncNode = SpeciesNode(None, None, None, self.structure_weight, self.bags_in_order,
                                   fakeSpeciesTree)  # corresponds to the super ancestor
        fakeSpeciesTree.root = superAncNode
        leftChild = SpeciesNode(superAncNode, None, None, self.structure_weight, self.bags_in_order,
                                fakeSpeciesTree)  # corresponds to Family 0 of the root
        rightChild = SpeciesNode(superAncNode, None, None, self.structure_weight, self.bags_in_order,
                                 fakeSpeciesTree)  # corresponds to Family 1 of the root
        superAncNode.leftChild = leftChild
        superAncNode.rightChild = rightChild

        # We need to check if we are dealing with di-nucleotides (which is the case for algos 3 and 4)
        if algo != SpeciesTree.ALGO(5):
            if algo == SpeciesTree.ALGO(3) or algo == SpeciesTree.ALGO(4):
                for familyIndex in range(self.root.getNbStructs()):
                    for pos in range(self.getSeqLength()):
                        posPair = self.root.getPairedPos(familyIndex, pos)

                        if posPair != -1:  # we are dealing with di-nucleotides in the costMatrix

                            # Stores the lowest cost for A,C,G and U over all the di-nucleotide score. (looking only
                            # at first position)
                            dictBestScoreForFirstNucleo = {}
                            for dinucleo in SpeciesNode.DINUCLEOS:

                                dinucleoCost = self.root.getCostOf(familyIndex, pos, dinucleo)

                                # first nucleo of the dinucleo is the one corresponding to the current position.
                                correspondingNucleo = dinucleo[0]

                                if correspondingNucleo not in dictBestScoreForFirstNucleo or dinucleoCost < \
                                        dictBestScoreForFirstNucleo[correspondingNucleo]:
                                    dictBestScoreForFirstNucleo[correspondingNucleo] = dinucleoCost

                            # We update the costMatrix, adding new key-values for single nucleotides.
                            for nucleo, minCost in dictBestScoreForFirstNucleo.items():
                                self.root.updatePosToNucleoCosts(familyIndex, pos, nucleo, minCost)

            # Here, we map the costMatrix of family 0 to the leftChild of the fakeTree and the costMatrix of family 1
            # to the rightChild of the fakeTree.
            leftChild.posToNucleoCosts[0] = self.root.posToNucleoCosts[0]
            rightChild.posToNucleoCosts[0] = self.root.posToNucleoCosts[1]

            for pos in range(fakeSpeciesTree.getSeqLength()):  # for each position
                # calling computeFitch for every position, in family 0 only (there is only one family in the fakeTree)
                superAncNode.computeFitch(0, pos)

            # Now the tricky part: we replace the cost matrices of the root (for fam 0 and fam 1) with the matrix of
            # the superAncestor... and we hope for the best.
            self.root.posToNucleoCosts[0] = superAncNode.posToNucleoCosts[0]
            self.root.posToNucleoCosts[1] = superAncNode.posToNucleoCosts[0]

        # For tree-decomposition approach, cost matrices are different: [fam][bag][nucleos -- chars in tuple]
        else:
            leftChild.bagsCostMatrix[0] = self.root.bagsCostMatrix[
                0]  # Here, we map the costMatrix of family 0 to the leftChild of the fakeTree
            rightChild.bagsCostMatrix[0] = self.root.bagsCostMatrix[
                1]  # and the costMatrix of family 1 to the rightChild of the fakeTree

            # TODO: merging of the two cost matrices. Do we need a different Fitch function?
            for bagIndex in range(len(self.bags_in_order)):  # for each position
                # calling computeFitch for every position, in family 0 only (there is only one family in the fakeTree)
                superAncNode.computeTreeDecomp(0, bagIndex, True)

            self.root.bagsCostMatrix[0] = superAncNode.bagsCostMatrix[0]
            self.root.bagsCostMatrix[1] = superAncNode.bagsCostMatrix[0]

    def convertBagsToSequences(self, family_id, cost_matrix):
        # The list contains at each bag, the nucleos with the minimum cost. (2-dimensional)
        optimal_sequences_unprocessed = []
        # The list contains all optimal sequences.
        optinal_sequences_processed = []
        self.getNbStructs()

        # Create a list containing all nucleo-combinations that gives the minimal cost for each bag.
        for bagIndex in range(len(self.bags_in_order)):  # for each position
            cost_nucleos = cost_matrix[family_id][bagIndex]
            optimal_nucleos = min(cost_nucleos, key=cost_nucleos.get)
            optimal_sequences_unprocessed.append(optimal_nucleos)

        sequence_template = ['_'] * self.sequenceLength
        assigned_positions = set()

        # Loop through all bags.
        for bagIndex in range(len(self.bags_in_order)):
            positions = self.bags_in_order[bagIndex]
            optimal_nucleos = optimal_sequences_unprocessed[bagIndex]
            sequence_template_tmp = sequence_template

            # Check intersection.
            intersection_positions = assigned_positions.intersection(set(positions))
            for combination in optimal_nucleos:
                sequence_template_tmp_tmp = sequence_template_tmp
                # first check if the intersection values mismatch.
                for index in range(len(combination)):
                    position = positions[index]

                    # Check intersection.
                    if position not in intersection_positions:
                        sequence_template_tmp_tmp[position] = combination[index]
                    elif combination[index] == sequence_template_tmp_tmp[position]:
                        pass
                    else:
                        break

                    # Commit the changes when reach the end.
                    if index == len(combination) - 1:
                        sequence_template_tmp = sequence_template_tmp_tmp


###
### Class SpeciesNode: information on a node of the SpeciesTree and methods to calculate scores
###
class SpeciesNode:
    BP_WEIGHT = 0.001  # The weight that will be used to multiply the values in the BASEPAIR_MATRIX_SIMPLE  --> the one we ended up using for ISMB
    # BP_WEIGHT = 0.05  #That was the best so far  (61 and 59 total error for 3 and 4 resp.)
    # BP_WEIGHT = 0.5  #A test for bio data --> does not seem to work with FinP-Traj
    NUCLEOS = ['A', 'C', 'G', 'U']  # Just a list of the nucleotides, will help with for loops
    DINUCLEOS = ['AA', 'AC', 'AG', 'AU', 'CA', 'CC', 'CG', 'CU', 'GA', 'GC', 'GG', 'GU', 'UA', 'UC', 'UG', 'UU']
    PAIRABLES = {'A': ['U'], 'U': ['A', 'G'], 'C': ['G'], 'G': ['C', 'U']}
    FITCH_MATRIX = [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]]
    SANKOFF_MATRIX = [[0, 2, 1, 2], [2, 0, 2, 1], [1, 2, 0, 2], [2, 1, 2, 0]]
    # BASEPAIR_MATRIX_SIMPLE = [[5, 5, 5, 1], [5, 5, 0, 5], [5, 0, 4, 3], [1, 5, 4, 5]]  #not bad, but 3 does better than 4
    BASEPAIR_MATRIX_SIMPLE = [[3, 3, 3, 1], [3, 3, 0, 3], [3, 0, 3, 2], [1, 3, 2, 3]]  # this one is not bad at all
    # BASEPAIR_MATRIX_SIMPLE = [[2.5, 2.5, 2.5, 0], [2.5, 2.5, 0, 2.5], [2.5, 0, 2.5, 2.5], [0, 2.5, 2.5, 2.5]]  #first matrix
    for i in range(len(BASEPAIR_MATRIX_SIMPLE)):
        for j in range(len(BASEPAIR_MATRIX_SIMPLE[i])):
            BASEPAIR_MATRIX_SIMPLE[i][j] *= BP_WEIGHT
    LIST_MUT_SEQS_MAX_SIZE = 1000  # maximum size of listMutatedSequences -> we want to order the list and keep the 1000 with lowest cost (most promising)

    _NODE_ID = 1  # Class variable; incrementing every time we call the constructor
    nucleoToIndexDict = {'A': 0, 'a': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3, 'U': 3,
                         'u': 3}  # mapping a nucleotide to an index, starting from 0
    nucToNucDict = {'A': 'A', 'a': 'A', 'C': 'C', 'c': 'C', 'G': 'G', 'g': 'G', 'T': 'U', 't': 'U', 'U': 'U',
                    'u': 'U'}  # translating nucleotides so that we deal only with A,C,G,U

    def addNucleoToSequencesList(self, sequencesList, position, listOfNucsToAdd):
        # New: to protect from memory explosion
        if len(sequencesList) > 25000000:
            print("Too many optimal sequences.", len(sequencesList), "position:", position)
            if options.rootOnly:
                return False  ####returns False if there are too many possible sequences and we have to stop trying
            sys.exit()

        setOfNucsToAdd = set()
        for nuc in listOfNucsToAdd:  # getting rid of duplicates (possibly coming from di-nucleotides) with a set
            setOfNucsToAdd.add(nuc)

        listOfNucsToAdd = list(setOfNucsToAdd)  # rewriting the list
        # print("--- after elimination of duplicates, listOfNucsToAdd = " + str(listOfNucsToAdd))

        initialNbSeqs = len(sequencesList)
        # No need to duplicate sequences if there is only one nucleo to add,
        # otherwise we need one additional copy for every additional nucleo to add.
        for x in range(len(setOfNucsToAdd) - 1):
            for y in range(initialNbSeqs):
                copy = sequencesList[y][:]
                sequencesList.append(copy)

        for x in range(len(sequencesList)):
            nucIndex = x // initialNbSeqs  # integer division to find the index of the nucleo to add in the sequences
            sequencesList[x][position] = listOfNucsToAdd[nucIndex]

        # print("At the end, sequencesList = " + str(sequencesList))
        return True  # if we get to this point, we have to return True now (May 30, 2016)

    def addOptimBagsToSequencesList(self, optim_bags_list):
        # all_fixed_bags_list = []  # list(list(tuple())) list of tuple lists, each represents an optimal sequence.
        sequences_list = []
        # Go through all possible combinations of nucleos in the bag.
        # keep the ones with the minimum cost.
        # # test
        # print("optim_bags_list length:", len(optim_bags_list), type(optim_bags_list[0][0]))
        args_bag = optim_bags_list.copy()
        # TODO: Create a list that contains number of bag combinations to skip at each bag.
        skip = 0
        tmp_skip = 1
        skip_list = list()
        optim_bags_list.reverse()
        # print("++++++++++++++++++++++++++++", optim_bags_list)
        # print("----------------------------", args_bag)

        for combs in optim_bags_list:
            tmp_skip *= len(list(combs))
            skip_list.append(tmp_skip)

        skip_list.reverse()

        for td_bags in itertools.product(*args_bag):
            # print(td_bags)
            if skip > 0:
                skip -= 1
                continue

            sequence_template = self.sequenceTemplate.copy()
            next_bag = True
            for bagID in range(len(td_bags)):

                positions = self.bags_ordered_list[bagID]
                fixed_nucleos = td_bags[bagID]
                for ind in range(len(fixed_nucleos)):
                    if sequence_template[positions[ind]] == '$':
                        sequence_template[positions[ind]] = fixed_nucleos[ind]
                    elif sequence_template[positions[ind]] != fixed_nucleos[ind]:
                        # # test
                        # print(ind, bagID, sequence_template[positions[ind]], fixed_nucleos[ind])
                        if bagID < (len(td_bags) - 1):
                            skip = skip_list[bagID + 1] - 1
                        next_bag = False
                        break
                if not next_bag:
                    break
            if '$' not in sequence_template:
                sequences_list.append(sequence_template.copy())

        # No valid sequence found, try the sub-optimal one.
        if not sequences_list:
            fixed_bags = [i[0] for i in optim_bags_list]
            sequence_template = self.sequenceTemplate.copy()
            for bagID in range(len(fixed_bags)):
                positions = self.bags_ordered_list[bagID]
                fixed_nucleos = fixed_bags[bagID]
                for ind in range(len(fixed_nucleos)):
                    if sequence_template[positions[ind]] == '$':
                        sequence_template[positions[ind]] = fixed_nucleos[ind]
            if '$' not in sequence_template:
                sequences_list.append(sequence_template.copy())

        return sequences_list

    # Constructor
    def __init__(self, parent, leftChild, rightChild, structure_w, bags_ordered_list, speciesTree, id=None):

        ### Data attributes
        self.structure_weight = structure_w
        self.parent = parent  # the parent node in the tree (root has parent = None)
        self.leftChild = leftChild  # the leftChild in the tree (leaves have leftChild = None)
        self.rightChild = rightChild  # the rightChild in the tree (leaves have rightChild = None)
        self.speciesTree = speciesTree  # a reference to the SpeciesTree that contains this node
        self.sequenceLength = speciesTree.getSeqLength()
        self.bags_ordered_list = bags_ordered_list  # Ordered bags from the tree-decomposition.
        self.bagsCostMatrix = [[defaultdict(lambda: float("inf")) for i in range(len(bags_ordered_list))]
                               for j in range(self.getNbStructs())] if bags_ordered_list else None
        # For each family (nbStructs), we have one list of all positions, and at each position a dict mapping a
        # nucleotide (or a dinucleotide) to its cost.
        self.posToNucleoCosts = [[defaultdict(lambda: float("inf")) for i in range(self.sequenceLength)] for j in
                                 range(self.getNbStructs())]
        # One *set* (*new, before it was a list) of optimal sequences (strings) for every family
        self.sequencesLists = [set() for i in range(speciesTree.getNbStructs())]
        self.sequenceTemplate = ["$"] * self.sequenceLength
        # self.structureTemplate = ["."] * self.sequenceLength

        if id is None:
            self.id = self.__class__._NODE_ID
            self.__class__._NODE_ID += 1
        else:
            self.id = id
            self.__class__._NODE_ID = id + 1

        self.depth = -1  # intializing to -1, don't forget to call setDepths() on the tree after it is contructed
        # list of most promising MutatedSequence objects (only one sequence for leaves) for each family
        self.listMutatedSequences = [[] for i in range(speciesTree.getNbStructs())]

    ###
    ### Accessors
    ###

    def isRoot(self):
        return self.speciesTree.isRoot(self)

    def isLeaf(self):
        return self.speciesTree.isLeaf(self)

    def getSeqLength(self):
        return self.speciesTree.getSeqLength()  # Warning: right now, we assume that all sequences have the same length

    def getNbStructs(self):
        return self.speciesTree.getNbStructs()

    def getPairedPos(self, familyIndex, position):
        return self.speciesTree.structuresList[familyIndex][position]

    def getNbSeqs(self, familyIndex):
        return len(self.sequencesLists[familyIndex])

    # Returns the list of paired positions for familyIndex
    def getStruct(self, familyIndex):
        return self.speciesTree.getStruct(familyIndex)

    # Returns the only sequence in the set for the specified familyIndex. There should be only one.
    def getOnlySequence(self, familyIndex):
        desiredList = list(self.sequencesLists[familyIndex])
        return desiredList[0]

    # Returns the set for the specified familyIndex
    def getSequencesListsSet(self, familyIndex):
        return self.sequencesLists[familyIndex]

    # For a specified familyIndex, returns the number of optimal sequences (for a leaf, there should be only the
    # input sequence, so just 1).
    def getNbOfSeqs(self, familyIndex):
        return len(self.sequencesLists[familyIndex])

    # To get a cost inside posToNucleoCosts for a specified familyIndex, position and nucleotide (char)
    def getCostOf(self, familyIndex, position, nucleotide):
        return self.posToNucleoCosts[familyIndex][position][nucleotide]

    def getCostOfBag(self, familyIndex, bagIndex, nucleotides):
        return self.bagsCostMatrix[familyIndex][bagIndex][nucleotides]

    # To return the index of the nucleo, using nucleoToIndexDict (shorter to write)
    def nucToInd(self, nucleo):
        index_of_nucleo = self.nucleoToIndexDict[nucleo]
        return index_of_nucleo

    # To return a translation of the nucleotide, using nucToNucDict
    def translateNuc(self, nucleo):
        return self.nucToNucDict[nucleo]

    # To return the value of the parameter T depending on depth of the node
    def getT(self):
        ## LINEAR GRADIENT ##

        # using a simple linear function now; depth 0 returns 50% (0.5)
        slope = (1 - 0.5) / self.speciesTree.maxDepth  # range is from 50% (for the root) to 100% (for the leaves)
        # --> 100% will never be returned, because this method will be called on the parent node
        t = 0.5 + self.depth * slope

        ## Checking and returning
        if t < 0.499 or t > 1.001:
            print("BIG problem: T is outside of [0.5, 1]")
            sys.exit()

        return t

    # To get the listMutatedSequences for a familyIndex
    def getListMutSeqs(self, familyIndex):
        return self.listMutatedSequences[familyIndex]

    # Prints a description of the node
    def printNode(self, printCosts):
        parentIdString = "None"
        if not self.isRoot():
            parentIdString = str(self.parent.id)

        lcIdString = "None"
        rcIdString = "None"
        if not self.isLeaf():
            lcIdString = str(self.leftChild.id)
            rcIdString = str(self.rightChild.id)

        print("\n--- Id = " + str(
            self.id) + "  Parent = " + parentIdString + "  LeftChild = " + lcIdString + "  RightChild = " + rcIdString + "  Depth = " + str(
            self.depth) + "  T = " + str(self.getT()))

        # NEW
        nbOptSeqs = 0

        for familyIndex in range(self.getNbStructs()):
            nbOptSeqs += self.getNbOfSeqs(familyIndex)
            print(". Family " + str(familyIndex) + ":")
            # test
            printed_seqs = 0
            for seq in self.sequencesLists[familyIndex]:
                if printed_seqs < 100000:
                    print(seq)
                    # The costs might not be printed in order, we might have to sort the keys
                    if printCosts:
                        for pos in range(self.getSeqLength()):
                            costString = "Pos " + str(pos) + ":  "
                            for key in sorted(self.posToNucleoCosts[familyIndex][pos]):
                                costString += key + ": " + str(self.posToNucleoCosts[familyIndex][pos][key]) + "   "
                            print(costString)
                    printed_seqs += 1
                else:
                    print("Too many sequences to print, only the first 100000 are printed.")
                    break

        return nbOptSeqs

    # To return a clone of this node, without the sequences unless it's a leaf node. Parent and child nodes are not
    # set here.
    def clone(self, clonedSpeciesTree):
        cloneNode = SpeciesNode(None, None, None, self.structure_weight, self.bags_ordered_list, clonedSpeciesTree,
                                self.id)  # keeping the same ids
        if self.isLeaf():
            for familyIndex in range(self.getNbStructs()):
                for seq in self.getSequencesListsSet(familyIndex):
                    cloneNode.addSequence(familyIndex, seq)  # copying the sequences, only for the leaves

        return cloneNode

    ###
    ### Mutators
    ###

    # To add a sequence to sequencesList at a certain familyIndex
    # Now checks if the sequence is already present; the method adds the sequence only if it is different
    def addSequence(self, familyIndex, sequence):
        if familyIndex >= self.speciesTree.getNbStructs():
            print("familyIndex : " + str(familyIndex) + " doesn't exist! Sequence cannot be added.")
            return

        if sequence not in self.sequencesLists[familyIndex]:
            self.sequencesLists[familyIndex].add(sequence)  # NEW: this is a set now, not a list

    # To update posToNucleoCosts for a specified familyIndex, a specified position and a specified nucleotide (char)
    # If addToValue is True, then we add theCost to the one that is there, otherwise we set the value
    def updatePosToNucleoCosts(self, familyIndex, position, nucleotide, theCost, addToValue=False):
        if addToValue:
            self.posToNucleoCosts[familyIndex][position][nucleotide] += theCost
        else:
            self.posToNucleoCosts[familyIndex][position][nucleotide] = theCost

    def updateBagCostMatrix(self, familyIndex, bagIndex, nucleotides, theCost, addToValue=False):
        if addToValue:
            if self.bagsCostMatrix[familyIndex][bagIndex][nucleotides] == float('inf'):
                self.bagsCostMatrix[familyIndex][bagIndex][nucleotides] = theCost
            else:
                self.bagsCostMatrix[familyIndex][bagIndex][nucleotides] += theCost
        else:
            self.bagsCostMatrix[familyIndex][bagIndex][nucleotides] = theCost

    # Recursive method to calculate the costs using the algo defined by the parameter 'algo' (Bottom-up step)
    def calculateCosts(self, algo, c_bar):
        if self.isLeaf():  # for the leaves, the costs are either 0 or infinity (infinity is the default cost)
            # test
            # print("-------------------------Leaf node-------------------------")

            for familyIndex in range(self.getNbStructs()):  # for each family
                if self.isLeaf() and self.getNbOfSeqs(familyIndex) != 1:
                    print("---> Huge problem, leaf node doesn't have exactly 1 sequence for family " + str(familyIndex))

                # sequence = self.getSequence(familyIndex, 0)  #There should be only one sequence, at index 0
                sequence = self.getOnlySequence(familyIndex)  # NEW: updated to deal with sets
                if algo != SpeciesTree.ALGO(5):
                    for pos in range(self.getSeqLength()):  # for each position
                        nucleo = sequence[pos]
                        posPair = self.getPairedPos(familyIndex, pos)

                        # If we have to consider the secondary structure
                        if (algo == SpeciesTree.ALGO(3) or algo == SpeciesTree.ALGO(4)) and posPair != -1:
                            pairedNuc = sequence[posPair]

                            # The pair of nucleos as the key.
                            self.updatePosToNucleoCosts(familyIndex, pos, nucleo + pairedNuc, 0)

                        else:  # If the algo doesn't consider the 2D struct or if this position is not paired.
                            # Updating only the cost of the nucleotide that is in the sequence at position pos.
                            self.updatePosToNucleoCosts(familyIndex, pos, nucleo, 0)
                else:
                    for bagIndex in range(len(self.bags_ordered_list)):
                        positions = self.bags_ordered_list[bagIndex]
                        nucleotides = []
                        for pos in positions:
                            nucleotides.append(str(sequence[pos]))
                        self.updateBagCostMatrix(familyIndex, bagIndex, tuple(nucleotides), 0)
            # print("Leaf done.")

        else:  # if the current node is an internal node
            # postorder
            # 147 nodes/clades for GLM-clan
            self.leftChild.calculateCosts(algo, c_bar)  # Calculate for left child first
            self.rightChild.calculateCosts(algo, c_bar)  # then for the right child
            # test
            # print("+++++++++++++++++++++++++Internal node+++++++++++++++++++++++++")

            # When the costs have been calculated for both children nodes, calculate for the current node
            for familyIndex in range(self.getNbStructs()):  # for each family
                # test
                # print("Family", familyIndex)
                if algo != SpeciesTree.ALGO(5):
                    for pos in range(self.getSeqLength()):  # for each position

                        if algo == SpeciesTree.ALGO(1):
                            self.computeFitch(familyIndex, pos)

                        elif algo == SpeciesTree.ALGO(2):
                            self.computeSankoff(familyIndex, pos)

                        elif algo == SpeciesTree.ALGO(3):
                            self.computeSeqAndStruct(familyIndex, pos)

                        elif algo == SpeciesTree.ALGO(4):
                            self.computeSeqAndStruct_2structsApprox(familyIndex, pos)
                else:
                    # For the tree decomposition method, bagIndex is required instead of position.
                    # To compute costs for each bag, we visit them backwardly. (forwardly for making selections)
                    # # test
                    # start_time = time.time()
                    # for bagIndex in reversed(range(len(self.bags_ordered_list))):
                    for bagIndex in range(len(self.bags_ordered_list)):
                        self.computeTreeDecomp(familyIndex, bagIndex)
                    # # test
                    # print("Time for one internal Node of one family:", time.time() - start_time)

            if c_bar is not None:
                c_bar.update(1)

            # test
            # print("Internal Node done.")

    # The top-down step, called after calculateCosts; this method selects the nucleotide of minimum cost for each
    # position and builds the ancestral sequences.
    # This method is called on a node for which we already have selected the optimal sequences, and we find the
    # optimal sequences for its children depending on that sequence.
    def selectNucleos(self, algo, s_bar):
        # Return if both children are leaves.
        if self.leftChild.isLeaf() and self.rightChild.isLeaf():
            return

        for familyIndex in range(self.getNbStructs()):
            # print("selectNucleos: working on familyIndex # " + str(familyIndex))
            # print("There are " + str(self.getNbSeqs(familyIndex)) + " optimal sequences at the parent node")
            for parentSequence in self.getSequencesListsSet(familyIndex):  # for every optimal sequence
                # if seqIndex % 100 == 0:
                # print("### selectNucleos: working on sequence # " + str(seqIndex))

                # using list representation of strings here, so we can modify easily characters
                sequenceLeftList = [[" " for x in range(self.getSeqLength())]]
                sequenceRightList = [[" " for x in range(self.getSeqLength())]]

                # parentSequence = self.getSequence(familyIndex, seqIndex)  #now going through all the optimal sequences  --> now parentSequence is in the loop
                if algo != SpeciesTree.ALGO(5):
                    for pos in range(self.getSeqLength()):
                        # print("selectNucleos: working on position # " + str(pos))
                        nucleo = parentSequence[pos]

                        if algo == SpeciesTree.ALGO(1) or algo == SpeciesTree.ALGO(2):
                            leftAndRightChildNucleos = None

                            if algo == SpeciesTree.ALGO(1):
                                leftAndRightChildNucleos = self.computeFitch(familyIndex, pos, nucleo)

                            elif algo == SpeciesTree.ALGO(2):
                                leftAndRightChildNucleos = self.computeSankoff(familyIndex, pos, nucleo)

                            # print(leftAndRightChildNucleos)
                            self.addNucleoToSequencesList(sequenceLeftList, pos, leftAndRightChildNucleos[0])
                            self.addNucleoToSequencesList(sequenceRightList, pos, leftAndRightChildNucleos[1])

                        elif algo == SpeciesTree.ALGO(3) or algo == SpeciesTree.ALGO(4):
                            posPair = self.getPairedPos(familyIndex, pos)

                            if algo == SpeciesTree.ALGO(3):
                                method = self.computeSeqAndStruct
                            else:
                                method = self.computeSeqAndStruct_2structsApprox

                            if posPair == -1:  # unpaired
                                leftAndRightChildNucleos = method(familyIndex, pos, nucleo)
                                self.addNucleoToSequencesList(sequenceLeftList, pos, leftAndRightChildNucleos[0])
                                self.addNucleoToSequencesList(sequenceRightList, pos, leftAndRightChildNucleos[1])
                            elif pos < posPair:  # paired in the 2D structure, and we did not consider the pos already
                                leftAndRightChildNucleos = method(familyIndex, pos, nucleo + parentSequence[posPair])

                                # A list of di-nucleotides is returned: first position in the list corresponds to the
                                # left child and second position corresponds to right child.
                                self.addNucleoToSequencesList(sequenceLeftList, pos,
                                                              [nuc[0] for nuc in leftAndRightChildNucleos[0]])
                                self.addNucleoToSequencesList(sequenceRightList, pos,
                                                              [nuc[0] for nuc in leftAndRightChildNucleos[1]])
                                self.addNucleoToSequencesList(sequenceLeftList, posPair,
                                                              [nuc[1] for nuc in leftAndRightChildNucleos[0]])
                                self.addNucleoToSequencesList(sequenceRightList, posPair,
                                                              [nuc[1] for nuc in leftAndRightChildNucleos[1]])
                            else:
                                continue
                else:
                    # TODO: need to figure out a way of enumerating sequences without knowing the bag information.
                    # Convert the parent sequence into bag instances.
                    fixed_bags = []
                    for bagID in range(len(self.bags_ordered_list)):
                        cNucleos = []
                        positions = self.bags_ordered_list[bagID]
                        for pos in positions:
                            cNucleos.append(parentSequence[pos])
                        fixed_bags.append(tuple(cNucleos))

                    # Create an empty template of a sequence for the recursive methods to fill and distribute.
                    # Start from the first bag, return at the last.
                    sequence_template_0 = ['_'] * self.sequenceLength
                    optim_seqs_dict_left = [math.inf, []]
                    self.selectNucleos4Bag_rec(familyIndex, 0, sequence_template_0, 0,
                                               optim_seqs_dict_left, fixed_bags=fixed_bags, child_option=0)
                    sequenceLeftList = optim_seqs_dict_left[1]

                    sequence_template_1 = ['_'] * self.sequenceLength
                    optim_seqs_dict_right = [math.inf, []]
                    self.selectNucleos4Bag_rec(familyIndex, 0, sequence_template_1, 0,
                                               optim_seqs_dict_right, fixed_bags=fixed_bags, child_option=1)
                    sequenceRightList = optim_seqs_dict_right[1]

                if len(sequenceLeftList) > 25000000 or len(sequenceRightList) > 25000000:
                    print("Too many sequences to enumerate.")
                    sys.exit()
                # If one of the children nodes is a leaf, we don't add the sequence for this node
                if not self.leftChild.isLeaf():
                    for seq in sequenceLeftList:
                        # Use join to transform the sequence back into a string.
                        self.leftChild.addSequence(familyIndex, "".join(seq))
                if not self.rightChild.isLeaf():
                    for seq in sequenceRightList:
                        self.rightChild.addSequence(familyIndex, "".join(seq))

        if s_bar is not None:
            s_bar.update(1)

        # Calling on the children nodes
        if not self.leftChild.isLeaf():
            self.leftChild.selectNucleos(algo, s_bar)
        if not self.rightChild.isLeaf():
            self.rightChild.selectNucleos(algo, s_bar)

    # The Fitch algorithm to get the costs for one position. If chosenNucleo is specified, then the method will return
    # a list of all left and right children nucleotides that gave rise to this nucleotide in the current node (for
    # top-down step)
    def computeFitch(self, familyIndex, position, chosenNucleo=None, costMatrix=None):

        if costMatrix is None:
            costMatrix = self.FITCH_MATRIX

        for nucleo in self.NUCLEOS:  # for every nucleotide, we want to find the minimum cost
            # only for top-down step; we want to do calculation only for specified nucleo
            if chosenNucleo is not None and nucleo != chosenNucleo:
                continue

            minLeft = float("inf")  # Initializing
            minRight = float("inf")
            minNucLeft = []  # for top-down step only -> now a list of all optimal choices
            minNucRight = []  # for top-down step only -> now a list of all optimal choices

            for nucleoChild in self.NUCLEOS:  # for every possible nucleotide in the left and right children

                costLeft = (self.leftChild.getCostOf(familyIndex, position, nucleoChild) +
                            costMatrix[self.nucToInd(nucleoChild)][self.nucToInd(nucleo)])

                if costLeft == minLeft:
                    minNucLeft.append(nucleoChild)  # appending
                elif costLeft < minLeft:
                    minLeft = costLeft
                    minNucLeft = [nucleoChild]  # starting a new list

                costRight = (self.rightChild.getCostOf(familyIndex, position, nucleoChild) +
                             costMatrix[self.nucToInd(nucleoChild)][self.nucToInd(nucleo)])

                if costRight == minRight:
                    minNucRight.append(nucleoChild)  # appending
                elif costRight < minRight:
                    minRight = costRight
                    minNucRight = [nucleoChild]  # starting a new list

            if chosenNucleo is not None:
                # Returning a list containing the list of optimal nucleos for the left child followed by the list of
                # nucleos for the right child.
                return [minNucLeft, minNucRight]

            # updating the minimum cost for nucleotide nucleo
            self.updatePosToNucleoCosts(familyIndex, position, nucleo, (minLeft + minRight))

    # The Sankoff algorithm, similar to Fitch, but with a different substitution matrix. If chosenNucleo is specified,
    # then the method will return the left and right children nucleotides that gave rise to this nucleotide in the
    # current node (for top-down step)
    def computeSankoff(self, familyIndex, position, chosenNucleo=None):

        # The Sankoff cost matrix:
        costMatrix = self.SANKOFF_MATRIX

        # Just calling the Fitch method with the Sankoff cost matrix
        return self.computeFitch(familyIndex, position, chosenNucleo, costMatrix)

    # The Sankoff algorithm for nucleotide substitutions + the 2D structure change cost. (considering only 1 structure
    # for each sequence here) ChosenNucleo here can be a list of 2 nucleos (1st corresponds to position,
    # 2nd corresponds to posPair).
    def computeSeqAndStruct(self, familyIndex, position, chosenNucleo=None):

        costMatrix = self.SANKOFF_MATRIX  # Using the Sankoff matrix for nucleotide substitution costs
        bpMatrix = self.BASEPAIR_MATRIX_SIMPLE

        # Checking if position is paired with another position posPair in the 2D structure
        posPair = self.getPairedPos(familyIndex, position)

        if posPair == -1:  # unpaired in the 2D structure
            # in this case chosenNucleo should be only one nucleo (not a list)
            return self.computeSankoff(familyIndex, position, chosenNucleo)

        # if the cost has been calculated already because of the base pair
        elif chosenNucleo is None and self.getCostOf(familyIndex, position, 'AA') < float("inf"):
            # print("Position already considered because of the base pair and chosenNucleo = " + str(chosenNucleo))
            return

        else:  # position not considered yet

            for nucleosBP in self.DINUCLEOS:  # for every possible di-nucleotide

                if chosenNucleo is not None and nucleosBP != chosenNucleo:
                    continue

                ### Calculating the minimum total cost to get nucleosBP ###
                minLeft = float("inf")  # Initializing
                minRight = float("inf")
                minNucleosLeft = []  # for top-down step only -> now a list of all optimal choices
                minNucleosRight = []  # for top-down step only -> now a list of all optimal choices

                for nucleosChild in self.DINUCLEOS:  # for every possible nucleotide base pair in the left and right children

                    costLeft = (self.leftChild.getCostOf(familyIndex, position, nucleosChild) +
                                costMatrix[self.nucToInd(nucleosChild[0])][self.nucToInd(nucleosBP[0])] +
                                costMatrix[self.nucToInd(nucleosChild[1])][self.nucToInd(nucleosBP[1])] +
                                bpMatrix[self.nucToInd(nucleosBP[0])][self.nucToInd(nucleosBP[1])])  # adding up the cost of both substitutions and bp cost

                    if costLeft == minLeft:
                        minNucleosLeft.append(nucleosChild)  # appending
                    elif costLeft < minLeft:
                        minLeft = costLeft
                        minNucleosLeft = [nucleosChild]  # starting a new list

                    costRight = (self.rightChild.getCostOf(familyIndex, position, nucleosChild) +
                                 costMatrix[self.nucToInd(nucleosChild[0])][self.nucToInd(nucleosBP[0])] +
                                 costMatrix[self.nucToInd(nucleosChild[1])][self.nucToInd(nucleosBP[1])] +
                                 bpMatrix[self.nucToInd(nucleosBP[0])][self.nucToInd(nucleosBP[1])])  # adding up the cost of both substitutions and bp cost

                    if costRight == minRight:
                        minNucleosRight.append(nucleosChild)  # appending
                    if costRight < minRight:
                        minRight = costRight
                        minNucleosRight = [nucleosChild]  # starting a new list

                if chosenNucleo is not None:
                    return [minNucleosLeft, minNucleosRight]

                seqAndStructTotalCost = minLeft + minRight

                # Adding the cost in posToNucleoCosts for the two positions, with the di-nucleotide as the key
                # in the dict. (first nucleo is the nucleo of the position and second nucleo is always the nucleo
                # at the base pair position)
                self.updatePosToNucleoCosts(familyIndex, position, nucleosBP, seqAndStructTotalCost)
                self.updatePosToNucleoCosts(familyIndex, posPair, nucleosBP[1] + nucleosBP[0],
                                            seqAndStructTotalCost)  # reverse order of the di-nucleotides

    def computeTreeDecomp(self, familyIndex, bagIndex, superAnc=False, chosenNucleos=None):
        # For every possible nucleos in a bag in its childs
        # Cost from substitution and MFE
        # Use SANKOFF_MATRIX for substitution cost at each position.
        # Structure cost however, need to be calculated for every permutation.

        # Walk through every possible permutation of length => bag size.
        positions = self.bags_ordered_list[bagIndex]
        args_bag = [self.NUCLEOS for x in range(len(positions))]

        # Go through all possible combinations of nucleos in the bag.
        # keep the ones with the minimum cost.
        for combination in itertools.product(*args_bag):
            # print(combination)
            # Calculate the cost of each permutation.
            # Accumulate the min-cost options onto the min-costs from left and right child.
            # Replace the positions given by the bag by the combination of nucleos we would like to try.
            minLeft = float("inf")  # Initializing
            minRight = float("inf")
            minNucleosLeft = []  # for top-down step only -> now a list of all optimal choices
            minNucleosRight = []  # for top-down step only -> now a list of all optimal choices

            for combination_child in itertools.product(*args_bag):

                costLeft = self.leftChild.findBagCost(familyIndex, bagIndex, combination_child, combination,
                                                      superAnc)

                if costLeft == minLeft:
                    minNucleosLeft.append(combination_child)  # appending
                elif costLeft < minLeft:
                    minLeft = costLeft
                    minNucleosLeft = [combination_child]  # starting a new list

                costRight = self.rightChild.findBagCost(familyIndex, bagIndex, combination_child, combination,
                                                        superAnc)

                if costRight == minRight:
                    minNucleosRight.append(combination_child)  # appending
                if costRight < minRight:
                    minRight = costRight
                    minNucleosRight = [combination_child]  # starting a new list

            if chosenNucleos:
                return [minNucleosLeft, minNucleosRight]

            totalCost = minLeft + minRight
            # print(totalCost)
            self.updateBagCostMatrix(familyIndex, bagIndex, combination, totalCost)

    # Recursion let's go!
    def selectNucleos4Bag_rec(self, fam, bag_id, sequence_template_orig, prev_cost, optim_sequences_dict, fixed_bags=None, child_option=None, order=0):
        sequence_template = sequence_template_orig.copy()
        cost_nucleos = None
        currentNode = None

        if fixed_bags is None:
            cost_nucleos = self.bagsCostMatrix[fam][bag_id]
        else:
            currentNode = self.leftChild if child_option == 0 else self.rightChild
        bag = self.bags_ordered_list[bag_id]
        args_bag = [['A', 'C', 'G', 'U'] for x in range(len(bag))]
        min_cost = float('inf')
        min_cost_list = []  # A list containing all combinations for a bag that gives the least cost.

        # Find positions in the sequence where values are assigned in previous bag(s).
        for ind in range(len(bag)):
            nucleotide = sequence_template[bag[ind]]
            if nucleotide != '_':
                args_bag[ind] = nucleotide

        for combination in itertools.product(*args_bag):
            if fixed_bags is None:
                cost_tmp = cost_nucleos[combination]
            else:
                cost_tmp = currentNode.findBagCost(fam, bag_id, combination, fixed_bags[bag_id])

            if cost_tmp < min_cost:
                min_cost_list = [combination]
                min_cost = cost_tmp
            elif cost_tmp == min_cost:
                min_cost_list.append(combination)

        new_cost = prev_cost + min_cost

        # TODO: 'min_cost_list' # of branches?
        templates_to_pass = []
        for comb in min_cost_list:
            sequence_template_tmp = sequence_template.copy()
            for ind in range(len(bag)):
                if comb[ind] != sequence_template_tmp[bag[ind]]:
                    sequence_template_tmp[bag[ind]] = comb[ind]
            templates_to_pass.append(sequence_template_tmp)

        # bag_id_next = (bag_id + 1) if order == 0 else (bag_id - 1)
        continue_recurse = True
        if (order == 0) and (bag_id < (len(self.bags_ordered_list) - 1)):
            bag_id_next = bag_id + 1
        elif (order == 1) and (bag_id > 0):
            bag_id_next = bag_id - 1
        else:
            continue_recurse = False

        if continue_recurse:
            for sequence_template_to_pass in templates_to_pass:
                rt_dict = [math.inf, []]
                if fixed_bags is None:
                    self.selectNucleos4Bag_rec(fam, bag_id_next, sequence_template_to_pass, new_cost, rt_dict, order=order)
                else:
                    self.selectNucleos4Bag_rec(fam, bag_id_next, sequence_template_to_pass, new_cost, rt_dict,
                                                         fixed_bags=fixed_bags, child_option=child_option, order=order)
                if rt_dict[0] < optim_sequences_dict[0]:
                    optim_sequences_dict[0] = rt_dict[0]
                    optim_sequences_dict[1] = rt_dict[1]
                elif rt_dict[0] == optim_sequences_dict[0]:
                    optim_sequences_dict[1] += rt_dict[1]
        else:
            # # test
            # print("New cost:", new_cost)
            optim_sequences_dict[0] = new_cost
            optim_sequences_dict[1] = templates_to_pass

    def findBagCost(self, family_index, bag_index, comb_child, comb, super_anc=False):
        subMatrix = self.SANKOFF_MATRIX
        T = self.parent.getT()  # getting the time/gradient parameter
        otherFamilyIndex = 1 if family_index == 0 else 0
        structure_cost = 0.0
        substitutionCost = 0.0

        # Walk through every possible permutation of length => bag size.
        positions = self.bags_ordered_list[bag_index]
        sequenceTemplate = self.sequenceTemplate.copy()
        # structureTemplate = self.structureTemplate
        # otherstructureTemplate = self.structureTemplate
        structsList = self.speciesTree.getStructs_dot_bracket()

        # Replace the positions given by the bag by the actual structure from the structures list.
        # print(comb, comb_child)
        comb_ind = 0
        # TODO: Might need to check the input structures, incomplete pairs are not accepted.
        for pos in positions:
            # structureTemplate[pos] = structsList[family_index][pos]
            # otherstructureTemplate[pos] = structsList[otherFamilyIndex][pos]
            sequenceTemplate[pos] = comb[comb_ind]
            comb_ind += 1

        with suppress_stdout_stderr():
            sequenceStr = "".join(str(x) for x in sequenceTemplate)
            mfe_cost = float(RNA.energy_of_structure(sequenceStr, structsList[family_index], 0))

            if not super_anc:
                mfe_cost_other = float(RNA.energy_of_structure(sequenceStr, structsList[otherFamilyIndex], 0))
                structure_cost = (T * mfe_cost + (1 - T) * mfe_cost_other)
            else:
                mfe_cost_other = float(RNA.energy_of_structure(sequenceStr, structsList[otherFamilyIndex], 0))
                structure_cost = (0.5 * mfe_cost + 0.5 * mfe_cost_other)

            for pos in range(len(comb)):
                substitutionCost += subMatrix[self.nucToInd(comb_child[pos])][self.nucToInd(comb[pos])]

        # Cost composition:
        # - Base cost of having the given combination at the current (child) node.
        # - Substitution cost (SANKOFF).
        # - MFE cost that considers structure information.
        totalCost = (self.getCostOfBag(family_index, bag_index, comb_child) + substitutionCost
                     + structure_cost * self.structure_weight)

        # print("##", self.getCostOfBag(family_index, bag_index, comb_child))
        # print("sub-cost:", substitutionCost, "struct-cost:", structure_cost * self.structure_weight)

        return totalCost

    # The Sankoff algorithm for nucleotide substitutions + the 2D structure change cost (considering the 2
    # structures, but not the full cycles) chosenNucleo here can be a list of 2 nucleos (1st corresponds to position,
    # 2nd corresponds to posPair)
    def computeSeqAndStruct_2structsApprox(self, familyIndex, position, chosenNucleo=None):

        # costMatrix = SpeciesNode.SANKOFF_MATRIX  # Using the Sankoff matrix for nucleotide substitution costs
        # bpMatrix = SpeciesNode.BASEPAIR_MATRIX_SIMPLE
        #
        # minNucsList = []  # A list that will contain in order: left and right minimum di-nucleotides

        posPair = self.getPairedPos(familyIndex,
                                    position)  # Checking if position is paired with another position posPair in the 2D structure

        if posPair == -1:  # unpaired in the 2D structure
            nucleosOrDinucleos = self.NUCLEOS  # will be used in the for loop
        else:
            nucleosOrDinucleos = self.DINUCLEOS
            if chosenNucleo is None and self.getCostOf(familyIndex, position, 'AA') < float(
                    "inf"):  # if the cost has been calculated already because of the base pair
                # print("Position already considered because of the base pair and chosenNucleo = " + str(chosenNucleo))
                return

        for nucleo_s in nucleosOrDinucleos:  # for every possible nucleo or dinucleo (dependending of the paired state) (called nucleo_s because it's either 1 or 2 nucleo(s))

            if chosenNucleo is not None and nucleo_s != chosenNucleo:
                continue

            ### Calculating the minimum total cost to get nucleo_s ###
            minLeft = float("inf")  # Initializing
            minRight = float("inf")
            minNucleosLeft = []  # for top-down step only -> now a list of all optimal choices
            minNucleosRight = []  # for top-down step only -> now a list of all optimal choices

            for nucleo_sChild in nucleosOrDinucleos:  # for every possible nucleotide or nucleotide base pair in the left and right children

                costLeft = self.leftChild.find2StructsApproxCost(familyIndex, position, posPair, nucleo_sChild,
                                                                 nucleo_s)

                if costLeft == minLeft:
                    minNucleosLeft.append(nucleo_sChild)  # appending
                elif costLeft < minLeft:
                    minLeft = costLeft
                    minNucleosLeft = [nucleo_sChild]  # starting a new list

                costRight = self.rightChild.find2StructsApproxCost(familyIndex, position, posPair, nucleo_sChild,
                                                                   nucleo_s)

                if costRight == minRight:
                    minNucleosRight.append(nucleo_sChild)  # appending
                if costRight < minRight:
                    minRight = costRight
                    minNucleosRight = [nucleo_sChild]  # starting a new list

            if chosenNucleo is not None:
                return [minNucleosLeft, minNucleosRight]

            seqAndStructTotalCost = minLeft + minRight

            # Adding the cost in posToNucleoCosts for the two positions, with the di-nucleotide as the key in the dict (first nucleo is the nucleo of the position
            # and second nucleo is always the nucleo at the base pair position)
            # If this position is not paired, we do only the first update
            self.updatePosToNucleoCosts(familyIndex, position, nucleo_s, seqAndStructTotalCost)
            if posPair != -1:
                self.updatePosToNucleoCosts(familyIndex, posPair, nucleo_s[1] + nucleo_s[0], seqAndStructTotalCost)  # reverse order of the di-nucleotides

    # Calculating the cost of going from nucleo_sChild to nucleo_s; **this method is called on the child node**.
    # Here we consider the 2 structures in the cost, but not cycles.
    def find2StructsApproxCost(self, familyIndex, position, posPair, nucleo_sChild, nucleo_s) -> float:

        subMatrix = self.SANKOFF_MATRIX
        bpMatrix = self.BASEPAIR_MATRIX_SIMPLE  # for now, using a very simple matrix
        T = self.parent.getT()  # getting the time parameter

        # intializing with the cost of nucleo_sChild of this node (which is a child node)
        costCurrentStruct = self.getCostOf(familyIndex, position, nucleo_sChild)
        costOtherStruct = 0.0

        # First dealing with the current structure
        if posPair == -1:  # unpaired in current struct
            costCurrentStruct += subMatrix[self.nucToInd(nucleo_sChild)][self.nucToInd(nucleo_s)]

        else:  # paired
            costCurrentStruct += (subMatrix[self.nucToInd(nucleo_sChild[0])][self.nucToInd(nucleo_s[0])] +
                                  subMatrix[self.nucToInd(nucleo_sChild[1])][self.nucToInd(nucleo_s[1])] + T *
                                  bpMatrix[self.nucToInd(nucleo_s[0])][self.nucToInd(nucleo_s[1])])  # T multiplies the bpCost ONLY

        # Second, we deal with the base pairs of the second structure directly connected to position and posPair
        if familyIndex == 0:
            otherFamily = 1
        else:
            otherFamily = 0

        # Whether position is paired or unpaired, we always have to check for a paired position for 'position' in the other structure
        otherFamilyPosPair = self.getPairedPos(otherFamily, position)

        if posPair != -1 and otherFamilyPosPair == posPair:  ## NEW June 9 2016: simpler than what I did before, we can return here if we have same bps in two structs
            return costCurrentStruct + (1 - T) * bpMatrix[self.nucToInd(nucleo_s[0])][self.nucToInd(nucleo_s[1])]  ##this corresponds to 100% of the bpcost (because same basepairs)

        if otherFamilyPosPair != -1:  # otherwise, if == -1, there is nothing to calculate (costOtherStruct stays the same)

            ##NEW### That didn't make a significant difference
            # if posPair != -1 and otherFamilyPosPair == posPair:
            #    print("Special case where we have the same basepairs in both structures")
            #    costOtherStruct += bpMatrix[self.nucToInd(nucleo_s[0])][self.nucToInd(nucleo_s[1])]
            # else:
            ########
            costSum = 0.0
            for nuc in self.NUCLEOS:  # we're not setting the nucleo at position: otherFamilyPosPair, so we take into account all possible nucleotides
                # not counting substitutions for otherFamilyPosPair for now:
                # costSum += bpMatrix[self.nucToInd(nucleo_sChild[0])][self.nucToInd(nuc)]  #we can always access nucleo_sChild[0], whether 'position' is paired or not
                # October 22nd 2015: Pretty sure that the line above had a bug. We need to look at all possible basepairs with the proposed nucleotide, not the child one
                costSum += bpMatrix[self.nucToInd(nucleo_s[0])][self.nucToInd(nuc)]  # we can always access nucleo_sChild[0], whether 'position' is paired or not

            costOtherStruct += costSum / 4  # adding the average cost for the 4 different nucleos possible

        if posPair != -1:  # paired in the current struct, so we have to consider 'posPair' in the other structure

            # checking for posPair
            otherFamilyPosPair = self.getPairedPos(otherFamily, posPair)
            if otherFamilyPosPair != -1:  # otherwise, if == -1, there is nothing to calculate (costOtherStruct stays the same)

                ###NEW###
                # if otherFamilyPosPair == position:
                #    print("Special case where we have the same basepairs in both structures -> we catch it a second time")
                # else:
                #########
                costSum = 0.0
                for nuc in self.NUCLEOS:  # we're not setting the nucleo at position: otherFamilyPosPair, so we take into account all possible nucleotides
                    # not counting substitutions for otherFamilyPosPair for now:
                    # costSum += bpMatrix[self.nucToInd(nucleo_sChild[1])][self.nucToInd(nuc)]  #here it is safe to access nucleo_sChild[1]
                    # October 22nd 2015: Pretty sure that the line above had a bug. We need to look at all possible basepairs with the proposed nucleotide, not the child one
                    costSum += bpMatrix[self.nucToInd(nucleo_s[1])][self.nucToInd(nuc)]  # here it is safe to access nucleo_s[1]

                costOtherStruct += costSum / 4  # adding the average cost for the 4 different nucleos possible
        # print(costCurrentStruct, costOtherStruct)
        return costCurrentStruct + (1 - T) * costOtherStruct  # (1-T) multiplies costOtherStruct, only because it contains only bp costs, as of now


# ======================================================= Tools ======================================================

# To read a string representing a newick tree with multiple sequences separated by a hyphen (-) and return a
# SpeciesTree. Each sequence represents a family (the sequences must be in the same order as the order of the
# structures in structureList).
def getSpeciesTreeFromNewick(newickString, structuresList, structure_w, bags_in_order):
    phyloTree = Phylo.read(StringIO(newickString), "newick")
    # Phylo.draw(phyloTree)
    # print(phyloTree)
    # print(len(phyloTree.root))
    print("Number of leaves:", phyloTree.count_terminals())
    print("Number of internal nodes:", len(phyloTree.get_nonterminals()))
    number_internal_nodes = len(phyloTree.get_nonterminals())

    # Reinitializing SpeciesNode.lastId
    SpeciesNode.lastId = 1
    speciesTree = SpeciesTree(structuresList, structure_w, bags_in_order)
    speciesTree.root = SpeciesNode(None, None, None, structure_w, bags_in_order, speciesTree)

    # Calling the recursive method on the root
    getSpeciesTreeFromNewickRec(speciesTree.root, phyloTree.root, structure_w, bags_in_order, speciesTree)

    # A boolean variable indicating if we are dealing with a tree representing simulated data (in which case all the
    # internal nodes will contain sequences).
    simulTreeBool = True
    if phyloTree.root.name is None:
        simulTreeBool = False  # No name for the root means no simulated sequence -> we are dealing with biological data

    return speciesTree, simulTreeBool, number_internal_nodes  # returning the boolean too


# Recursive method to build the SpeciesTree from the newick tree
def getSpeciesTreeFromNewickRec(speciesNodeInternal, phyloNodeInternal, structure_w, bags_in_order, speciesTree):
    # print("Inside getSpeciesTreeFromNewickRec")
    # print("TEST confidence = " + str(phyloNodeInternal.confidence))
    leftSpeciesNode = None
    rightSpeciesNode = None
    phyloLeftChild = None
    phyloRightChild = None

    # First adding the sequences in the SpeciesNodeInternal
    if phyloNodeInternal.name is not None:
        sequences = phyloNodeInternal.name.split("-")
        for familyIndex in range(len(sequences)):
            speciesNodeInternal.addSequence(familyIndex, sequences[familyIndex])

    # Then creating the SpeciesNode objects representing the children nodes
    counter = 0
    for child in phyloNodeInternal:
        counter += 1
        if counter == 1:
            leftSpeciesNode = SpeciesNode(speciesNodeInternal, None, None, structure_w, bags_in_order,
                                          speciesTree)  # Specifying the parent in the constructor
            speciesNodeInternal.leftChild = leftSpeciesNode  # Creating the parent-child link
            phyloLeftChild = child
        elif counter == 2:
            rightSpeciesNode = SpeciesNode(speciesNodeInternal, None, None, structure_w, bags_in_order,
                                           speciesTree)  # Specifying the parent in the constructor
            speciesNodeInternal.rightChild = rightSpeciesNode  # Creating the parent-child link
            phyloRightChild = child
        else:  # More than 2 children??!!
            print("---> BIG PROBLEM: newick tree is not binary.")
            counter = 2  # fuck the third one

    if counter == 2:  # Should be 0 or 2. Counter == 0: we're working on a leaf, and there is nothing more to do
        # Calling the recursive method on the left and right children nodes
        getSpeciesTreeFromNewickRec(leftSpeciesNode, phyloLeftChild, structure_w, bags_in_order, speciesTree)
        getSpeciesTreeFromNewickRec(rightSpeciesNode, phyloRightChild, structure_w, bags_in_order, speciesTree)


# Reading a file containing simulated trees and returns a list of SpeciesTree objects
def readFile(filename, tw_l, algo, structure_weight=0.0, d_tree=None, backboned='0'):
    file = open(filename, 'r')
    structuresRead = False
    structuresParenthesisList = []  # contains the parenthesis representations of the structures
    structuresPairedPosList = []  # The list of pairedPos lists representing the structures
    speciesTrees = []  # the list of SpeciesTree objects that will be returned
    bags_in_order = None
    simulTreeBool = True
    numberOfInternalNodes = None

    for line in file:

        if line[0] == "#":  # reading a tree index, we don't really need to store this information
            structuresRead = True  # we're done reading structures
            if not structuresPairedPosList:
                for structure in structuresParenthesisList:
                    structuresPairedPosList.append(structure)
                    # structuresPairedPosList.append(parenthesisToPairedPos(structure))
                # print(structuresPairedPosList)
            if (algo == SpeciesTree.ALGO(5)) and not bags_in_order:
                if d_tree:
                    pickled_tree_dir = 'tree_bags' if backboned == '0' else 'backboned_tree_bags'
                    PATH_DTREE = os.path.join('..', 'Resources', pickled_tree_dir, 'dieted_treeDecomp_' + d_tree + '_' + tw_l + '.data')
                    with open(PATH_DTREE, 'rb') as f:
                        bags_in_order = pickle.load(f)
                else:
                    tree_decomposition_object = StructureMerge(structuresParenthesisList)
                    # # test
                    # print("td object created.")
                    # tree = tree_decomposition_object.root_undirected_tree()[0]
                    tw, tree = tree_decomposition_object.get_tree_decomp()
                    # tree_decomposition_object.print_and_plot_outputs()
                    bags_list = list(tree.nodes)

                    # # test
                    # print("tree-decomposition calculated.")
                    # print("tree-width:", tw)
                    # exit()
                    # Tree-diet application.
                    # st_time = time.time()
                    OPT, real_edges, color_dictionary = tree_decomposition_object.get_tree_diet(int(tw_l), False)
                    # test
                    # print("tree-diet done in", time.time() - st_time, ".")
                    for tag in color_dictionary:
                        if tag != -1:
                            color_inst = color_dictionary[tag]
                            bag_o = bags_list[tag]
                            bag_n = []
                            for pos in bag_o:
                                if color_inst[pos] == 1:
                                    bag_n.append(pos)
                            bags_list[tag] = tuple(bag_n)

                    bags_in_order = bags_list
                # test
                # print("number of bags in total:", len(bags_in_order))
                # print(bags_in_order)
                # exit()
                # print("tw before:", tree_decomposition_object.get_tree_decomp())
                # print("tw after:", max(len(x) - 1 for x in bags_in_order))
                # print(bags_in_order)

        elif not structuresRead:
            structuresParenthesisList.append(line.strip())

        else:  # reading a newick string
            speciesTree, simulTreeBool, numberOfInternalNodes = getSpeciesTreeFromNewick(line.strip(),
                                                                                         structuresPairedPosList,
                                                                                         structure_weight, bags_in_order)
            speciesTrees.append(speciesTree)

    return speciesTrees, simulTreeBool, numberOfInternalNodes


# ===================================================== Execution ====================================================

def reconstruction(args):
    algorithm_number = int(args[0])
    if 1 <= algorithm_number <= 5:
        algo = SpeciesTree.ALGO(algorithm_number)
    else:
        print("\nError: invalid algoNumber.")
        sys.exit()

    tree_width_limit = args[1]
    structure_weight = float(args[3])
    dieted_tree = args[4] if len(args) > 4 else None
    bb_or_not = args[5] if len(args) > 5 else '0'
    speciesTreesList, simultTreeBool, numberOfInternalNodes = readFile(args[2], tree_width_limit, algo,
                                                                       structure_weight,
                                                                       dieted_tree, bb_or_not)  # reading the simulFile argument
    startTime = time.time()  # to measure time of computation
    # counter = 1
    pbar_total = (numberOfInternalNodes + 1) if options.rootOnly else (2 * numberOfInternalNodes)
    simulTree = speciesTreesList[0]
    # for simulTree in speciesTreesList:
    #     print("\nTree # " + str(counter))
        # simulTree.printTree()  #printing the original tree (the one from the input file) --> not sure if necessary
        # exit()
    progress_bar = tqdm(total=pbar_total, desc="Reconstruction", miniters=1, position=0,
                        leave=False) if not options.experimentOnly else None

    if simultTreeBool:  # If we are dealing with simulated data
        inferenceTree = simulTree.cloneWithoutAncSeqs()

    else:  # bio data
        inferenceTree = simulTree
    inferenceTree.setDepths()  # necessary

    # with suppress_stdout_stderr():
    inferenceTree.reconstructAncestralSeqs(algo, progress_bar)

    # --> useful to look at the solutions for every node of the tree
    # In order to limit the stdout buffer size, only print the tree when it is Tree-decomp for now.
    # test
    # totalNbOptSeqs = inferenceTree.printTree() if algo == SpeciesTree.ALGO(5) else None
    # inferenceTree.printTree() if (not options.rootOnly) or SpeciesTree.ALGO(5) else None
    # inferenceTree.reconstructAncestralSeqs(algo)
    if (not options.experimentOnly) or (algo == SpeciesTree.ALGO(5)):
        totalNbOptSeqs = inferenceTree.printTree()  # --> useful to look at the solutions for every node of the tree

    if not simultTreeBool:
        if not options.experimentOnly:
            print("\nThere are " + str(totalNbOptSeqs) + " optimal sequences in the inferred tree")
    else:
        (totalNumberOfErrors, totalNumberOfNucleos, totalNumberOfOptimalSeqs, averageNumberOfErrors,
         totalNumberOfNucleos_withStruct, totalNumberOfErrors_withStruct, totalNumberOfNucleos_withoutStruct,
         totalNumberOfErrors_withoutStruct, averageNumberOfErrors_withStruct, averageNumberOfErrors_withoutStruct,
         numberOfNodesInTheTree, totalNumberOfErrors_root, totalNumberOfOptimalSeqs_root,
         averageNumberOfErrors_root, averageNumberOfErrors_root_withoutStruct,
         averageNumberOfErrors_root_withStruct) = inferenceTree.compareTrees(simulTree, inferenceTree)

        for i in range(len(totalNumberOfErrors)):
            if i == 0:
                print("\nFAMILY 0 ONLY:\n")

            elif i == 1:
                print("\nFAMILY 1 ONLY:\n")

            elif i == 2:
                print("\nALL FAMILIES (STRUCTS):\n")

            if not options.rootOnly:
                ###TOTAL
                print("Total: " + str(totalNumberOfErrors[i]) + " errors / " + str(totalNumberOfNucleos[i]) +
                      " nucleotides (" + str(totalNumberOfErrors[i] / totalNumberOfNucleos[i] * 100) + " % error)")
                print("Total for unstructured positions: " + str(totalNumberOfErrors_withoutStruct[i]) +
                      " errors / " + str(totalNumberOfNucleos_withoutStruct[i]) + " nucleotides (" +
                      str(totalNumberOfErrors_withoutStruct[i] / totalNumberOfNucleos_withoutStruct[
                          i] * 100) + " % error)")

                print("Total for structured positions: " + str(
                    totalNumberOfErrors_withStruct[i]) + " errors / " + str(totalNumberOfNucleos_withStruct[i]) +
                      " nucleotides (" + str(
                    totalNumberOfErrors_withStruct[i] / totalNumberOfNucleos_withStruct[i] * 100) + " % error)")

                ###AVERAGE
                print("Average (over " + str(totalNumberOfOptimalSeqs[i]) + " optimal sequences): " + str(
                    averageNumberOfErrors[i]) +
                      " aver. errors per optimal sequence of length = " + str(
                    totalNumberOfNucleos[i] / numberOfNodesInTheTree[i]) +
                      " nucleotides (" + str(averageNumberOfErrors[i] / (
                        totalNumberOfNucleos[i] / numberOfNodesInTheTree[i]) * 100) + " % error)")

                print("Average for unstructured positions: " + str(averageNumberOfErrors_withoutStruct[i]) +
                      " aver. errors per optimal sequence, considering " + str(
                    totalNumberOfNucleos_withoutStruct[i] / numberOfNodesInTheTree[i]) +
                      " nucleotides per seq. (" + str(averageNumberOfErrors_withoutStruct[i] / (
                        totalNumberOfNucleos_withoutStruct[i] / numberOfNodesInTheTree[i]) * 100) + " % error)")

                print("Average for structured positions: " + str(averageNumberOfErrors_withStruct[i]) +
                      " aver. errors per optimal sequence, considering " + str(
                    totalNumberOfNucleos_withStruct[i] / numberOfNodesInTheTree[i]) +
                      " nucleotides per seq. (" + str(averageNumberOfErrors_withStruct[i] / (
                        totalNumberOfNucleos_withStruct[i] / numberOfNodesInTheTree[i]) * 100) + " % error)\n")

            ###ROOT
            print("-For the root only:")
            print("Best: " + str(totalNumberOfErrors_root[i]) + " errors")
            print("Average (over " + str(totalNumberOfOptimalSeqs_root[i]) + " optimal sequences): " +
                  str(averageNumberOfErrors_root[i]) + " aver. errors per optimal sequence of length = " +
                  str(totalNumberOfNucleos[i] / numberOfNodesInTheTree[i]) + " nucleotides (" +
                  str(averageNumberOfErrors_root[i] / (
                          totalNumberOfNucleos[i] / numberOfNodesInTheTree[i]) * 100) + " % error)")

            print("Average for unstructured positions: " + str(averageNumberOfErrors_root_withoutStruct[i]) +
                  " aver. errors per optimal sequence, considering " + str(
                totalNumberOfNucleos_withoutStruct[i] / numberOfNodesInTheTree[i]) +
                  " nucleotides per seq. (" + str(averageNumberOfErrors_root_withoutStruct[i] / (
                    totalNumberOfNucleos_withoutStruct[i] / numberOfNodesInTheTree[i]) * 100) + " % error)")

            print("Average for structured positions: " + str(averageNumberOfErrors_root_withStruct[i]) +
                  " aver. errors per optimal sequence, considering " + str(
                totalNumberOfNucleos_withStruct[i] / numberOfNodesInTheTree[i]) +
                  " nucleotides per seq. (" + str(averageNumberOfErrors_root_withStruct[i] / (
                    totalNumberOfNucleos_withStruct[i] / numberOfNodesInTheTree[i]) * 100) + " % error)")

    elapsedTime = time.time() - startTime
    print("\n### Computation took " + str(elapsedTime) + " seconds.")


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-e", action="store_true", dest="experimentOnly", default=False,
                      help="Infer the ancestral sequences for the experiment only.")
    parser.add_option("-r", action="store_true", dest="rootOnly", default=False,
                      help="Infer the ancestral sequences for the root only.")

    (options, arguments) = parser.parse_args()

    if len(arguments) < 2:
        print("\n-USAGE: python fitchAndCompany.py  [-e]  [-r]  algoNumber treewidth  treeFile  structure-weight")
        print("\n -e: optional parameter to suppress the progress bar.")
        print("\n -r: optional parameter to find ancestral sequences at the root only.")
        print("\nalgoNumber:  1=Fitch  2=Sankoff  3=achARNement 1.0 - one structure  4=achARNement 1.0  5=achARNement 2.0")
        print("treewidth: the size limit for the largest bag of the tree-decomposition based on the structures.")
        print("treeFile: the file containing the 2 structures and a list of trees (simulated or real) in Newick format.")
        print("structure-weight: the weight for the structure cost.")
        sys.exit()

    reconstruction(arguments)
