import SpeciesNode
###
### Class MutatedSequence: contains a sequence and has methods to find all mutated sequences that can be reached from this sequence
###
class MutatedSequence:

    # Constructor
    def __init__(self, sequence=None, leftChild=None, rightChild=None, cost=None):

        self.sequence = sequence  # the sequence representing this object (string)
        self.leftChild = leftChild  # the left child from which this sequence was obtained by a mutation process (bottom-up)
        self.rightChild = rightChild  # the right child from which this sequence was obtained by a mutation process (bottom-up)
        self.cost = cost  # the cost of childSequence + the cost of going from childSequence to this sequence

    ###
    ### Accessors
    ###

    # Returns a set of sequences (strings) that can be reached by doing a maximum of nbMuts mutations from this MutatedSequence.
    # isLeftChild is true if we're looking for sequences on the left branch of the SpeciesTree, false if it's for the right branch
    def getSetReachableMutSeqsStrings(self, nbMuts):

        # print(self.sequence)
        setMutSeqsStrings = set(
            [self.sequence])  # a set of strings, self.sequence represents the sequence obtained after 0 mutations

        for x in range(nbMuts):  # repeating this nbMuts times

            newListMut = []  # temporary list

            for mutSeqString in setMutSeqsStrings:

                # print(mutSeqString)

                for pos in range(len(mutSeqString)):  # going through all the positions

                    nucAtPos = SpeciesNode.nucToNucDict[mutSeqString[pos]]  # translating just in case
                    for nuc in SpeciesNode.NUCLEOS:
                        if nucAtPos != nuc:  # not changing a nucleotide for the same nucleotide
                            newSequence = mutSeqString[:pos] + nuc + mutSeqString[pos + 1:]
                            # print("Before mutation, sequence was = " + mutSeqString + "  and new sequence = " + newSequence)
                            newListMut.append(newSequence)
            setMutSeqsStrings.update(newListMut)  # adding the whole list to the set

        return setMutSeqsStrings  # the set of sequences (strings) that can be reached by a max of nbMut mutations from this MutatedSequence Object

    # Returns the cost of going from the sequence represented by this object, to the mutSeq parameter (string)
    # We don't count the cost of the whole sequence, just the effect of the modified nucleotides, and we add that to the cost of this object
    # currentStruct and otherStruct are lists of paired positions, T is the T parameter coming from the SpeciesNode
    def getCostOfMutSeq(self, mutSeq, currentStruct, otherStruct, T):

        subMatrix = SpeciesNode.SANKOFF_MATRIX
        bpMatrix = SpeciesNode.BASEPAIR_MATRIX_SIMPLE

        if len(mutSeq) != len(self.sequence):
            print("Big problem, sequences are not the same length.")
            print("mutSeq: " + mutSeq + " and self.sequence: " + self.sequence)

        listMutatedPos = []
        cost = self.cost  # initializing with the cost of this MutatedSequence
        for pos in range(len(mutSeq)):

            if mutSeq[pos] != self.sequence[pos]:
                listMutatedPos.append(pos)  # keeping track of the mutated positions

                # print("TEST positional argument: " + self.sequence[pos] + "   " + mutSeq[pos])
                cost += subMatrix[SpeciesNode.nucleoToIndexDict[self.sequence[pos]]][
                    SpeciesNode.nucleoToIndexDict[mutSeq[pos]]]  # substitutions are not multiplied by T

                # basepair cost, currentStruct
                posPair = currentStruct[pos]

                if posPair != -1 and posPair not in listMutatedPos:  # only if it is paired and if we didn't already consider this position in the cost because it was mutated too
                    bpCostAfter = bpMatrix[SpeciesNode.nucleoToIndexDict[mutSeq[pos]]][SpeciesNode.nucleoToIndexDict[
                        mutSeq[posPair]]]  # we just want the bpCost of the positions in mutSeq

                    cost += T * bpCostAfter

                # basepair cost, otherStruct
                posPair = otherStruct[pos]

                if posPair != -1 and posPair not in listMutatedPos:  # only if it is paired in otherStruct and the position wasn't already considered in the cost
                    bpCostAfter = bpMatrix[SpeciesNode.nucleoToIndexDict[mutSeq[pos]]][SpeciesNode.nucleoToIndexDict[
                        mutSeq[posPair]]]  # we just want the bpCost of the positions in mutSeq

                    cost += (1 - T) * bpCostAfter

        return cost  # we return cost, which corresponds to the transition cost

    # Clone
    def clone(self):

        return MutatedSequence(self.sequence, self.leftChild, self.rightChild, self.cost)
