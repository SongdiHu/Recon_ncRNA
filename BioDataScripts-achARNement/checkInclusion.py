import sys

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("USAGE: python checkInclusion ancestralSeqsFile1 ancestralSeqsFile2")
        sys.exit()

    setSeqsFile1 = set()
    setSeqsFile2 = set()


    file1 = open(sys.argv[1], 'r')
    file2 = open(sys.argv[2], 'r')

    lineCounter = 0
    for line in file1:

        if lineCounter < 2:  #First two lines should represent structures
            lineCounter += 1
            continue
        
        setSeqsFile1.add(line.strip())

    lineCounter = 0
    for line in file2:

        if lineCounter < 2:  #First two lines should represent structures
            lineCounter += 1
            continue
        
        setSeqsFile2.add(line.strip())

    totalNbSeqsFile1 = len(setSeqsFile1)
    seqsIncluded = 0
    #We check if the sequences in file 1 are included in the file 2
    for seq in setSeqsFile1:

        if seq in setSeqsFile2:
            seqsIncluded += 1

    percentageInc = seqsIncluded / totalNbSeqsFile1 * 100

    print(str(percentageInc) + "% of the sequences in " + sys.argv[1] + " are included in " + sys.argv[2])
