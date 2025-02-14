import os
import sys
import re
from math import exp
from subprocess import check_output
from tempfile import NamedTemporaryFile as NTF
from operator import itemgetter

R = 0.0019872041 
T = 310.15

RNAFOLD = 'RNAfold'

def weight_in_ensemble(seq, ss):


    def frequency_ensemble(weight, en_ensemble):
        return exp((-en_ensemble+weight)/(R*T)) 

    def energy_freq(out_rnafold):

        out_rnafold = out_rnafold.split(b'\n')
        en = out_rnafold[2]
        en = float(en[en.rfind(b'(') + 1:en.rfind(b')')])
        en_ensemble = float(re.match(b'[^\d]+(-\d+\.\d+)[^\d]+', out_rnafold[3]).groups()[0])

        return en, en_ensemble

    tmp_file = NTF(dir='.', delete=False)
    tmp_file.write(bytes('>\n%s' % seq, 'UTF-8'))
    tmp_file.close()

    try:
        mfe = check_output([RNAFOLD, '--noPS', '-p0'], stdin=open(tmp_file.name, "r"))
    except:
        print("WTF")
        print(seq)
        print(ss)
        sys.exit(1)
    finally:
        os.remove(tmp_file.name)

    tmp_file = NTF(dir='.', delete=False)
    tmp_file.write(bytes('>\n%s\n%s' % (seq, ss.replace('.', 'x')), 'UTF-8'))
    tmp_file.close()

    try:
        const = check_output([RNAFOLD, '--noPS', '-p0', '-C'], stdin=open(tmp_file.name, "r"))
    except:
        print("WTF")
        print(seq)
        print(ss)
        sys.exit(1)
    finally:
        os.remove(tmp_file.name)

    _, mfe_ensemble = energy_freq(mfe)
    constr_en, _ =  energy_freq(const)
    #print(mfe_ensemble)
    return frequency_ensemble(mfe_ensemble, constr_en), constr_en


if __name__ == '__main__':
    #seq = 'AUUAAACCGCACAGAAAUUGAACUCUGUGCGCUGCCGGUGGGUACACCGGGGGGUGUAGACUAUUCGUGACGAGGGACGCUGUCGCGCAGGCGAAUCUAC'
    #ss = '.......((((((((........))))))))((.((((((....)))))).))..(((((...(((((.....(.(((...))).)....))))))))))'
    #print(weight_in_ensemble(seq, ss))


    if len(sys.argv) < 2:
        print("USAGE: python freq_in_ensemble.py file_structsAndOptSeqs")
        sys.exit()
        
    theFile = open(sys.argv[1], 'r')

    lineCounter = 0
    structuresList = []
    listOfDict = []
    seqID = 1
    for line in theFile:

        if lineCounter < 2:  #The first two lines should contain the two structures
            structuresList.append(line.strip())
            lineCounter += 1

        else:  #the rest of the lines represent optimal sequences
            theSeq = line.strip()
            freqS1, energyS1 = weight_in_ensemble(theSeq, structuresList[0])
            freqS2, energyS2 = weight_in_ensemble(theSeq, structuresList[1])
            #print(str(energyS1) + "  " + str(energyS2))
            ratio = (2 * freqS1 * freqS2) / (freqS1 + freqS2)
            listOfDict.append({"seq":theSeq, "freqS1":freqS1, "freqS2":freqS2, "ratio":ratio, "ID":seqID, "energyS1":energyS1, "energyS2":energyS2})
            seqID += 1

    print("Optimal sequences ordered by freqS1:")

    newlist = sorted(listOfDict, key=itemgetter('freqS1'), reverse=True) 

    for i in range(len(newlist)):
        print(str(i+1) + ") " + newlist[i]["seq"] + " (" + str(newlist[i]["ID"]) + ")  (freqS1 = " + str(newlist[i]["freqS1"]) + ", enS1 = " + str(newlist[i]["energyS1"]) + ")")

    print("\n")

    print("Optimal sequences ordered by freqS2:")

    newlist = sorted(listOfDict, key=itemgetter('freqS2'), reverse=True) 

    for i in range(len(newlist)):
        print(str(i+1) + ") " + newlist[i]["seq"] + " (" + str(newlist[i]["ID"]) + ")  (freqS2 = " + str(newlist[i]["freqS2"]) + ", enS2 = " + str(newlist[i]["energyS2"]) + ")")

    print("\n")

    print("Optimal sequences ordered by ratio (2 * freqS1 * freqS2) / (freqS1 + freqS2):")

    newlist = sorted(listOfDict, key=itemgetter('ratio'), reverse=True) 

    for i in range(len(newlist)):
        print(str(i+1) + ") " + newlist[i]["seq"] + " (" + str(newlist[i]["ID"]) + ")  (ratio = " + str(newlist[i]["ratio"]) + ")")


    #For table output
        
    print("\n")

    print("Optimal sequences ordered by ratio (freqS1 * freqS2) / (freqS1 + freqS2), table output:\n")

    print("Sequence & Ratio & EnergyS1 & FreqS1 & EnergyS2 & FreqS2")
    
    for i in range(len(newlist)):
        print(newlist[i]["seq"] + " & " + str('{:.2e}'.format(newlist[i]["ratio"])) + " & " + str(newlist[i]["energyS1"]) + " & " + str('{:.2e}'.format(newlist[i]["freqS1"])) + " & " + str(newlist[i]["energyS2"]) + " & " + str('{:.2e}'.format(newlist[i]["freqS2"])))
