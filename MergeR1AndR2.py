from __future__ import division
import sys, os

def getFastaLengths(fastaDirList, sampleName):
    result = {}
    for fastaDir in fastaDirList:
        fasta1 = fastaDir + sampleName + ".unassembled.forward.fasta"
        fasta2 = fastaDir + sampleName + ".unassembled.reverse.fasta"
        with open(fasta1) as f1, open(fasta2) as f2:
            for l1, l2 in zip(f1, f2):
                if l1.startswith(">"):
                    name = l1[1:].split(" ")[0]
                else:
                    result[name] = len(l1.rstrip("\n")) + len(l2.rstrip("\n"))
    return result

def parseBlastout(blastout):
    result = {}
    with open(blastout) as b:
        for l in b:
            led = l.rstrip("\n").split("\t")
            read = led[0]
            acc = led[1]
            data = led[2:]
            accDict = result.get(read, {})
            accDict[acc] = data
            result[read] = accDict
    return result

def recalculateData(data1, data2):
    id1 = float(data1[1])
    qlen1 = int(data1[2])
    alignement1 = int(data1[3])
    nrMM1 = round((id1 / 100) * alignement1)

    id2 = float(data2[1])
    qlen2 = int(data2[2])
    alignement2 = int(data2[3])
    nrMM2 = round((id2 / 100) * alignement2)

    avEval = round((float(data1[4]) + float(data2[4])) / 2, 3)
    avBit = round((float(data1[5]) + float(data2[5])) / 2, 3)

    taxid = data1[0]
    lineage = data1[-1]

    newqlen = qlen1 + qlen2
    newID = round(((nrMM1 + nrMM2) / newqlen) * 100, 3)
    newalignement = alignement1 + alignement2

    newData = [str(i) for i in [taxid, newID, newqlen, newalignement, avEval, avBit, lineage]]
    return newData


def recalculateIDSingleRead(accData, length):
    newAccData = {}
    for acc, data in accData.items():
        data[2] = str(length)
        newAccData[acc] = data
    return newAccData

def dict_compare(d1, d2):
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    onlyd1 = d1_keys - d2_keys
    onlyd2 = d2_keys - d1_keys
    return shared_keys, onlyd1, onlyd2

def compare_blastouts(blastoutF1, blastoutF2, paired_lengths):
    result = {}
    print("processing " + blastoutF1)
    blastout1 = parseBlastout(blastoutF1)
    print("processing " + blastoutF2)
    blastout2 = parseBlastout(blastoutF2)
    print("comparing both files")
    ReadComparison = dict_compare(blastout1, blastout2)
    ReadCommon = ReadComparison[0]
    print(str(len(ReadCommon)) + " reads in common")
    onlyBlast1 = ReadComparison[1]
    print(str(len(onlyBlast1)) + " reads only in Read1")
    onlyBlast2 = ReadComparison[2]
    print(str(len(onlyBlast2)) + " reads only in Read2")

    for read in ReadCommon:
        totalReadLength = paired_lengths[read]
        accData1 = blastout1[read]
        accData2 = blastout2[read]
        newAccData = {}

        accComparison = dict_compare(accData1, accData2)
        accCommon = accComparison[0]
        onlyAcc1 = accComparison[1]
        onlyAcc2 = accComparison[2]

        for acc in accCommon:
            newData1 = accData1[acc]
            newData2 = accData2[acc]

            newAccData[acc] = recalculateData(newData1, newData2)

        for acc in onlyAcc1:
            newData1 = accData1[acc]
            newData1[2] = str(totalReadLength)
            newAccData[acc] = newData1

        for acc in onlyAcc2:
            newData2 = accData2[acc]
            newData2[2] = str(totalReadLength)
            newAccData[acc] = newData2

        result[read] = newAccData

    for read in onlyBlast1:
        totalReadLength = paired_lengths[read]
        result[read] = recalculateIDSingleRead(blastout1[read], totalReadLength)

    for read in onlyBlast2:
        totalReadLength = paired_lengths[read]
        result[read] = recalculateIDSingleRead(blastout2[read], totalReadLength)

    return result

lineageR1 = sys.argv[1]
lineageR2 = sys.argv[2]
fastaDirs = sys.argv[3].split(",")
outDir = sys.argv[4]


for b in os.listdir(lineageR1):
    blastoutF1 = lineageR1 + b
    sampleName = b.split(".")[0]
    blastoutF2 = lineageR2 + b.replace("forward", "reverse")

    #fasta1 = fastaDir + sampleName + ".unassembled.forward.fasta"
    #fasta2 = fastaDir + sampleName + ".unassembled.reverse.fasta"
    print("geting lengths for parired reads")
    paired_lengths = getFastaLengths(fastaDirs, sampleName)
    print("merging " + blastoutF1 + " and " + blastoutF2)
    newBlastData = compare_blastouts(blastoutF1, blastoutF2, paired_lengths)
    outFile = sampleName + ".R1AndR2Combined_ntAllunpaired.NACorrected.Lineage.blastout"
    print("writing results to " + outFile)
    with open(outDir + b.split(".")[0] + outFile, "w") as out:
        for read, accData in newBlastData.items():
            for acc, data in accData.items():
                out.write(read + "\t" + acc + "\t" + "\t".join(data) + "\n")
