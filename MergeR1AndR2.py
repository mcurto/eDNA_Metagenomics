from __future__ import division
import sys, os

#make a dictionary with the read name and corresponding length to be used in blast parameters re-estimation
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

#Save blast outputs information in dictionary
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

#Extract identity, query cover, evalue, and bit score from both reads and recalculate these values
def recalculateData(data1, data2):

    #Get number of mismatches from read 1
    id1 = float(data1[1])
    qlen1 = int(data1[2])
    alignement1 = int(data1[3])
    nrMM1 = round((id1 / 100) * alignement1)

    #Get number of mismatches from read 2
    id2 = float(data2[1])
    qlen2 = int(data2[2])
    alignement2 = int(data2[3])
    nrMM2 = round((id2 / 100) * alignement2)

    #recalclutae evalue and bit score by avering out these values across both reads
    avEval = round((float(data1[4]) + float(data2[4])) / 2, 3)
    avBit = round((float(data1[5]) + float(data2[5])) / 2, 3)

    taxid = data1[0]

    # the new query length will be the sum of the length of both reads
    newqlen = qlen1 + qlen2

    #new identity will be the sum of the mismatches devided by the total reads length
    newID = round(((nrMM1 + nrMM2) / newqlen) * 100, 3)

    # new alignement length will be the sum of the alignement length of both reads
    newalignement = alignement1 + alignement2

    newData = [str(i) for i in [taxid, newID, newqlen, newalignement, avEval, avBit]]
    return newData

# recalculate blast parameters in case only one read matches. In this case only query coverage is recalculated
def recalculateIDSingleRead(accData, length):
    newAccData = {}
    for acc, data in accData.items():
        data[2] = str(length)
        newAccData[acc] = data
    return newAccData

# compare parameters from both reads and find which reads matches were found only for R1, only for R2 and for both.
def dict_compare(d1, d2):
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    onlyd1 = d1_keys - d2_keys
    onlyd2 = d2_keys - d1_keys
    return shared_keys, onlyd1, onlyd2

# parse each one of the blast outputs, compare the results from read 1 and read 2, recalculate the parameters and save them into a dictionary
def compare_blastouts(blastoutF1, blastoutF2, paired_lengths):
    result = {}
    # parsing the blast outputs
    print("parsing " + blastoutF1)
    blastout1 = parseBlastout(blastoutF1)
    print("parsing " + blastoutF2)
    blastout2 = parseBlastout(blastoutF2)

    # geting which matches are in common between read 1 and read 2 results
    print("comparing both files")
    ReadComparison = dict_compare(blastout1, blastout2)
    ReadCommon = ReadComparison[0]
    print(str(len(ReadCommon)) + " reads in common")
    onlyBlast1 = ReadComparison[1]
    print(str(len(onlyBlast1)) + " reads only in Read1")
    onlyBlast2 = ReadComparison[2]
    print(str(len(onlyBlast2)) + " reads only in Read2")


    # iterate through results and check if they are matching the same accession
    for read in ReadCommon:
        totalReadLength = paired_lengths[read]
        accData1 = blastout1[read]
        accData2 = blastout2[read]
        newAccData = {}

        # check if reads are matching the same accession
        accComparison = dict_compare(accData1, accData2)
        accCommon = accComparison[0]
        onlyAcc1 = accComparison[1]
        onlyAcc2 = accComparison[2]

        # recalculate parameters for reads with the same match
        for acc in accCommon:
            newData1 = accData1[acc]
            newData2 = accData2[acc]

            newAccData[acc] = recalculateData(newData1, newData2)

        # recalculate parameters for matches uniq to read 1
        for acc in onlyAcc1:
            newData1 = accData1[acc]
            newData1[2] = str(totalReadLength)
            newAccData[acc] = newData1

        # recalculate parameters for matches uniq to read 2
        for acc in onlyAcc2:
            newData2 = accData2[acc]
            newData2[2] = str(totalReadLength)
            newAccData[acc] = newData2

        result[read] = newAccData

    # recalculate parameters for reads where only read 1 has a match
    for read in onlyBlast1:
        totalReadLength = paired_lengths[read]
        result[read] = recalculateIDSingleRead(blastout1[read], totalReadLength)

    # recalculate parameters for reads where only read 2 has a match
    for read in onlyBlast2:
        totalReadLength = paired_lengths[read]
        result[read] = recalculateIDSingleRead(blastout2[read], totalReadLength)

    return result

# Getting positional arguments
blastoutR1 = sys.argv[1]
blastoutR2 = sys.argv[2]
fastaDirs = sys.argv[3].split(",")
outDir = sys.argv[4]

# Go through multiple files and merge results from read 1 and read 2
for b in os.listdir(blastoutR1):
    # File containing read 1 blast results
    blastoutF1 = blastoutR1 + b
    # Get prefix that serves as sample name
    sampleName = b.split(".")[0]
    print("processing sample " + sampleName)
    # File containing read 2 blast results
    blastoutF2 = blastoutR2 + b.replace("forward", "reverse")
    # Get read lengths
    print("getting lengths from paired reads")
    paired_lengths = getFastaLengths(fastaDirs, sampleName)
    print("merging " + blastoutF1 + " and " + blastoutF2)
    # combine results from both reads
    newBlastData = compare_blastouts(blastoutF1, blastoutF2, paired_lengths)
    # save output file into a tab seperated format
    outFile = sampleName + ".R1AndR2Combined_ntAllunpaired.NACorrected.blastout"
    print("writing results to " + outFile)
    with open(outDir + b.split(".")[0] + outFile, "w") as out:
        for read, accData in newBlastData.items():
            for acc, data in accData.items():
                out.write(read + "\t" + acc + "\t" + "\t".join(data) + "\n")
