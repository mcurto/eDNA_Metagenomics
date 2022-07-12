import sys, os

def ParseAccSpecTaxID(AccSpecTaxID):
    result = {}
    with open(AccSpecTaxID) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            result[led[0]] = led[-1]
    return result

def ParseTaxIDLineage(TaxIDLineage):
    result = {}
    with open(TaxIDLineage) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            result[led[0]] = led[1]
    return result


def AddTaxID_LineageInfo(blastout, outFile, LineageInfo, TaxIDInfo):
    with open(outFile, "w") as out:
        with open(blastout) as f:
            for l in f:
                led = l.rstrip("\n").split("\t")
                acc_list = [acc for acc in led[1].split("&")]
                acc = led[1]
                TaxID = TaxIDInfo[acc]
                lineage = LineageInfo[TaxID]
                led[2] = TaxID
                led.append(lineage)
                out.write("\t".join(led) + "\n")

print("processing Acc TaxID info")
TaxIDInfo = ParseAccSpecTaxID(sys.argv[1])

print("processing Lineage info")
LineageInfo = ParseTaxIDLineage(sys.argv[2])

inDir = sys.argv[3]
outDir = sys.argv[4]

for blastout in os.listdir(inDir):
    print("processing file {}".format(blastout))
    outFile = outDir + blastout.replace(".blastout", ".lineage.blastout")
    AddTaxID_LineageInfo(inDir + blastout, outFile, LineageInfo, TaxIDInfo)
