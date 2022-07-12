import sys, os


def parseTaxIDLineage(TaxIDLineage):
    result = {}
    with open(TaxIDLineage) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            if led[1] != "":
                result[led[0]] = led[1]
    return result


def AddLineage(LineageData, FiltIn, FiltOut):
    with open(FiltOut, "w") as out:
        with open(FiltIn) as f:
            LineageMissing = []
            lineNr = 0
            for l in f:
                led = l.rstrip("\n").split("\t")
                TaxIDList = [TaxID.split(";")[-1] for TaxID in led[2].split("&")]
                LineageList = []
                for TaxID in TaxIDList:
                    Lineage = LineageData.get(TaxID, "N/A")
                    if Lineage == "N/A":
                        LineageMissing.append(TaxID)
                    else:
                        LineageList.append(Lineage)
                led.append("&".join(LineageList))
                out.write("\t".join(led) + "\n")

            print("Lineage missing for TaxIDs: " + ";".join(list(set(LineageMissing))))

def processFiles(TaxIDData, inDir, outDir):
    count = 0
    for f in os.listdir(inDir):
        print("processing " + f)
        out = outDir + f.replace(".blastout", ".NA_corr.blastout")
        AddLineage(TaxIDData, inDir + f, out)

inDir = sys.argv[1]
outDir = sys.argv[3]


print("Parse TaxID to lineage file")
TaxIDtoLineage = parseTaxIDLineage(sys.argv[2])

print("Replace data")
processFiles(TaxIDtoLineage, inDir, outDir)
