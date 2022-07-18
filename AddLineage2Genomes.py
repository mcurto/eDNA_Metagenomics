import sys, os

# Parse file containing the correspondance between accession number and taxon ID
def ParseAccSpecTaxID(AccSpecTaxID):
    result = {}
    with open(AccSpecTaxID) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            result[led[0]] = led[-1]
    return result

# Parse file containing the correspondance between taxon ID and lineage
def ParseTaxIDLineage(TaxIDLineage):
    result = {}
    with open(TaxIDLineage) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            result[led[0]] = led[1]
    return result

# Iterate throuhg blast output
def AddTaxID_LineageInfo(blastout, outFile, LineageInfo, TaxIDInfo):
    with open(outFile, "w") as out:
        with open(blastout) as f:
            for l in f:
                led = l.rstrip("\n").split("\t")
                # Get accession
                acc = led[1]
                # Get taxon ID rom accession
                TaxID = TaxIDInfo[acc]
                # Get lineage from Taxon ID
                lineage = LineageInfo[TaxID]
                # Add Taxon ID information
                led[2] = TaxID
                # Add lineage to last collumon
                led.append(lineage)
                out.write("\t".join(led) + "\n")

# import accession Taxon ID correspondance
print("processing Acc TaxID info")
TaxIDInfo = ParseAccSpecTaxID(sys.argv[1])

# import Taxon ID lineage correspondance
print("processing Lineage info")
LineageInfo = ParseTaxIDLineage(sys.argv[2])

#Get information from positional arguments
inDir = sys.argv[3]
outDir = sys.argv[4]

# Iterate through blast outputs, process them and save the resutls into new files
for blastout in os.listdir(inDir):
    print("processing file {}".format(blastout))
    outFile = outDir + blastout.replace(".blastout", ".lineage.blastout")
    AddTaxID_LineageInfo(inDir + blastout, outFile, LineageInfo, TaxIDInfo)
