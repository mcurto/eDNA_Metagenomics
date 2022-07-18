import sys, os

#parse information of the correspondance between taxon ID and taxonomic lineage
def parseTaxIDLineage(TaxIDLineage):
    result = {}
    with open(TaxIDLineage) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            if led[1] != "":
                result[led[0]] = led[1]
    return result

# Iterate through filtered blast output file and add lineage information
def AddLineage(LineageData, FiltIn, FiltOut):
    with open(FiltOut, "w") as out:
        with open(FiltIn) as f:
            # list to save taxID for which a lineage is not found
            LineageMissing = []
            lineNr = 0
            for l in f:
                led = l.rstrip("\n").split("\t")
                #Get taxID matching per read
                TaxIDList = [TaxID.split(";")[-1] for TaxID in led[2].split("&")]
                #Make a list of lineages to be added to the record
                LineageList = []
                #Getting taxID, retrive lineage and add it to the list of lineages. In case the lineage is missing add it to the list of missing lineages.
                for TaxID in TaxIDList:
                    Lineage = LineageData.get(TaxID, "N/A")
                    if Lineage == "N/A":
                        LineageMissing.append(TaxID)
                    else:
                        LineageList.append(Lineage)
                # Add the content of the list of lineages to the record
                led.append("&".join(LineageList))
                # output record
                out.write("\t".join(led) + "\n")

            print("Lineage missing for TaxIDs: " + ";".join(list(set(LineageMissing))))

#Process all files in the input directory
def processFiles(TaxIDData, inDir, outDir):
    count = 0
    for f in os.listdir(inDir):
        print("processing " + f)
        out = outDir + f.replace(".blastout", ".NA_corr.blastout")
        AddLineage(TaxIDData, inDir + f, out)

#Import positional arguments
inDir = sys.argv[1]
outDir = sys.argv[3]

#Read taxID lineage correspondance
print("Parse TaxID to lineage file")
TaxIDtoLineage = parseTaxIDLineage(sys.argv[2])

#Adding lineage information
print("Replace data")
processFiles(TaxIDtoLineage, inDir, outDir)
