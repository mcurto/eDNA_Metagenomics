import sys

#Extract ASVs from table
def extractHead(inTable, outFasta):
    with open(outFasta, 'w') as fasta:
        #Get first line of file and make list of ASVs
        with open(inTable, 'r') as t:
            OTUs = t.readline().rstrip("\r\n").split("\t")
        OTUNr = 1
        #Iterate through ASVs and save them in fasta format
        for OTU in OTUs:
            fasta.write(">ASV_{}\n{}\n".format(OTUNr, OTU.strip('"')))
            OTUNr += 1
        #Report number of ASVs
        print("{} OTUs saved".format(OTUNr - 1))

#Run all
extractHead(sys.argv[1], sys.argv[2])
