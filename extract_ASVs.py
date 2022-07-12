import sys

def extractHead(inTable, outFasta):
    with open(outFasta, 'w') as fasta:
        with open(inTable, 'r') as t:
            OTUs = t.readline().rstrip("\r\n").split("\t")
        OTUNr = 1
        for OTU in OTUs:
            fasta.write(">ASV_{}\n{}\n".format(OTUNr, OTU.strip('"')))
            OTUNr += 1
        print("{} OTUs saved".format(OTUNr - 1))

extractHead(sys.argv[1], sys.argv[2])
        
