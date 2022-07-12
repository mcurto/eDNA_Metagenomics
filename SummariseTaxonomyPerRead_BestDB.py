import sys, os

# Python program to get average of a list
def Average(lst):
    return sum(lst) / len(lst)

def compareTaxonomyList(tax1, tax2):
    result = []
    for i in range(0, min(len(tax1), len(tax2))):
        if tax1[i] == tax2[i]:
            result.append(tax1[i])
        else:
            break
    return result

# combine taxonomy
def combineTaxonomy(taxonomy):
    # parse tax
    taxParsed = []
    taxED = taxonomy.split("&")
    for tax in taxED:
        if tax != "N/A":
            taxParsed += [tax.split(';')]
    # combine taxonomy with common parents
    if len(taxParsed) > 0:
        TaxFinal = taxParsed[0]
        for i in range(1,len(taxParsed)):
            TaxFinal = compareTaxonomyList(TaxFinal, taxParsed[i])
            if len(TaxFinal) == 0:
                 TaxFinal = ["unclassified"]
        return ';'.join(TaxFinal)
    else:
        return "N/A"

def summarizeTaxonomyPerReads(FiltFile, outFile):
    with open(outFile, "w") as out:
        with open(FiltFile) as f:
            for l in f:
                led = l.rstrip("\n").split("\t")
                data = led[:6]
                lineage = led[-1]
                newLineage = combineTaxonomy(lineage)
                out.write("\t".join(data + [newLineage]) + "\n")

inDir = sys.argv[1]
outDir = sys.argv[2]

for Filt in os.listdir(inDir):
    out = outDir + Filt.rstrip(".blastout") + ".TaxPerRead.blastout"
    summarizeTaxonomyPerReads(inDir + Filt, out)
