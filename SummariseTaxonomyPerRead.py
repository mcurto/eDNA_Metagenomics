import sys, os

# Get mean from elements from a list
def Average(lst):
    return sum(lst) / len(lst)

#  Compare two lineages saved in list format and get the most common ancestor among them
def compareTaxonomyList(tax1, tax2):
    result = []
    # Iterate through both lists and while the taxon from both lists is the same save it into the list and continue. At the point that the taxon are different stop iterating.
    for i in range(0, min(len(tax1), len(tax2))):
        if tax1[i] == tax2[i]:
            result.append(tax1[i])
        else:
            break
    return result

# combine taxonomy from multiple equaly good matches of the same read
def combineTaxonomy(taxonomy):
    # parse tax
    taxParsed = []
    taxED = taxonomy.split("&")
    for tax in taxED:
        if tax != "N/A":
            taxParsed += [tax.split(';')]
    # combine taxonomy with common parents
    # in case not lineage is recovered the list of lineages will be empty. If that is the case the final lineage will be outputed as "N/A". Otherwise we combine the taxonomy
    if len(taxParsed) > 0:
        # use the first lineage as a reference to start the lineage comparison
        TaxFinal = taxParsed[0]
        # iterate through the lineages and compare it to the reference
        for i in range(1,len(taxParsed)):
            #make the reference as the the result of the comparison
            TaxFinal = compareTaxonomyList(TaxFinal, taxParsed[i])
            #if there is not taxon in common among all lineages the result will be an empty list. If that is the case the resulting lineage will be "unclassified"
            if len(TaxFinal) == 0:
                 TaxFinal = ["unclassified"]
        # save the resulting lineage as a ";" seperated list
        return ';'.join(TaxFinal)
    else:
        return "N/A"

# Go through all matches, get the most common ancestor per read and save it into the output file
def summarizeTaxonomyPerReads(FiltFile, outFile):
    with open(outFile, "w") as out:
        with open(FiltFile) as f:
            for l in f:
                led = l.rstrip("\n").split("\t")
                read = led[0]
                params = led[3:8]
                lineage = led[-1]
                newLineage = combineTaxonomy(lineage)
                out.write("\t".join([read] + params + [newLineage]) + "\n")

#Get positional arguments
inDir = sys.argv[1]
outDir = sys.argv[2]

#Iterate through all files in the input directory and save the results into a new output file
for Filt in os.listdir(inDir):
    out = outDir + Filt.rstrip(".blastout") + ".TaxPerRead.blastout"
    summarizeTaxonomyPerReads(inDir + Filt, out)
