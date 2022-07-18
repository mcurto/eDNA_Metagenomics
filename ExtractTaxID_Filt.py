import sys

#iterate through records from filtered blast outputs and etract TaxID information
def ExtractTaxID(blastout, fileOut):
    with open(fileOut, "w") as out:
        with open(blastout) as f:
            result = set()
            for l in f:
                TaxIDinfo = l.split("\t")[2].split("&")
                for TID in TaxIDinfo:
                    TIDed = TID.split(";")
                    result.add(TIDed[-1])
            out.write("\n".join(list(result)) + "\n")

ExtractTaxID(sys.argv[1], sys.argv[2])
