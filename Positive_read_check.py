import sys

def parse_blastout(blastout):
    result = {}
    with open(blastout) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            read = led[0]
            acc = led[1]
            lineage = led[-1]
            temp = result.get(read, [[],[]])
            temp_acc = temp[0] + [acc]
            temp_lineage = temp[1] + [lineage]
            temp_acc = list(set(temp_acc))
            temp_lineage = list(set(temp_lineage))
            result[read] = [temp_acc, temp_lineage]
    return result

def parse_read_final_data(TaxSum_per_read, blastout_parsed, taxon, out_prefix):
    with open(out_prefix + ".check_data.blastout", "w") as out_check:
        with open(out_prefix + ".filt.blastout", "w") as out_filt:
            with open(TaxSum_per_read) as f:
                for l in f:
                    led = l.rstrip("\n").split("\t")
                    read = led[0]
                    lineage = led[-1]
                    blastout_data = blastout_parsed[read]
                    acc_data = blastout_data[0]
                    lineages = blastout_data[1]
                    nr_lin = 0
                    for lin in lineages:
                        if lin != "":
                            lin_ed = lin.split(";")
                            if not lin_ed[-1] in lineage.split(";") and taxon in lin_ed:
                                nr_lin += 1
                    nr_acc = len(acc_data)
                    out_check.write("\t".join(led + [str(nr_acc) ,str(nr_lin)]) + "\n")
                    if nr_lin > 0:
                        out_filt.write("\t".join(led + [str(nr_acc) ,str(nr_lin)]) + "\n")


blastou_data = parse_blastout(sys.argv[1])
TaxSum_per_read = sys.argv[2]
taxon = sys.argv[3]
out_file = sys.argv[4]

parse_read_final_data(TaxSum_per_read, blastou_data, taxon, out_file)
