import sys

#Get two lists per read: 1 - all acession numbers from that match; 2 - all lineages from that match
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

# Check if there are multiple matches per read to multiple taxa of the same taxonomic group. Also check if it is finding matches to different acessions.
def parse_read_final_data(TaxSum_per_read, blastout_parsed, taxon, out_prefix):
    # Save two outputs:
    #1- summary of number of specific lineages matching per read and number of acessions
    with open(out_prefix + ".check_data.blastout", "w") as out_check:
        #2- filterd results where reads not following the criteria are excluded
        with open(out_prefix + ".filt.blastout", "w") as out_filt:
            # Read merged results from multiple databases
            with open(TaxSum_per_read) as f:
                for l in f:
                    led = l.rstrip("\n").split("\t")
                    read = led[0]
                    lineage = led[-1]
                    #Get lists of acessions and lineages
                    blastout_data = blastout_parsed[read]
                    #List of acessions
                    acc_data = blastout_data[0]
                    #List of lineages
                    lineages = blastout_data[1]
                    nr_lin = 0
                    #Go through lineages
                    for lin in lineages:
                        if lin != "":
                            #make list out of lineage
                            lin_ed = lin.split(";")
                            #if lineage is not the same to the main lineage and it is specific (from the defined taxonomic group) count the record
                            if not lin_ed[-1] in lineage.split(";") and taxon in lin_ed:
                                nr_lin += 1
                    #Get number of acessions
                    nr_acc = len(acc_data)
                    #Save results from number of acessions and lineages to sum file
                    out_check.write("\t".join(led + [str(nr_acc) ,str(nr_lin)]) + "\n")
                    #If multiple macthes to specific group output the match to final file
                    if nr_lin > 0:
                        out_filt.write("\t".join(led + [str(nr_acc) ,str(nr_lin)]) + "\n")

#Import positional argument
blastou_data = parse_blastout(sys.argv[1])
TaxSum_per_read = sys.argv[2]
taxon = sys.argv[3]
out_file = sys.argv[4]

#Run all
parse_read_final_data(TaxSum_per_read, blastou_data, taxon, out_file)
