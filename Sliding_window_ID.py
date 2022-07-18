import sys
import math
from Bio import SeqIO

# Used coverage file to define the silding windows to use. Only sliding windowns where all positions have a coverade above the minimum defined are used
def get_nr_s_windows(cov_file, s_window, min_cov):
    print("reading coverage file")
    #Make an empty dictionary to save the sliding windows per chromossome and the number of positions within a sliding window with a coverage above the threashold
    s_windows = {}
    current_chr = ""
    #Open coverage file
    with open(cov_file) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            #Get chromossome, postion and coverage information
            chr = led[0]
            pos = int(led[1])
            cov = int(led[2])
            #Define current cromossome
            if chr != current_chr:
                current_chr = chr
                print("Processing {} ...".format(chr))
            #Print information of the progress since it can take some time
            if pos % 1000000 == 0:
                print("{}M positions processed".format(pos / 1000000), end = "\r")

            #The curent sliding window will be the whole part of the decimal number resulting from the division between the position and the size of the sliding window
            sw = math.floor(pos / s_window)
            # If the position has a coverage above the minimum count it
            if cov >= min_cov:
                chr_sw = s_windows.get(chr, {})
                chr_sw[sw] = chr_sw.get(sw, 0) + 1
                s_windows[chr] = chr_sw
    #Filtering sliding windows based on the number of positions above the threashold
    print("filtering sliding windows based on minimum depht of {}".format(min_cov))
    #Final dict with sliding windows per chromossome
    final = {}
    #Iterate through the previsou dict and if the number of positions counted is equal to the size of the sliding window save it in final
    for chr, chr_sw in s_windows.items():
        for sw, nr_pos in chr_sw.items():
            if nr_pos == s_window:
                new_chr_sw = final.get(chr, []) + [sw]
                final[chr] = new_chr_sw
    return  final

#Parse vcf. Since the vcf file has been filered just to contain positions with SNPs the number of positions found per sliding window will be the number of mismatches per sliding window. This function retrieves this information
def parse_vcf(vcf, s_window):
    result = {}
    nr_lines = 0
    #Read vcf file
    with open(vcf) as f:
        for l in f:
            #Skip header
            if not l.startswith("#"):
                nr_lines += 1
                #Print information of the progress since it can take some time
                if nr_lines % 500000 == 0:
                    print("{}K records processed".format(nr_lines / 1000), end = "\r")
                #Get chromossome, postion and coverage information
                led = l.split("\t")
                chrm = led[0]
                pos = int(led[1])
                #Get current sliding window
                sw = math.floor(pos / s_window)
                #Get the number of positions (mismatches) per sliding window
                chrm_result = result.get(chrm, {})
                chrm_result[sw] = chrm_result.get(sw, 0) + 1
                result[chrm] = chrm_result
    return result

#Calculate the percentage of identity per sliding window
def get_ID(parsed_vcf, sw_data, s_window, out_file):
    with open(out_file, "w") as out:
        #Iterate through chromossomes
        for chr, sliding_windows in sw_data.items():
            chrm_result = parsed_vcf.get(chr, {})
            #Iterate throuhg the filtered sliding windows
            for sw in sliding_windows:
                #Get the number of mismatches calculate identity and save into a file the chromossome, sliding window and identity information
                nr_mm = chrm_result.get(sw, 0)
                ID = round(((s_window - nr_mm) / s_window) * 100, 2)
                out.write("{}\t{}\t{}\n".format(chr, sw * s_window, ID))

# Get sliding window size from positional argument 1
s_window = int(sys.argv[1])
print("Getting sliding windows")
#Get and filter all possible sliding windows
s_windows_data = get_nr_s_windows(sys.argv[3], s_window, int(sys.argv[4]))
print("Reading vcf")
#Parese vcf and get number of mismatches per sliding window
vcf_parsed = parse_vcf(sys.argv[2], s_window)
#Calculate and save identity information
print("Saving result" + " "*20)
get_ID(vcf_parsed, s_windows_data, s_window, sys.argv[5])
