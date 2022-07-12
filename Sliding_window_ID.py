import sys
import math
from Bio import SeqIO


def get_nr_s_windows(cov_file, s_window, min_cov):
    print("reading coverage file")
    s_windows = {}
    current_chr = ""
    with open(cov_file) as f:
        for l in f:
            led = l.rstrip("\n").split("\t")
            chr = led[0]
            pos = int(led[1])
            cov = int(led[2])
            if chr != current_chr:
                current_chr = chr
                print("Processing {} ...".format(chr))
            if pos % 1000000 == 0:
                print("{}M positions processed".format(pos / 1000000), end = "\r")

            sw = math.floor(pos / s_window)
            if cov >= min_cov:
                chr_sw = s_windows.get(chr, {})
                chr_sw[sw] = chr_sw.get(sw, 0) + 1
                s_windows[chr] = chr_sw
    print("filtering sliding windows based on minimum depht of {}".format(min_cov))
    final = {}
    for chr, chr_sw in s_windows.items():
        for sw, nr_pos in chr_sw.items():
            if nr_pos == s_window:
                new_chr_sw = final.get(chr, []) + [sw]
                final[chr] = new_chr_sw
    return  final

def parse_vcf(vcf, s_window):
    result = {}
    nr_lines = 0
    with open(vcf) as f:
        for l in f:
            if not l.startswith("#"):
                nr_lines += 1
                if nr_lines % 500000 == 0:
                    print("{}K records processed".format(nr_lines / 1000), end = "\r")
                led = l.split("\t")
                chrm = led[0]
                pos = int(led[1])
                sw = math.floor(pos / s_window)
                chrm_result = result.get(chrm, {})
                chrm_result[sw] = chrm_result.get(sw, 0) + 1
                result[chrm] = chrm_result
    return result

def get_ID(parsed_vcf, sw_data, s_window, out_file):
    with open(out_file, "w") as out:
        for chr, sliding_windows in sw_data.items():
            chrm_result = parsed_vcf.get(chr, {})
            for sw in sliding_windows:
                nr_mm = chrm_result.get(sw, 0)
                ID = round(((s_window - nr_mm) / s_window) * 100, 2)
                out.write("{}\t{}\t{}\n".format(chr, sw * s_window, ID))


s_window = int(sys.argv[1])
print("Getting sliding windows")
s_windows_data = get_nr_s_windows(sys.argv[3], s_window, int(sys.argv[4]))
print("Reading vcf")
vcf_parsed = parse_vcf(sys.argv[2], s_window)
print("Saving result" + " "*20)
get_ID(vcf_parsed, s_windows_data, s_window, sys.argv[5])
