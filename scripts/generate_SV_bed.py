#!/usr/bin/python3

'''
version 1.0

created by Rafael Della Coletta
2019-11-27
'''


import argparse as ap
import pandas as pd
from natsort import natsorted


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script generates a BED file containing deletions and
             translocations smaller than 1 Mb from a hapmap file containing
             structural variant calls.''')
# add positional arguments
parser.add_argument("hmp_SVs", type=str,
                    help="hapmap file with structural variant calls")
parser.add_argument("output_name", type=str,
                    help="name of output file")
# optional argument for which type of variants to exclude?

# pass arguments into variables
args = parser.parse_args()
hmp_SVs = args.hmp_SVs
output_name = args.output_name


# funtion to extract positions of SNPs or SVs
def get_variant_info(hmp_file):

    # create dictionary to store information extracted
    variant_info = {}
    # skip hapmap header
    hmp_file.readline()
    # read hapmap line by line
    for line in hmp_file:
        line = line.strip()
        line = line.split("\t")
        # get chrom number
        chr = line[2]
        # get sv type
        sv_type = line[0].split(".")[0]
        # keep only deletions or translocations
        if sv_type == "tra" or sv_type == "del":
            # get sv location
            if sv_type == "tra":
                # get position on reference genome for TRA
                # NOTE: still need to get END position with Arun...
                sv_start = line[3]
                sv_end = str(int(sv_start) + 1)
            if sv_type == "del":
                sv_start = line[0].split(".")[2]
                sv_end = line[0].split(".")[3]
            # add sv info to dictionary
            sv_info = sv_start + "," + sv_end + "," + sv_type
            if chr not in variant_info:
                variant_info[chr] = [sv_info]
            else:
                variant_info[chr].append(sv_info)

    return variant_info


# open input file
infile_SVs = open(hmp_SVs, "r")

# parse hapmap file to get coordinates of variants
print("Extracting SV coordinates...")

SVs_info = get_variant_info(infile_SVs)

with open(output_name, "w") as bedfile:
    # print header
    print("chrom", "chromStart", "chromEnd", sep="\t", file=bedfile)
    for chr in SVs_info.keys():
        for coord in SVs_info[chr]:
            sv_start = coord.split(",")[0]
            sv_end = coord.split(",")[1]
            # don't consider SVs bigger than 1Mb
            if abs(int(sv_end) - int(sv_start)) < 1000000:
                print(chr, sv_start, sv_end, sep="\t", file=bedfile)
# sort bed file by chrom and position
bed_table = pd.read_table(output_name, sep="\t",
                          keep_default_na=False, dtype="unicode")
# make sure columns are the correct type
bed_table.chrom = bed_table.chrom.astype('category')
bed_table.chrom.cat.reorder_categories(natsorted(set(bed_table.chrom)),
                                       inplace=True, ordered=True)
bed_table.chromStart = bed_table.chromStart.astype('int32')
# sort by chromosome and then by start position
bed_sorted = bed_table.sort_values(["chrom", "chromStart"])
# write sorted hapmap
bed_sorted.to_csv(output_name, sep="\t", index=False)

print("Done!")

# close file
infile_SVs.close()
