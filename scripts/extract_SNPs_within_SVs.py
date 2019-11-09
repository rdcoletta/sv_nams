#!/usr/bin/python3

'''
version 1.1

created by Rafael Della Coletta
2019-09-26
'''


import argparse as ap
import pandas as pd
from natsort import natsorted


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script removes SNPs from a hapmap file that fall within the
             boundaries of a structural variant in another hapmap file.''')
# add positional arguments
parser.add_argument("hmp_SNPs", type=str,
                    help="hapmap file with SNPs calls")
parser.add_argument("hmp_SVs", type=str,
                    help="hapmap file with structural variant calls")
parser.add_argument("output_name", type=str,
                    help="name of output file")
parser.add_argument("--SVs_pos", action="store_true",
                    help="output BED file with SVs positions")
# optional argument for which type of variants to exclude?

# pass arguments into variables
args = parser.parse_args()
hmp_SNPs = args.hmp_SNPs
hmp_SVs = args.hmp_SVs
output_name = args.output_name


# funtion to extract positions of SNPs or SVs
def get_variant_info(hmp_file, is_SV=False):

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
        # if hapmap file is from SNPs...
        if is_SV is False:
            # get snp position
            pos = line[3]
            # add coordinates to dictionary
            if chr not in variant_info:
                variant_info[chr] = [pos]
            else:
                variant_info[chr].append(pos)
        # if hapmap file is from SVs...
        if is_SV is True:
            # get sv type
            sv_type = line[0].split(".")[0]
            # get sv location
            if sv_type == "tra":
                # get position on reference genome for TRA
                # NOTE: still need to get END position with Arun...
                sv_start = line[3]
                sv_end = str(int(sv_start) + 1)
            else:
                sv_start = line[0].split(".")[2]
                sv_end = line[0].split(".")[3]
            # add sv info to dictionary
            sv_info = sv_start + "," + sv_end + "," + sv_type
            if chr not in variant_info:
                variant_info[chr] = [sv_info]
            else:
                variant_info[chr].append(sv_info)

    return variant_info


# open input files
infile_SNPs = open(hmp_SNPs, "r")
infile_SVs = open(hmp_SVs, "r")

# parse hapmap file to get coordinates of variants
print("Extracting SNPs and SVs coordinates...")

SNPs_info = get_variant_info(infile_SNPs, is_SV=False)
SVs_info = get_variant_info(infile_SVs, is_SV=True)

if args.SVs_pos:
    with open(output_name + ".bed", "w") as bedfile:
        # print header
        print("chrom", "chromStart", "chromEnd", sep="\t", file=bedfile)
        for chr in SVs_info.keys():
            for coord in SVs_info[chr]:
                sv_start = coord.split(",")[0]
                sv_end = coord.split(",")[1]
                print(chr, sv_start, sv_end, sep="\t", file=bedfile)
    # sort bed file by chrom and position
    bed_table = pd.read_table(output_name + ".bed", sep="\t",
                              keep_default_na=False, dtype="unicode")
    # make sure columns are the correct type
    bed_table.chrom = bed_table.chrom.astype('category')
    bed_table.chrom.cat.reorder_categories(natsorted(set(bed_table.chrom)),
                                           inplace=True, ordered=True)
    bed_table.chromStart = bed_table.chromStart.astype('int32')
    # sort by chromosome and then by start position
    bed_sorted = bed_table.sort_values(["chrom", "chromStart"])
    # write sorted hapmap
    bed_sorted.to_csv(output_name + ".bed", sep="\t", index=False)

print("Done!")

# close files
infile_SNPs.close()
infile_SVs.close()


# the following code is a bit cumbersome, but it's a quicker way to find which
# SNPs are in the range of a SV than using a bunch of loops inside loops


# create a dictionary to store all bases that are in a SV (e.g., if SV ranges
# from position 100 to 105, I will store the bases 100,101,102,103,104,105)
SVs_range = {}
# create list to store positions of SNPs to be excluded from analysis
SNPs_within_SVs = {}

print("\nIdentifying SNPs within a SV...")

# for each SV in each chromosome
for chr in SVs_info.keys():
    # make sure chromosome is also in SNPs dictionary
    if chr in SNPs_info.keys():
        for SV in SVs_info[chr]:
            # get SV range
            sv_start = int(SV.split(",")[0])
            sv_end = int(SV.split(",")[1])
            sv_range = list(range(sv_start, sv_end + 1, 1))
            # add all bases of the SV range into the dictionary
            if chr not in SVs_range:
                SVs_range[chr] = [sv_range]
            else:
                # SVs_range[chr] will be a list with many sublists inside
                SVs_range[chr].append(sv_range)
        # transform all sublists into one flat list
        SVs_range[chr] = [base for range in SVs_range[chr] for base in range]
        # remove redundancy
        SVs_range[chr] = set(SVs_range[chr])
        # keep only SNPs that are within the range
        if chr not in SNPs_within_SVs:
            SNPs_within_SVs[chr] = [SNP for SNP in SNPs_info[chr] if int(SNP) in SVs_range[chr]]
        else:
            [SNPs_within_SVs[chr].append(SNP) for SNP in SNPs_info[chr] if int(SNP) in SVs_range[chr]]
        # erase all values for that dictionary key to save memory
        SVs_range[chr] = 0
        print("  ", len(SNPs_within_SVs[chr]), "SNPs within a SV in chr",
              chr + "!")
        # open output file to write table with SNPs inside SVs
        with open(output_name + "_" + chr + ".txt", "w") as outfile:
            for pos in SNPs_within_SVs[chr]:
                print(chr, pos, sep="\t", file=outfile)

print("Done!")
