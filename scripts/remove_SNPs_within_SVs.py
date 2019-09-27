#!/usr/bin/python3

'''
version 1.0

created by Rafael Della Coletta
2019-09-26
'''


import argparse as ap


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script removes SNPs from a VCF file that fall within the
             boundaries of a structural variant in another VCF file.

note: VCF header from original and filtered files are the same.''')
# add positional arguments
parser.add_argument("vcf_SNPs", type=str,
                    help="VCF file with SNPs calls")
parser.add_argument("vcf_SVs", type=str,
                    help="VCF file with structural variant calls")
parser.add_argument("output_name", type=str,
                    help="name of output file")
# optional argument for which type of variants to exclude?

# pass arguments into variables
args = parser.parse_args()
vcf_SNPs = args.vcf_SNPs
vcf_SVs = args.vcf_SVs
output_name = args.output_name


# funtion to extract positions of SNPs or SVs
def get_variant_info(vcf_file, is_SV=False):

    # define which information to extract from VCF header
    info_list = ["CHROM", "POS", "INFO"]
    # create list to store indices of each info above
    header_idx = []
    # create dictionary to store information extracted
    variant_info = {}

    # read vcf line by line
    for line in vcf_file:
        line = line.strip()
        # if line starts with ## -- skip
        if line[0:2] == "##":
            continue
        # if line starts with # -- save as header
        elif line[0:1] == "#":
            # save header
            header = line.split("\t")
            # remove # from #CHROM
            header[0] = header[0][1:]
            # get indices of columns with important information from header
            for info in info_list:
                header_idx.append(header.index(info))
            # print(header_idx)
        # otherwise, extract info from line
        else:
            line = line.split("\t")
            # get chrom number
            chr = line[header_idx[0]]
            # get snp position
            pos = line[header_idx[1]]

            # if vcf file is from SNPs...
            if is_SV is False:
                # add coordinates to dictionary
                if chr not in variant_info:
                    variant_info[chr] = [pos]
                else:
                    variant_info[chr].append(pos)

            # if vcf file is from SVs...
            if is_SV is True:
                # get sv type
                sv_type = line[header_idx[2]].split(";SVTYPE=")
                sv_type = sv_type[1].split(";")[0]
                # exclude translocations
                if sv_type != "TRA":
                    # get end position of sv
                    sv_end = line[header_idx[2]].split(";END=")
                    sv_end = sv_end[1].split(";")[0]
                    # add sv info to dictionary
                    sv_info = pos + "," + sv_end + "," + sv_type
                    if chr not in variant_info:
                        variant_info[chr] = [sv_info]
                    else:
                        variant_info[chr].append(sv_info)

    return variant_info


# open input files
infile_SNPs = open(vcf_SNPs, "r")
infile_SVs = open(vcf_SVs, "r")

# parse vcf file to get coordinates of variants
print("Extracting SNPs and SVs coordinates...")

SNPs_info = get_variant_info(infile_SNPs, is_SV=False)
SVs_info = get_variant_info(infile_SVs, is_SV=True)
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
    print("  ", len(SNPs_within_SVs[chr]), "SNPs within a SV in", chr + "!")


# open SNP vcf file again
infile_SNPs = open(vcf_SNPs, "r")

# also open output file to write the filtered vcf
outfile = open(output_name, "w")

print("\nGenerating new VCF file only with SNPs not within a SV range...")

# read vcf line by line
for line in infile_SNPs:
    line = line.strip()
    # if line starts with ## -- print line to new vcf
    if line[0:1] == "#":
        print(line, file=outfile)
    # otherwise, extract info from line
    else:
        # get chrom number and snp position
        chr = line.split("\t")[0]
        pos = line.split("\t")[1]
        # make sure the chromosome being parsed has SNPs within a SV
        if chr in SNPs_within_SVs.keys():
            # only print lines to new vcf if positions are not within a SV
            if pos not in SNPs_within_SVs[chr]:
                print(line, file=outfile)

print("Done!")

# close files
infile_SNPs.close()
outfile.close()
