#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script filter SNPs out of a hapmap file based on coordinates of SNPs
             (chromosome and position) in another tab-delimited file.

Usage: Rscript remove_SNPs_within_SVs [hmp_file] [file_with_SNPs_coord] [output_name]\n\n")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 3 arguments
if (length(args) != 3) {
  stop("incorrect number of arguments provided.
  
Usage: Rscript remove_SNPs_within_SVs [hmp_file] [file_with_SNPs_coord] [output_name]
       ")
}

# assign arguments to variables
hmp.file <- args[1]
snps.file <- args[2]
output.name <- args[3]



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### filtering ----

# load data
cat("Loading data...\n")
hmp.file <- fread(hmp.file, header = TRUE, data.table = FALSE)
snps.file <- fread(snps.file, header = FALSE, data.table = FALSE)
cat("Done!\n\n")

# get chromosomes in the hmp file
chrs.to.parse <- unique(hmp.file[, "chrom"])

# create an empty data frame to store output
output.hmp <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.file)))
colnames(output.hmp) <- colnames(hmp.file)


for (chr in chrs.to.parse) {
  
  cat("Removing SNPs from chr ", chr, "...\n", sep = "")
  
  # make sure snps to remove are in the same chr as the ones in the hapmap
  snps.to.remove <- snps.file[which(snps.file[, 1] == chr), 2]
  hmp.file.filtered <- hmp.file[which(hmp.file[, "chrom"] == chr), ]
  
  # remove SNPs
  hmp.file.filtered <- hmp.file.filtered[which(!hmp.file.filtered[, "pos"] %in% snps.to.remove), ]
  
  # append to final hapmap
  output.hmp <- rbind(output.hmp, hmp.file.filtered)
  
  cat("Done!\n\n")
}


# count how many SNPs were removed
for (chr in chrs.to.parse) {
  
  if (chr %in% unique(snps.file[, 1])) {
    n.snps.before <- NROW(hmp.file[which(hmp.file[, "chrom"] == chr), ])
    n.snps.after <- NROW(output.hmp[which(output.hmp[, "chrom"] == chr), ])
    n.snps.removed <- n.snps.before - n.snps.after
    
    cat(n.snps.removed, "SNPs removed in chromosome", chr, "\n")
  }
}

# write output
fwrite(output.hmp, file = output.name, sep = "\t", quote = FALSE, na = "NA")
