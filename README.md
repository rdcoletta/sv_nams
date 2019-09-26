# Projecting SVs to NAM lines

by Rafael Della Coletta and Candice Hirsch (September, 2019)

> The goal of this analysis is to project/impute the structural variants (SVs) called for the NAM founders into the lines of each NAM population. To do this, we need both SNP and SV calls for the founders, and SNP data for all NAM lines.




## Project folder

All data, scripts, and output of analyses are located on the folder `/home/hirschc1/della028/projects/sv_nams/` from my account at the Minnesota Supercomputing Institute (MSI).




## Transfering data from CyVerse to local folder

On September 26th, Dr. Arun Seetharam shared the data needed for SV projection via CyVerse. The CyVerse path to the data is `/iplant/home/shared/NAM/PANDA/SVs-impute`. The following commands were used to transfer this data to my folder at MSI so I can do my analyses.


```bash
# go to data folder of the project
cd ~/projects/sv_nams/data

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/PANDA/SVs-impute
# check if files match what Arun described
ils
# download data
iget -K NAM_sv_sniffles_v1.vcf.gz
iget -K GBS-output.tar.gz
iget -K B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz

# decompress files
gunzip NAM_sv_sniffles_v1.vcf.gz
tar xvzf GBS-output.tar.gz
gunzip B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz
```

After downloading and decompressing the files, these are the data that I will be using:

* `NAM_sv_sniffles_v1.vcf`: file with SV calls for NAM founders.
* `B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf`: file with SNP calls for NAM founders.
* `GBS-output/populations.snps.vcf`: file with SNP calls (GBS) for all NAM lines.
