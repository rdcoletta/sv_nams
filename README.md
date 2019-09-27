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




## Quality control

Steps:

1. Remove snps in SV before dividing into populations
2. Separate file by pop
3. Percent missing/het data for each population
4. Karyotype to see distribution of markers (SNPs, SVs, SNPs + SVs)
4. Allele frequency for each ril pop


### Remove SNPs that are within the boundaries of a SV (OR JUST DELETIONS??)

Test script on smaller datasets first:

```bash
cd ~/projects/sv_nams/data

head -n 5000 NAM_sv_sniffles_v1.vcf > test_SV_founders.vcf
tail -n 5000 NAM_sv_sniffles_v1.vcf >> test_SV_founders.vcf
head -n 5000 B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > test_SNP_founders.vcf
tail -n 5000 B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf >> test_SNP_founders.vcf
```

Questions:

* Should I exclude SNPs that fall into any SV? Or only if it's in a DEL or INS? What about INV and TRA?
  > I will start by removing any SNP that falls into any type of SV, except translocations (will not count them).

I wrote `scripts/remove_SNPs_within_SVs.py` to create a new vcf file only with SNPs that are not within SV boundaries. **Translocations were not considered**. This script uses a lot of memory and, so far, it's very slow to write such big files. Thus, submitted a `qsub` job using the bash script `scripts/remove_SNPs_within_SVs.sh`, which just runs the python script, to generate the filtered file `data/B73v5.NAM-gatk-snps.not-in-SVs.vcf`:

```bash
# go to project folder
cd ~/projects/sv_nams/scripts
# run bash script
qsub remove_SNPs_within_SVs.sh
```


<mark>TO DO:</mark>
* Expand writing on QC section
* Increase efficiency of `scripts/remove_SNPs_within_SVs.py`


### Creating VCF files for each NAM population
