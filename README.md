# Projecting SVs to NAM lines

by Rafael Della Coletta and Candice Hirsch (September-November, 2019)

> The goal of this analysis is to project structural variants (SVs) indentified for the NAM founders onto the RILs of each NAM population. To do this, we need both SNP and SV calls for the founders, and SNP data for all NAM lines.




## Project folder

All data, scripts, and output of analyses are located on the folder `/home/hirschc1/della028/projects/sv_nams/` from my account at the Minnesota Supercomputing Institute (MSI).

```bash
cd ~/projects/

mkdir -p sv_nams/{analysis,data,scripts}
```




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




## Methodology overview


<mark>ADJUST THIS!</mark>

1. Transform VCF files into Hapmap format
2. Remove snps in SV before dividing into populations
2. Separate file by pop
3. Percent missing/het data for each population
4. Karyotype to see distribution of markers (SNPs, SVs, SNPs + SVs)
4. Allele frequency distribution for each ril pop


## Requirements


<mark>ADJUST THIS! AND ADD VERSIONS OF R PACKAGES</mark>

| Software | Version | Libraries / Packages / Plugins |
| -------- | ------- | ------------------------------ |
| R        | 3.6     | `data.table`, `ggplot2`        |
| Python   |         | `argparse`, `pandas`           |
| TASSEL   | 5       |                                |
| vcftools |         |                                |

> Note: most of the bash `for` loops below can all be parallelized for better perfomance. I'm still learning the best way to do it, so I will just show a sequential way of doing that, and you can parallelize the way you think it's best.



## Resequencing data

### Transform VCF into Hapmap format

I transformed VCF files to hapmap files before any kind of analysis because they are much smaller, easier to parse and quicker to analyze. Additionally, hapmap format would be used for projections anyways. The software [TASSEL 5](https://www.maizegenetics.net/tassel) can do this relatively fast with SNPs. However, a small complication arises when doing that for structural variants, since hapmap files were originally designed to store genotypic information as nucleotides. The way we got around that was writing a custom script (`scripts/vcf2hapmap.py`) to consider each SV as binary data (i.e. either present or not present) and code them as "nucleotides". Thus, if a SV is `A`bsent in a genotype, it was coded as `AA`, but if the SV is `T`here, it was coded as `TT`.

Converting SNPs from original VCF to hapmap would a take lot of time using TASSEl because the file is ~23 Gb. Therefore, I performed a series of UNIX commands to split the VCF into chromosomes and scaffolds abd transform each of them into hapmap format. These `.hmp.txt` files will be used to filter out SNPs within SVs, and only then it will be merged back into a single hapmap file (the intermediate files will be deleted).


```bash
# go to project folder
cd ~/projects/sv_nams

# create folder to store intermediate files
mkdir data/tmp

# filter vcf files
for i in {1..10}; do
  grep "^#"  data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > data/tmp/NAM_founders_SNPs.chr$i.vcf
  grep -w "^chr$i" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf >> data/tmp/NAM_founders_SNPs.chr$i.vcf
done
# get scaffolds as well
grep "^#" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > data/tmp/NAM_founders_SNPs.scaffs.vcf
grep "^scaf_" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf >> data/tmp/NAM_founders_SNPs.scaffs.vcf

# transform each filtered vcf in hapmap
for i in {1..10}; do
  run_pipeline.pl -Xmx10g -importGuess data/tmp/NAM_founders_SNPs.chr$i.vcf -export data/tmp/NAM_founders_SNPs.chr$i.hmp.txt -exportType HapmapDiploid
done

# do the same for scaffolds, but this needs to be sorted first (according to TASSEL)
run_pipeline.pl -SortGenotypeFilePlugin -inputFile data/tmp/NAM_founders_SNPs.scaffs.vcf -outputFile data/tmp/NAM_founders_SNPs.scaffs.sorted.vcf -fileType VCF
# transform to diploid format (e.g. "AA" instead of "A")
run_pipeline.pl -Xmx10g -importGuess data/tmp/NAM_founders_SNPs.scaffs.sorted.vcf -export data/tmp/NAM_founders_SNPs.scaffs.hmp.txt -exportType HapmapDiploid
```

> TASSEL throws this error when sorting vcf files `ERROR net.maizegenetics.dna.map.PositionListBuilder - validateOrdering: Position	Chr:SCAF_100	Pos:79721	InsertionPos:0	Name:SSCAF_100_79721	Variants:A/C	MAF:NaN	Ref:A and Position	Chr:SCAF_99	Pos:109272	InsertionPos:0	Name:SSCAF_99_109272	Variants:T/A	MAF:NaN	Ref:T out of order`. However, I think it's just a warning showing which positions were in the wrong position. I'm able to load the sorted vcf file and transform it into hapmap format without problems. Also, no SNP is lost when sorting the file.

Then, I converted the VCF file of SV calls for all NAM founders into hapmap with the following commands:

```bash
# go to project folder
cd ~/projects/sv_nams

# for explanation on how to use the script...
python scripts/vcf2hapmap.py -h

# convert SVs vcf to hmp
python scripts/vcf2hapmap.py data/NAM_sv_sniffles_v1.vcf data/NAM_founders_SVs.not-sorted.hmp.txt

# sort hmp file
run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile data/NAM_founders_SVs.not-sorted.hmp.txt -outputFile data/NAM_founders_SVs.sorted.hmp.txt -fileType Hapmap
# convert to diploid format
run_pipeline.pl -Xmx10g -importGuess data/NAM_founders_SVs.sorted.hmp.txt -export data/NAM_founders_SVs.hmp.txt -exportType HapmapDiploid

# correct typo in a genotype: it's supposed to be IL14H and not IL14
sed -i 1s/IL14.bam/IL14H.bam/ ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt
```

Additional information about the SV is displayed on its ID, since hapmap format doesn't have fields available for adding such information. For example, the ID `del.chr1.51711.71809` on the first column of the hapmap file means that the SV is a deletion on chr1 that starts at 51,711 and ends at 71,809. The second column will also have the chromosome location for that deletion, but the third column will contain the **midpoint position** for that SV (i.e. 61,760). These two columns will always be the coordinates according to the reference genome. Although there will be somewhat redundant information on IDs of most SVs (like DELs, DUPs, INSs, and INVs), the ID will contain very important info about translocations, as it will show the respective location of the TRA in the **non-reference chromosome**.

Importantly, any SV called as heterozygous in the VCF file (i.e. `0/1`) was considered as **not** having a SV, therefore they were coded as `AA`.


### Remove SNPs that are within the boundaries of a SV

SNPs that are found inside SVs are problematic, because they will have segregation issues when you compare multiple lines that have or not that SV. Thus, I wrote `scripts/extract_SNPs_within_SVs.py` to extract SNPs that are within SV boundaries and `scripts/remove_SNPs_within_SVs.R` to remove them. I made two different scripts because it is faster to filter a huge data table than it is to write one line by line.

```bash
# go to project folder
cd ~/projects/sv_nams

# get list of SNP positions to be removed
for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  python scripts/extract_SNPs_within_SVs.py data/tmp/NAM_founders_SNPs.$i.hmp.txt data/NAM_founders_SVs.hmp.txt data/tmp/SNPs_to_remove --SVs_pos
done

# merge SNPs to remove from scaffolds into one file
for scaf_file in $(ls data/tmp/SNPs_to_remove_SCAF_*.txt); do
  cat $scaf_file >> data/tmp/SNPs_to_remove_scaffs.txt
done

# remove SNPs within SVs for each chromosome in parallel
for i in {1..10}; do
  Rscript scripts/remove_SNPs_within_SVs.R data/tmp/NAM_founders_SNPs.chr$i.hmp.txt data/tmp/SNPs_to_remove_$i.txt data/NAM_founders_SNPs.chr$i.not-in-SVs.hmp.txt &
done
# do the same for scaffolds
Rscript scripts/remove_SNPs_within_SVs.R data/tmp/NAM_founders_SNPs.scaffs.hmp.txt data/tmp/SNPs_to_remove_scaffs.txt data/NAM_founders_SNPs.scaffs.not-in-SVs.hmp.txt

# correct typo in a genotype: it's supposed to be M37W and not MS37W
sed -i 1s/MS37W/M37W/ ~/projects/sv_nams/data/NAM_founders_SNPs.*

wc -l data/tmp/NAM_founders_SNPs.*.hmp.txt
# 27,272,341 SNPs total
wc -l data/NAM_founders_SNPs.*.not-in-SVs.hmp.txt
# 21,949,085 SNPs not in SVs


# # remove intermediate files
# rm data/tmp/*
# rmdir data/tmp/
```


## GBS data

### Creating hapmap files for each NAM population

After removing SNPs within SVs from parental data, I had to make sure the GBS data with genotypic data on RILs had the same SNPs as the parents. However, GBS files are very very large, and parsing it would take way too much time. Thus, I split the vcf file with GBS data by chromosome and transformed each file into the hapmap format:

```bash
# go to project folder
cd ~/projects/sv_nams

# create folder to save temporary files
mkdir data/GBS-output/tmp

# print header to each new chromosome vcf file
# awk will print only lines that start with # and quit after first mismatch
# (this is faster than using grep, which will go through the entire file and then quit)
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  awk '{if(/^#/)print;else exit}' data/GBS-output/populations.snps.vcf > data/GBS-output/tmp/NAM_rils_SNPs.$chr.vcf
done

# split vcf file by chromosome using parallel
for i in {1..10}; do
  grep -w '^chr$i' data/GBS-output/populations.snps.vcf >> data/GBS-output/tmp/NAM_rils_SNPs.chr$i.vcf
done
grep "^scaf_" data/GBS-output/populations.snps.vcf >> data/GBS-output/tmp/NAM_rils_SNPs.scaffs.vcf
```

The bottleneck of performance now is that I have 5000+ lines in the vcf file, and transforming this file into hapmap format with TASSEL will take a lot of time. Thus I used [vcftools](https://vcftools.github.io/index.html) to create a vcf file for each NAM family (i.e. 25 vcf files with ~200 RILs each). But first, I had to create a file telling which RILs belong to which family based on the information on http://maizecoop.cropsci.uiuc.edu/nam-rils.php using `scripts/create_file_with_nam_rils_info.py`.

```bash
# go to project folder
cd ~/projects/sv_nams

# get list of all NAM RIL names (and parents) with GBS data
head -n 1 data/GBS-output/populations.sumstats.tsv | cut -f 2 > data/nam_ril_populations.txt
# rearrange information in a table
python scripts/create_file_with_nam_rils_info.py data/nam_ril_populations.txt

# create a file with SNP positions from NAM parents to keep in NAM RILs
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  sed 1d data/NAM_founders_SNPs.$chr.not-in-SVs.hmp.txt | cut -f 3,4 > data/GBS-output/tmp/SNPs-to-keep_from_NAM_founders.$chr.txt
  sed -i "y/SCAF/scaf/" data/GBS-output/tmp/SNPs-to-keep_from_NAM_founders.$chr.txt
  sed -i "s/^[0-9][0-9]*/chr&/1" data/GBS-output/tmp/SNPs-to-keep_from_NAM_founders.$chr.txt
done

# read "data/nam_ril_populations.txt" file line by line for each chromosome
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  {
    # skip header of "nam_ril_populations.txt" file
    read
    # set delimeter to tab
    IFS="\t"
    # read file line by line
    while read -r line; do
      # get name of the cross being parsed
      cross=$(echo $line | cut -f 1)
      # check if directory exists; if it doesnt, create one to store results
      [[ -d data/GBS-output/tmp/$cross ]] || mkdir -p data/GBS-output/tmp/$cross
      # transform line of the file into multiploe lines so that vcftools recognize 1 genotype to keep per line
      echo $line |  tr "\t" "\n" | tr "," "\n" > data/GBS-output/tmp/$cross/genotypes_to_keep.txt
      # use vcftools to filter a vcf file
      vcftools --vcf data/GBS-output/tmp/NAM_rils_SNPs.$chr.vcf \
               --keep data/GBS-output/tmp/$cross/genotypes_to_keep.txt \
               --positions data/GBS-output/tmp/SNPs-to-keep_from_NAM_founders.$chr.txt \
               --out data/GBS-output/tmp/$cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs \
               --recode \
               --recode-INFO-all
    done
  } < "data/nam_ril_populations.txt"
done

# once the above is done, i have to sort each vcf file and export to hapmap format
cd ~/projects/sv_nams/data/GBS-output/tmp/

# commands for sorting
for cross in $(ls -d B73x*); do
  for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  vcf-sort $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.recode.vcf > $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.sorted.vcf
  done
done

# commands for transforming to hapmap
for cross in $(ls -d B73x*); do
  for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  run_pipeline.pl -Xmx10g -importGuess $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.sorted.vcf -export $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.hmp.txt -exportType HapmapDiploid
  done
done
```


### Collapsing duplicated SNPs

After taking a closer look at the hapmap files generated, I noticed that about many SNP positions were duplicated, i.e. there were multiple calls for the same SNP as seen in this example:

| chr | pos | B73 | Parent 2 | RIL 1 | RIL 2 | RIL 3 |
| --- | --- | --- | -------- | ----- | ----- | ----- |
| 1   | 100 | TT  | NN       | NN    | NN    | TA    |
| 1   | 100 | TT  | NN       | NN    | AA    | AA    |
| 1   | 100 | NN  | AA       | NN    | NN    | NN    |

To correct that, I collapsed these duplicates with `scripts/collapse_GBS_markers.R` by letting only the unique call for each RIL and converting to `NN` in case there was a among the calls for that RIL. The duplicated SNPs showed in the previous table would be collapsed into:

| chr | pos | B73 | Parent 2 | RIL 1 | RIL 2 | RIL 3 |
| --- | --- | --- | -------- | ----- | ----- | ----- |
| 1   | 100 | TT  | AA       | NN    | AA    | NN    |


```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# collapse duplicated SNPs
for cross in $(ls -d B73x*); do
  for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  Rscript ~/projects/sv_nams/scripts/collapse_GBS_markers.R $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.hmp.txt $cross
  done
done
```

After collapsing the duplicated SNPs, I merged the all hapmap files from each chromosome into one file.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# merge hapmap files
for cross in $(ls -d B73x*); do
  cat $cross/NAM_rils_SNPs.$cross.chr1.not-in-SVs.collapsed.hmp.txt > $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt
  for chr in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
    sed 1d $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.collapsed.hmp.txt >> $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt
  done
done
```



### Overlaying resequencing data into parental GBS data

The NAM parents were genotyped by both resequencing and GBS. Thus, it is possible that some SNP calls disagree between the two methods. To remove such SNPs, I overlayed the resequencing data into the GBS data and turned a SNP call into `NN` if there was a disagreement using `scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R`. This script produces a hapmap for each cross containing only the parental data for that cross.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt ~/projects/sv_nams/data/NAM_founders_SNPs.chr1.not-in-SVs.hmp.txt $cross $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt
done
```

> Note: Apparently, the parent Tzi8 doesn't have gbs data. Therefore, I used the entire resequencing data for that parent.













Summary raw gbs data:

write final table
plot final distribution of missing data all crosses
plot distribution -- missing data from rils of all families
plot distribution -- missing snps from all families

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

Rscript ~/projects/sv_nams/scripts/summary_raw_gbs.R ~/projects/sv_nams/data/GBS-output/tmp ~/projects/sv_nams/analysis/qc/raw_gbs
```






Selected only polymorphic, then SNPs present in more than 30% of RILs, and use sliding window approach to remove incorrect calls (resulted in about 30k SNPs) with `select_best_SNPs_per_pop.R`

```bash
mkdir ~/projects/sv_nams/analysis/qc

# # filter snps
# Rscript ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.R B73xB97 \
#                                                               B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.not-imputed.hmp.txt \
#                                                               B73xB97/NAM_gbs-parents_SNPs.B73xB97.not-in-SVs.reseq-overlay.hmp.txt \
#                                                               ~/projects/sv_nams/analysis/qc/filter_best_SNPs \
#                                                               --max_missing=0.3 --window_size=15 --window_step=1 --min_snps_per_window=5

# cd ~/projects/sv_nams/data/GBS-output/tmp/
# for cross in $(ls -d B73x*); do
#   echo "Rscript ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.R $cross $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt ~/projects/sv_nams/analysis/qc/filter_best_SNPs --max_missing=0.3 --window_size=15 --window_step=1 --min_snps_per_window=5"
# done > ~/projects/sv_nams/scripts/commands_for_select-best-markers.txt
#
# qsub select_best_SNPs_per_pop_large.sh


qsub select_best_SNPs_per_pop.sh


wc -l B73x*/*.not-imputed.best-markers.hmp.txt
#

```

> Note: The reason why there is `not-imputed` in the filename is because during preliminary tests I tried to impute GBS SNPs using FSFHap from TASSEL to decrease the number of missing data. But after selecting the best markers, I found the imputing SNPs didn't reduce much the missing data and it was actually causing some troubles later during SV projection. So we decided not to impute SNPs at this stage.





**My guess is that projection will need to be tuned for different crosses...**
**Probably will have to select top 30k markers, instead of percent present in RILs and/or include monomorphic as well**




Impute missing SNPs with FSFHap for each cross:


```bash
# create pedigree file for FSFHap plugin from tassel
cd ~/projects/sv_nams
python scripts/id2pedigree.py data/nam_ril_populations.txt data/nam_rils_pedigree.txt


# impute missing data without fill gaps
cd ~/projects/sv_nams/data/GBS-output/tmp/

# raw gbs
# run_pipeline.pl -Xmx10g -h B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.not-imputed.hmp.txt -FSFHapImputationPlugin -pedigrees ~/projects/sv_nams/data/nam_rils_pedigree.txt -logfile B73xB97/fsfhap_raw-gbs.log -cluster true -bc false -fillgaps false -endPLugin -export B73xB97/NAM_rils_SNPs -exportType HapmapDiploid
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -h $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt -FSFHapImputationPlugin -pedigrees ~/projects/sv_nams/data/nam_rils_pedigree.txt -logfile $cross/fsfhap_raw-gbs.log -cluster true -bc false -fillgaps false -endPLugin -export $cross/NAM_rils_SNPs -exportType HapmapDiploid
done

# best markers
# run_pipeline.pl -Xmx10g -h B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.not-imputed.best-markers.hmp.txt -FSFHapImputationPlugin -pedigrees ~/projects/sv_nams/data/nam_rils_pedigree.txt -logfile B73xB97/fsfhap_best-markers.log -cluster true -bc false -fillgaps false -endPLugin -export B73xB97/NAM_rils_SNPs_best-markers -exportType HapmapDiploid
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -h $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt -FSFHapImputationPlugin -pedigrees ~/projects/sv_nams/data/nam_rils_pedigree.txt -logfile $cross/fsfhap_best-markers.log -cluster true -bc false -fillgaps false -endPLugin -export $cross/NAM_rils_SNPs_best-markers -exportType HapmapDiploid
done


# change file name of output from tassel
cd ~/projects/sv_nams/data/GBS-output/tmp/

# raw gbs
# mv B73xB97/NAM_rils_SNPs2.hmp.txt B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.FSFHap-imputed.hmp.txt
for cross in $(ls -d B73x*); do
  mv $cross/NAM_rils_SNPs2.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.hmp.txt
done

# best makers
# mv B73xB97/NAM_rils_SNPs_best-markers2.hmp.txt B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.FSFHap-imputed.best-markers.hmp.txt
for cross in $(ls -d B73x*); do
  mv $cross/NAM_rils_SNPs_best-markers2.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.best-markers.hmp.txt
done


```






QC raw GBS and filtered set before and after imputation with `plot_rils_karyotypes.R`:


```bash

### CREATE COMMAND FILE TO RUN THIS IN PARALLEL (USE ECHO TO PRINT ALL COMMANDS)
# note: make sure to that you enable X11 forwarding before running this script

module load R/3.6.0

# Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt \
#                                                          ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed \
#                                                          B73xB97 \
#                                                          ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_raw-gbs \
#                                                          B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.not-imputed.hmp.txt \
#                                                          B73xB97/NAM_gbs-parents_SNPs.B73xB97.not-in-SVs.reseq-overlay.hmp.txt \
#                                                          --rils=random --parents_in_data=TRUE --overlay_reseq=TRUE

# karyotypes before FSFHap imputation

# raw-gbs
cd ~/projects/sv_nams/data/GBS-output/tmp/
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_raw-gbs $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=random --parents_in_data=TRUE --overlay_reseq=TRUE
done

for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_before-FSFHap_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_before-FSFHap_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/missing_data_before_imputation_raw-gbs.txt
done


# best markers
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_best-markers $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=random --parents_in_data=TRUE --overlay_reseq=FALSE
done

for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_before-FSFHap_$cross\_best-markers_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_before-FSFHap_$cross\_best-markers_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/missing_data_before_imputation_best-markers.txt
done

# 0.22392 (window size 10, minimum SNPs 5)
# 0.19517 (window size 10, minimum SNPs 3)
# 0.55202 (window size 15, minimum SNPs 8)
# 0.52677 (no sliding window)


# karyotypes after FSFHap imputation

# raw gbs
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_raw-gbs/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/after_imputation_raw-gbs $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=$rils --parents_in_data=FALSE --overlay_reseq=FALSE
done

for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_after-FSFHap_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_after-FSFHap_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/missing_data_after_imputation_raw-gbs.txt
done



# best markers
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_raw-gbs/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/after_imputation_best-markers $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.best-markers.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt --rils=$rils --parents_in_data=FALSE --overlay_reseq=FALSE
done

for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.best-markers.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_after-FSFHap_$cross\_best-markers_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_after-FSFHap_$cross\_best-markers_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/missing_data_after_imputation_best-markers.txt
done


#### B73xHp301 and B73xM162W actually increased number of missing data


# 0.19249 (window size 10, minimum SNPs 5)
# 0.16625 (window size 10, minimum SNPs 3)
# 0.47576 (window size 15, minimum SNPs 8)
# 0.41132 (no sliding window)
```



Make sure SNPs from FSFHap imputed RILs have the same name as the ones from parents

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/correct_SNP-names_rils.R $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.best-markers.hmp.txt
done
```






## Merge SNPs with SVs


```bash
module load R/3.6.0

cd ~/projects/sv_nams/data/GBS-output/tmp/

# merge svs

# Rscript ~/projects/sv_nams/scripts/merge_SVs_and_SNPs.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt B73xB97/NAM_gbs-parents_SNPs.B73xB97.not-in-SVs.reseq-overlay.hmp.txt B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.FSFHap-imputed.best-markers.consensus.correct-marker-names.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.B73xB97.hmp.txt ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.B73xB97.best-markers.consensus.not-projected.hmp.txt

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/merge_SVs_and_SNPs.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.best-markers.correct-marker-names.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt
done


for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt
done


```



Quick summary:


```bash
# cd ~/projects/sv_nams/data/GBS-output/tmp/
#
# cross="B73xB97"
#
# # reseq svs
# cut -f 1,2,3,4,5,6,7,8,9,10,11,12,14 ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt > $cross/$cross\_parents_SVs.hmp.txt
# run_pipeline.pl -Xmx6g -importGuess $cross/$cross\_parents_SVs.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_reseq-parents_SVs_$cross\_OverallSummary
# grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_reseq-parents_SVs_$cross\_OverallSummary1.txt
#
# # gbs parents
# run_pipeline.pl -Xmx6g -importGuess $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_gbs-parents_SNPs_$cross\_OverallSummary
# grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_gbs-parents_SNPs_$cross\_OverallSummary1.txt
#
# # gbs rils before fsfhap
# run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_before-FSFHap_$cross\_OverallSummary
# grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_before-FSFHap_$cross\_OverallSummary1.txt
#
# # gbs rils after fsfhap
# run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.correct-marker-names.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_after-FSFHap_$cross\_OverallSummary
# grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/NAM_rils_SNPs_after-FSFHap_$cross\_OverallSummary1.txt
#
#
# wc -l ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt
# wc -l $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt
# wc -l $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt
# wc -l $cross/NAM_rils_SNPs.$cross.not-in-SVs.FSFHap-imputed.correct-marker-names.hmp.txt
```





## Projection


```bash

cd ~/projects/sv_nams/data/GBS-output/tmp/

# sort parents sv+snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt -outputFile ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -export ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -exportType HapmapDiploid
done

# sort rils sv+snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt -outputFile ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -export ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -exportType HapmapDiploid
done

# just number of rows
for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt
done



# create new folder
mkdir ~/projects/sv_nams/analysis/imputation

# create command files
for cross in $(ls -d B73x*); do
  # create haplotypes from parents
  echo "run_pipeline.pl -Xmx10g -FILLINFindHaplotypesPlugin  -hmp ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt  -o ~/projects/sv_nams/analysis/imputation/donors_$cross -hapSize 1000 -minTaxa 1"
done > ~/projects/sv_nams/scripts/commands_for_create_haplotypes_for_projection.txt

for cross in $(ls -d B73x*); do
  # impute ril genotypes based on
  echo "run_pipeline.pl -Xmx10g -FILLINImputationPlugin -hmp ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -d ~/projects/sv_nams/analysis/imputation/donors_$cross -o ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -hapSize 1000 -accuracy"
done > ~/projects/sv_nams/scripts/commands_for_project_SVs.txt


qsub create_haplotypes_for_projection.sh
qsub project_SVs.sh


for cross in $(ls -d B73x*); do
  # convert to hapmap diploid
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -export ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -exportType HapmapDiploid
done

```



QC SV projection:

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

Rscript ~/projects/sv_nams/scripts/count_projected_SVs.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt ~/projects/sv_nams/analysis/imputation


# best markers
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/after_imputation_best-markers/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/SV_projection ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt --rils=$rils
done



for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/imputation/missing_data_after_SV-projection_best-markers.txt
done


# 0.82248 (window size 10, minimum SNPs 5)
# 0.81671 (window size 10, minimum SNPs 3)
# 0.88475 (window size 15, minimum SNPs 8)
# 0.86701 (no sliding window)

```

Count SVs with `count_projected_SVs.R`:
* window size 10, minimum SNPs 5 = 70k
* window size 10, minimum SNPs 3 = 70k
* window size 15, minimum SNPs 8 = 50k
* no sliding window: 22k = 20k
* no filtering at all = 4k





## Redo projection without imputation

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/
B73x*/*.not-imputed.best-markers.hmp.txt



#### make sure SNPs from FSFHap imputed RILs have the same name as the ones from parents
cd ~/projects/sv_nams/data/GBS-output/tmp/
module load R/3.6.0

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/correct_SNP-names_rils.R $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt
done


#### merge svs
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/merge_SVs_and_SNPs.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.hmp.txt
done


for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.hmp.txt
done



#### projection
cd ~/projects/sv_nams/data/GBS-output/tmp/

# sort parents sv+snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt -outputFile ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -export ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -exportType HapmapDiploid
done

# sort rils sv+snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.hmp.txt -outputFile ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.sorted.hmp.txt -export ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.sorted.hmp.txt -exportType HapmapDiploid
done

# just number of rows
for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.sorted.hmp.txt
done



# create new folder
mkdir ~/projects/sv_nams/analysis/projection_no-imp

# create command files
for cross in $(ls -d B73x*); do
  # create haplotypes from parents
  echo "run_pipeline.pl -Xmx10g -FILLINFindHaplotypesPlugin  -hmp ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt  -o ~/projects/sv_nams/analysis/projection_no-imp/donors_$cross -hapSize 1000 -minTaxa 1"
done > ~/projects/sv_nams/scripts/commands_for_create_haplotypes_for_projection_no-imp.txt

for cross in $(ls -d B73x*); do
  # impute ril genotypes based on
  echo "run_pipeline.pl -Xmx10g -FILLINImputationPlugin -hmp ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.sorted.hmp.txt -d ~/projects/sv_nams/analysis/projection_no-imp/donors_$cross -o ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -hapSize 1000 -accuracy"
done > ~/projects/sv_nams/scripts/commands_for_project_SVs_no-imp.txt


qsub ~/projects/sv_nams/scripts/create_haplotypes_for_projection_no-imp.sh
qsub ~/projects/sv_nams/scripts/project_SVs_no-imp.sh



for cross in $(ls -d B73x*); do
  # convert to hapmap diploid
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -export ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -exportType HapmapDiploid
done



#### QC projection
cd ~/projects/sv_nams/data/GBS-output/tmp/

Rscript ~/projects/sv_nams/scripts/count_projected_SVs.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt ~/projects/sv_nams/analysis/projection_no-imp

for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_best-markers/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/SV_projection_no-imp ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt --rils=$rils
done


for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/projection_no-imp/NAM_rils_SVs-SNPs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/projection_no-imp/missing_data_after_SV-projection_not-imputed_best-markers.txt
done

```



## Re-run projection without combining haplotypes in het regions (`-hybNN false`)


```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# create new folder
mkdir ~/projects/sv_nams/analysis/projection_no-imp_corrected

# create command files
for cross in $(ls -d B73x*); do
  # create haplotypes from parents
  echo "run_pipeline.pl -Xmx10g -FILLINFindHaplotypesPlugin  -hmp ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt  -o ~/projects/sv_nams/analysis/projection_no-imp_corrected/donors_$cross -hapSize 1000 -minTaxa 1"
done > ~/projects/sv_nams/scripts/commands_for_create_haplotypes_for_projection_no-imp.txt

for cross in $(ls -d B73x*); do
  # impute ril genotypes based on
  echo "run_pipeline.pl -Xmx10g -FILLINImputationPlugin -hmp ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.not-projected.sorted.hmp.txt -d ~/projects/sv_nams/analysis/projection_no-imp_corrected/donors_$cross -o ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -hapSize 1000 -accuracy -hybNN false"
done > ~/projects/sv_nams/scripts/commands_for_project_SVs_no-imp.txt


qsub ~/projects/sv_nams/scripts/create_haplotypes_for_projection_no-imp.sh
qsub ~/projects/sv_nams/scripts/project_SVs_no-imp.sh



for cross in $(ls -d B73x*); do
  # convert to hapmap diploid
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -export ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -exportType HapmapDiploid
done



#### QC projection
cd ~/projects/sv_nams/data/GBS-output/tmp/

Rscript ~/projects/sv_nams/scripts/count_projected_SVs.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt ~/projects/sv_nams/analysis/projection_no-imp_corrected

for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_best-markers/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/SV_projection_no-imp_corrected ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt --rils=$rils
done


for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.$cross.not-imputed.best-markers.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/projection_no-imp_corrected/missing_data_after_SV-projection_not-imputed_best-markers.txt
done

```


## Put back SNPs that were filtered out into hapmap of each population

`add_monomorphic_SVs.R`



##### Filter ONLY SVs:

```bash
# make sure it's sorted and let TASSEL correct the alleles' column
run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs.hmp.txt -outputFile ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs.hmp.txt -fileType Hapmap
run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs.hmp.txt -export ~/projects/sv_nams/analysis/projection_no-imp_corrected/NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs.hmp.txt -exportType HapmapDiploid

cp NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs.hmp.txt NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs_v2.hmp.txt

```


Upload data to Cyverse

```bash
# go to data folder of the project
cd ~/projects/sv_nams/data

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/Misc
# check if files match what Arun described
ils
# upload data
iput -K NAM_rils_SVs-SNPs.all-crosses.only-projected.SVs_v2.hmp.txt
# exit iRods
iexit full
```





## Panzea SNPs

Plot karyotypes: p1, p2 and hets for each of the 3 RILs used as example in previous steps.

Downloaded `NAM_map_and_genos-120731.zip` at <http://cbsusrv04.tc.cornell.edu/users/panzea/download.aspx?filegroupid=8>. Unziped at `data/NAM_SNPchip`.

<mark>TO DO:</mark>:
* Use `wget` with the link to download directly to MSI and do the below analysis there.

Merge chromosome files into one hapmap and convert to diploid format using tassel:

```bash

cd ~/OneDrive/University\ of\ Minnesota/PhD/hirsch_lab/projects/sv_nams/data/NAM_SNPchip/hapmap

cat NAM_SNP_genos_raw_20090921_chr1.hmp > NAM_SNP_genos_raw_20090921.hmp.txt
for chr in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10; do
  sed 1d NAM_SNP_genos_raw_20090921_$chr.hmp >> NAM_SNP_genos_raw_20090921.hmp.txt
done

run_pipeline.pl -Xmx1g -importGuess NAM_SNP_genos_raw_20090921.hmp.txt -export NAM_SNP_genos_raw_20090921.hmp.txt -exportType HapmapDiploid

```
