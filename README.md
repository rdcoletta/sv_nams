# Projecting SVs to NAM lines

by Rafael Della Coletta and Candice Hirsch (September-December, 2019)

> The goal of this analysis is to project structural variants (SVs) indentified in the NAM founders onto the RILs of each NAM population. To do this, we need both SNP and SV calls for the founders, and SNP data for all NAM lines.




## Project folder

All data, scripts, and output of analyses are located on the folder `/home/hirschc1/della028/projects/sv_nams/` from my account at the Minnesota Supercomputing Institute (MSI).

```bash
cd ~/projects/

mkdir -p sv_nams/{analysis,data,scripts}
```




## Transfering data from CyVerse to local folder

On November 27th, Dr. Arun Seetharam shared the data needed for SV projection via CyVerse. The CyVerse path to the data is `/iplant/home/shared/NAM/PANDA/SVs-impute`. The following commands were used to transfer this data to my folder at MSI so I can do my analyses.


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
iget -K NAM-structural-variations-v2.0.vcf.gz
iget -K GBS-output.tar.gz
iget -K B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz

# decompress files
gunzip NAM-structural-variations-v2.0.vcf.gz
tar xvzf GBS-output.tar.gz
gunzip B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz
```

After downloading and decompressing the files, these are the data that I will be using:

* `NAM-structural-variations-v2.0.vcf`: file with SV calls for NAM founders.
* `B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf`: file with SNP calls for NAM founders.
* `GBS-output/populations.snps.vcf`: file with SNP calls (GBS) for all NAM lines.

**UPDATE:** on February 5th, Arun sent me the newest version of the SV calls via Slack. This version corrects the high amount of missing data called on previous versions. The filename is `data/NAM-structural-variations-v3.0.txt` and will be used for the rest of the analysis.



## Requirements

| Software | Version | Additional libraries / modules                                                         |
| -------- | ------- | -------------------------------------------------------------------------------------- |
| R        | 3.3.3   | `data.table (v1.12.4)`, `ggplot2 (v3.2.0)`, `foreach (v1.4.7)`, `doParallel (v1.0.15)` |
| Python   | 3.6.6   | `argparse (v1.1)`, `pandas (v0.23.4)`, `natsort (v6.0.0)`                              |
| TASSEL   | 5.2.56  | -                                                                                      |
| vcftools | 0.1.17  | -                                                                                      |
| plink    | 1.9     | -                                                                                       |

> Note: most of the bash `for` loops below can all be parallelized for better perfomance using [GNU parallel](https://www.gnu.org/software/parallel/). I'm show a sequential way of doing that because it's easier to understand.



## Resequencing data

### Transform VCF into Hapmap format

I transformed VCF files to hapmap files before any kind of analysis because they are much smaller, easier to parse and quicker to analyze. Additionally, hapmap format would be used for projections anyways. The software [TASSEL 5](https://www.maizegenetics.net/tassel) can do this relatively fast with SNPs. However, a small complication arises when doing that for structural variants, since hapmap files were originally designed to store genotypic information as nucleotides. The way we got around that was writing a custom script (`scripts/vcf2hapmap.py`) to consider each SV as binary data (i.e. either present or not present) and code them as "nucleotides". Thus, if a SV is `A`bsent in a genotype, it was coded as `AA`, but if the SV is `T`here, it was coded as `TT`.

Converting SNPs from original VCF to hapmap would a take lot of time using TASSEl because the file is ~23 Gb. Therefore, I performed a series of UNIX commands to split the VCF into chromosomes and scaffolds and transform each of them into hapmap format.


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

# correct typo in a genotype: it's supposed to be M37W and not MS37W
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  sed -i 1s/MS37W/M37W/ data/tmp/NAM_founders_SNPs.$chr.hmp.txt
done
```

> TASSEL throws this error when sorting vcf files `ERROR net.maizegenetics.dna.map.PositionListBuilder - validateOrdering: Position	Chr:SCAF_100	Pos:79721	InsertionPos:0	Name:SSCAF_100_79721	Variants:A/C	MAF:NaN	Ref:A and Position	Chr:SCAF_99	Pos:109272	InsertionPos:0	Name:SSCAF_99_109272	Variants:T/A	MAF:NaN	Ref:T out of order`. However, I think it's just a warning showing which positions were in the wrong position. I'm able to load the sorted vcf file and transform it into hapmap format without problems. Also, no SNP is lost when sorting the file.

Then, I converted the file with SV calls for all NAM founders into hapmap with the following commands:

```bash
# go to project folder
cd ~/projects/sv_nams

# for explanation on how to use the script...
python scripts/variants2hapmap.py -h

# convert SVs vcf to hmp
python scripts/variants2hapmap.py data/NAM-structural-variations-v3.0.txt data/NAM_founders_SVs.not-sorted.hmp.txt

# sort hmp file
run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile data/NAM_founders_SVs.not-sorted.hmp.txt -outputFile data/NAM_founders_SVs.sorted.hmp.txt -fileType Hapmap
# convert to diploid format
run_pipeline.pl -Xmx10g -importGuess data/NAM_founders_SVs.sorted.hmp.txt -export data/NAM_founders_SVs.hmp.txt -exportType HapmapDiploid
```

Additional information about the SV is displayed on its ID, since hapmap format doesn't have fields available for adding such information. For example, the ID `del.chr1.51711.71809` on the first column of the hapmap file means that the SV is a deletion on chr1 that starts at 51,711 and ends at 71,809. The second column will also have the chromosome location for that deletion, but the third column will contain the **midpoint position** for that SV (i.e. 61,760). These two columns will always be the coordinates according to the reference genome. Although there will be somewhat redundant information on IDs of most SVs (like DELs, DUPs, INSs, and INVs), the ID will contain very important info about translocations, as it will show the respective location of the TRA in the **non-reference chromosome**.

Importantly, any SV called as heterozygous in the VCF file (i.e. `0/1`) was considered as **not** having a SV, therefore they were coded as `AA`.


### Identifying SNPs that are within the boundaries of a SV

SNPs that are found inside deletions are problematic, because they will have segregation issues when you compare multiple lines that have or not that SV. Thus, I wrote `scripts/generate_SV_bed.py` to create a BED file with start and end positions of deletions smaller than 100kb. I set up a 100 kb threshold because there were some extremely large deletions (>100 Mb) that would make me remove nearly all SNPs in this step and more than ~95% of the deletions were within that range. Importantly, each NAM population will be filtered separately.

> Translocations can cause SNP segregation issues as well. However, dealing with translocations is even more complicated, especially for SV projection and downstream GWAS, and we will ignore them here.

```bash
# go to project folder
cd ~/projects/sv_nams

# column numbers corresponding to NAM parents range from 13 to 37
for i in {13..37}; do
  # get NAM name
  NAM=$(head -n 1 data/NAM_founders_SVs.hmp.txt | cut -f $i)
  # create SV files for each NAM cross
  cut -f 1-12,$i data/NAM_founders_SVs.hmp.txt > data/NAM_founders_SVs_B73x$NAM.hmp.txt
  # create bed file with deletion boundaries for each population
  python scripts/generate_SV_bed.py data/NAM_founders_SVs_B73x$NAM.hmp.txt data/tmp/SNPs_to_remove_B73x$NAM.bed
  # # remove bed file's header for TASSEL compatibility
  # sed -i 1d data/tmp/SNPs_to_remove_B73x$NAM.bed
done

# after running the above, I noticed that 3 NAM parent names were slighlty different
# from the parental names in the GBS data (upper and lower case problme)
# so, I decided to change names of crosses now to avoid mismatches downstream
mv data/NAM_founders_SVs_B73xHP301.hmp.txt data/NAM_founders_SVs_B73xHp301.hmp.txt
mv data/NAM_founders_SVs_B73xIL14H.hmp.txt data/NAM_founders_SVs_B73xIl14H.hmp.txt
mv data/NAM_founders_SVs_B73xOh7b.hmp.txt data/NAM_founders_SVs_B73xOh7B.hmp.txt

mv data/tmp/SNPs_to_remove_B73xHP301.bed data/tmp/SNPs_to_remove_B73xHp301.bed
mv data/tmp/SNPs_to_remove_B73xIL14H.bed data/tmp/SNPs_to_remove_B73xIl14H.bed
mv data/tmp/SNPs_to_remove_B73xOh7b.bed data/tmp/SNPs_to_remove_B73xOh7B.bed
```


## GBS data

### Creating hapmap files for each NAM population and removing SNPs within SVs

The GBS file in VCF format contains information about millions of SNPs in all NAM lines (+5000), making it a very large file. Parsing it as it is would take way too much time. Thus, before doing any filtering, I split the VCF file by chromosome so I can run the filtering process for each chromosome in parallel.

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

# split vcf file by chromosome
module load parallel
for i in {1..10}; do
  sem -j +0 "grep -w '^chr$i' data/GBS-output/populations.snps.vcf >> data/GBS-output/tmp/NAM_rils_SNPs.chr$i.vcf"
done
sem --wait
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
      # transform line of the file into multiple lines so that vcftools recognize 1 genotype to keep per line
      echo $line |  tr "\t" "\n" | tr "," "\n" > data/GBS-output/tmp/$cross/genotypes_to_keep.txt
      # use vcftools to filter a vcf file
      vcftools --vcf data/GBS-output/tmp/NAM_rils_SNPs.$chr.vcf \
               --keep data/GBS-output/tmp/$cross/genotypes_to_keep.txt \
               --exclude-bed data/tmp/SNPs_to_remove_$cross.bed \
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
  echo "Rscript ~/projects/sv_nams/scripts/collapse_GBS_markers.R $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.hmp.txt $cross"
  done
done > ~/projects/sv_nams/scripts/commands_for_collapse-GBS-SNPs.txt

qsub ~/projects/sv_nams/scripts/collapse_GBS_markers.sh
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

# check number of SNPs per pop
wc -l B73*/*not-in-SVs.not-imputed.hmp.txt
# ~2.7M (with slightly different number per population)
```



### Overlaying resequencing data into parental GBS data

The NAM parents were genotyped by both resequencing and GBS. Thus, it is possible that some SNP calls disagree between the two methods. To remove such SNPs, I overlayed the resequencing data into the GBS data and turned a SNP call into `NN` if there was a disagreement using `scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R`. This script produces a hapmap for each cross containing **only the parental** data for that cross. Also, note that a number of SNPs in GBS data are not represented in the resequecing data, which can decrease even more the number of SNPs used for projection.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  echo "Rscript ~/projects/sv_nams/scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt ~/projects/sv_nams/data/tmp/NAM_founders_SNPs.chr1.hmp.txt $cross $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt"
done > ~/projects/sv_nams/scripts/commands_for_overlay-parental-SNPs.txt

qsub ~/projects/sv_nams/scripts/overlay_reseq-parental-SNPs_onto_GBS-data.sh

# check number of SNPs per pop
wc -l B73*/*not-in-SVs.reseq-overlay.hmp.txt
# ~1M (with slightly different number per population)
```

> Note: Apparently, the parent Tzi8 doesn't have gbs data. Therefore, I used the entire resequencing data for that parent.



### Select best GBS markers

GBS data contain a lot of missing data and also a lot of redundant information (many SNPs tightly linked). Besides increasing computation time with such big dataset, I also found in my preliminary analysis that using this raw GBS data had a strong negative impact on projections. Thus, we decided to filter this dataset by selecting only polymorphic SNPs, SNPs present in at least 30% of RILs, using a sliding window approach to remove incorrect calls (see [Huang et al, Genome Research, 2009](https://genome.cshlp.org/content/19/6/1068.abstract)) and removing SNPs with allele frequency < 0.4 or > 0.6. To do that, I ran `scripts/select_best_SNPs_per_pop.R`.

```bash
# go to data folder
cd ~/projects/sv_nams/data/GBS-output/tmp/

# # filter SNPs
# for cross in $(ls -d B73x*); do
#   Rscript ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.R $cross $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt ~/projects/sv_nams/analysis/qc/filter_best_SNPs --max_missing=0.3 --window_size=15 --window_step=1 --min_snps_per_window=5
# done

qsub ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.sh

# check how many SNPs remained
wc -l B73x*/*.not-imputed.best-markers.hmp.txt
# range from 13k to 50k
```

**In summary**, from initial ~2.7 million SNPs for each population (after collapsing duplicated markers), ~1 million SNPs were present in the parental resequecing data, and the number of best SNPs selected varied from ~13k to ~50k depending on the population (mean of ~32k SNPs and median of ~31k).

> Note: The reason why there is `not-imputed` in the filename is because during preliminary tests I tried to impute GBS SNPs using FSFHap from TASSEL to decrease the number of missing data. But after selecting the best markers, I found the imputing SNPs didn't reduce much the missing data and it was actually causing some troubles later during SV projection. So we decided not to impute SNPs at this stage.



### Summary

**Before filtering**

I wrote `scripts/summary_raw_gbs.R` to plot some basic statistics of the raw GBS SNPs (after collapsing duplicates) such as total number of RILs, total number of SNPs, percentage of missing data and percentage of polymorphic SNPs for each population. It also plots the distribution of missing data per population, missing data per RIL, and missing data per SNP.

| Average number SNPs | Average missing data | Average polymorphic |
| ------------------- | -------------------- | ------------------- |
| 2,730,043           | 85.48%               | 29.32%              |

```bash
cd ~/projects/sv_nams/

# create new folder to store qc results
mkdir -p analysis/qc

# summarize data
Rscript scripts/summary_raw_gbs.R data/GBS-output/tmp analysis/qc/raw_gbs
```

In order to visualize how the markers are distributed along the chromosomes, I ploted karyotypes of 3 random RILs for each population with `scripts/plot_ril_karyotypes.R`. I used the chromosome coordinates from `analysis/qc/B73_RefGen_V4_chrm_info.txt` and the centromere positions from `analysis/qc/centromeres_Schneider-2016-pnas_v4.bed`, which are from B73 v4 reference genome.


```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# create karyotypes
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/raw-gbs $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=random --parents_in_data=TRUE --overlay_reseq=TRUE
done

# # run extra qc with TASSEL
# for cross in $(ls -d B73x*); do
#   run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/raw_gbs/NAM_rils_SNPs_raw-gbs_$cross\_OverallSummary
#   (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/raw_gbs/NAM_rils_SNPs_raw-gbs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/raw_gbs/missing_data_raw-gbs.txt
# done
```


**After filtering**

When I filtered the raw gbs, `scripts/select_best_SNPs_per_pop.R` already produces some summary data. Thus, I just needed to generate the karyotypes for the same RILs used before filtering.

| Average number SNPs | Average missing data | Average polymorphic |
| ------------------- | -------------------- | ------------------- |
| 32,190              | 21.88%               | 100%                |

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# create karyotypes
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/raw-gbs/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/best-markers $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=$rils --parents_in_data=TRUE --overlay_reseq=FALSE
done

# # run extra qc with TASSEL
# for cross in $(ls -d B73x*); do
#   run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/best-markers/NAM_rils_SNPs_best-markers_$cross\_OverallSummary
#   (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/best-markers/NAM_rils_SNPs_best-markers_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/best-markers/missing_data_best-markers.txt
# done
```

Finally, I wrote `scripts/correct_SNP-names_rils.R` to make sure SNPs from RILs have the same name as the ones from parents. This script will generate the file `NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt`, which will be used later when merging SNPs with SVs.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/correct_SNP-names_rils.R $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt
done
```


## Merge SNPs with SVs

Now that I have selected which SNPs will be used to anchor projections, it's time to merge this SNP data with the SV calls into the same file for each population. I wrote `scripts/merge_SVs_and_SNPs.R` to do this task, and it generates a hapmap for the two parents of a cross (both SNP and SV calls come from resequecing), and another hapmap with the RILs of that cross (SNP calls come from GBS, and SV calls were set to `NN` since that's what will be projected later).

```bash
module load R/3.6.0

cd ~/projects/sv_nams/data/GBS-output/tmp/

# merge svs
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/merge_SVs_and_SNPs.R ~/projects/sv_nams/data/NAM_founders_SVs_$cross.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt
done

# check that parents and rils have the same number of markers
for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt
done
```

> Note: the B73 parent should not have any SV present, because all SVs were called against the B73 reference genome. However, since in practice there might be some miscalls, Arun used B73Ab10 as a negative control. However, this parent has a lot of missing data. This missing data can be problematic for projection. Since it's expected that B73Ab10 has very few SVs, I converted all `NN` into `AA` and left all non-missing calls as they were.




## Projection

All the projections of parental SVs into RILs will be done for each NAM population separately using the FILLIN plugin from TASSEL 5. In short, this program generate haplotypes for each parent of a cross (`-FILLINFindHaplotypesPlugin`), and then it imputes missing data in the RILs (i.e. SVs) using those parental haplotypes (`-FILLINImputationPlugin`).

But, before that, I need to make sure the files are sorted:

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

# make sure number of rows of parental and RIL data matches in each population
for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt
done
```

After some preliminary tests, I found that using a haplotype block size of 2,000 sites (`-hapSize 2000`) gives the best projection rate without compromising accuracy too much. Also, since I'm creating haplotypes for each parent individually, I need to set the minimum number of taxa to create a haplotype to 1 (`-minTaxa 1`). Importantly, the option `-hybNN` should be turned to `false`, otherwise the algorithm will combine haplotypes in recombination breakpoints if it thinks that region is heterozygous. If this option is not turned off, some projected regions will be messy and can end up with more than two alleles. It it's off, nothing is projected for that region.

```bash
# create new folder
mkdir ~/projects/sv_nams/analysis/projection

# create haplotypes from parents
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -FILLINFindHaplotypesPlugin  -hmp ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt  -o ~/projects/sv_nams/analysis/projection/donors_$cross -hapSize 2000 -minTaxa 1
done

# impute ril genotypes based on parental haplotypes
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -FILLINImputationPlugin -hmp ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -d ~/projects/sv_nams/analysis/projection/donors_$cross -o ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -hapSize 2000 -accuracy -hybNN false
done

# convert projected hapmap to diploid format
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -export ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -exportType HapmapDiploid
done
```

After projection, I wrote `scripts/count_projected_SVs.R` to compute how many SVs were projected per population and make summary plots about projections. I also wrote `scripts/plot_ril_karyotypes_SVs.R` to plot karyotypes of few RILs showing from which parent the SVs come from. The RILs selected were the same ones previously used to plot karyotypes of raw GBS data and best selected SNP markers.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# calculate amount of projected SVs
Rscript ~/projects/sv_nams/scripts/count_projected_SVs.R ~/projects/sv_nams/data ~/projects/sv_nams/analysis/projection

# plot karyotypes of SVs present in each parent of a cross
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/projection ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt --rils=$rils --expected-SVs=TRUE
done

# now plot karyotypes with projected SVs for few RILs of each cross
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/best-markers/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/projection ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt --rils=$rils --expected-SVs=FALSE
done

# additional QC about missing data
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/projection/missing_data_best-markers_after_SV-projection.txt
done
# average missing data 0.09
```

The average percentage of projected SVs across all populations was **87%** (~180k SVs) with average accuracy of **93%**. The amount projected was a bit higher (88%) if considering only polymorphic SVs between the two parents of a cross. Only five crosses had projection rate below 75% (B73xCML322, B73xCML333, B73xOh7B, B73xP39, and B73xTzi8).



## Merge all projected SVs of each population in one file

The final file that will be sent to Jianming's group for GWAS will be a hapmap file **only with SVs** with all RILs of all NAM populations. I wrote `scripts/merge_SVs_after_projection.R` to filter and merge all SVs in one file. After combining multiple files, I have to correct the `alleles` column on the final hapmap to make sure that all possible alleles are represented in that column.

```bash
cd ~/projects/sv_nams/analysis/projection

# make sure files are sorted
for file in NAM_rils_projected-SVs-only.B73x*; do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile $file -outputFile $file -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess $file -export $file -exportType HapmapDiploid
done

# merge hapmaps
Rscript ~/projects/sv_nams/scripts/merge_SVs_after_projection.R

# make sure final file is sorted and let TASSEL correct the alleles' column
run_pipeline.pl -Xmx10g \
                -SortGenotypeFilePlugin \
                -inputFile ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.hmp.txt \
                -outputFile ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -fileType Hapmap
run_pipeline.pl -Xmx10g \
                -importGuess ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -export ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -exportType HapmapDiploid
```



## Upload final hapmap to Cyverse

Lastly, I uploaded the final hapmap file `~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt` to the shared folder on Cyverse.

```bash
# go to data folder of the project
cd ~/projects/sv_nams/analysis/projection

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/Misc
# check if files match what Arun described
ils
# upload data
iput -K NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt
# exit iRods
iexit full
```
