#!/bin/bash
#PBS -l walltime=48:00:00,nodes=1:ppn=10,mem=30gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis
#PBS -V
#PBS -N remove_SNPs_within_SVs
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams
# run script
python scripts/remove_SNPs_within_SVs.py data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf data/NAM_sv_sniffles_v1.vcf data/B73v5.NAM-gatk-snps.not-in-SVs.vcf
