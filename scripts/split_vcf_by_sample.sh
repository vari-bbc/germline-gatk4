#!/bin/bash

set -e
set -u
set -o pipefail

module load bbc2/bcftools/bcftools-1.17


file="analysis/final/07_snpEff/all.merged.filt.PASS.snpEff.vcf.gz"
basenm=$(basename $file)

mkdir -p homo_alt
mkdir -p per_sample_vcf

for sample in `bcftools query -l $file`; do
  bcftools view --min-ac 2 --genotype hom -Ov -s $sample -o ./homo_alt/${basenm/.vcf*/.$sample.vcf} $file # extract sample and keep only sites with 2 ALT allele counts and is homozygous
  bcftools view --min-ac 1 -Ov -s $sample -o ./per_sample_vcf/${basenm/.vcf*/.$sample.vcf} $file # extract sample and keep only sites with at least one ALT allele count
done
