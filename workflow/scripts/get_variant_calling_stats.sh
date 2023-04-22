#!/bin/bash

# Set the VCF file paths
unfiltered_vcf="output/6_variant_filtering/c_gatk_variants_merged/unfiltered.vcf.gz"
filtered_vcf="output/6_variant_filtering/c_gatk_variants_merged/filtered.vcf.gz"

# Extract the counts of SNPs and indels for each VCF file
snps_raw=$(bcftools view -v snps "$unfiltered_vcf" | grep -v "^#" | wc -l)
snps_pass=$(bcftools view -v snps "$filtered_vcf" | grep -v "^#" | wc -l)
indels_raw=$(bcftools view -v indels "$unfiltered_vcf" | grep -v "^#" | wc -l)
indels_pass=$(bcftools view -v indels "$filtered_vcf" | grep -v "^#" | wc -l)

# Save the results to a file
echo -e "snps_raw\tsnps_pass\tindels_raw\tindels_pass" >stats_4_variant_calling.tsv
echo -e "${snps_raw}\t${snps_pass}\t${indels_raw}\t${indels_pass}" >>stats_4_variant_calling.tsv
