#!/bin/bash

# Script for editing sample IDs in the vcf 

snps=killerwhale2_allvariantcalls

# Zip and index vcf if not done already
#bgzip $snps.vcf
tabix -p vcf $snps.vcf.gz

# Make list of current sample ID names
vcf-query -l $snps.vcf.gz > vcf_sampleID

# Shorten sample IDs by removing text up to the last "/"
awk '{print $NF}' FS=/ vcf_sampleID > vcf_sampleID_shortened

# Also want to remove characters in sample ID that have "_S43" or something like that which was added from the sequencing process, unrelated to sample ID
# Use sed to remove "_S" and all the characters afterwards
sed 's/\_S.*//' vcf_sampleID_shortened > vcf_sampleID_shortened2

# Modify header in vcf file to change sample ID to the shortened ID list
bcftools reheader -s vcf_sampleID_shortened2 $snps.vcf.gz -o $snps.ID.vcf.gz

# Index new vcf file
tabix -p vcf $snps.ID.vcf.gz
