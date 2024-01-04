#!/bin/bash

# Script for doing the SNP filters.
# Including:
# Step 1: Indel and PASS filter
# Step 2: QUAL filter, MQ filter, QD filter (after this step, renamed to "fitler1")
# Step 3: Missingness filter, biallelic filter
# Step 4: Removing small scaffolds (<100kb)
# Step 5: Removing sex-linked snps
# Step 6: Additional HWE, MAF filter
# Step 7: LD pruning

# The resulting file through Step 7 is for population structure analyses. Specific anaylses may use different filters in the end (specifically referring to Step 6 & 7; e.g. need to exclude the hwe filter for selection analyses, need to exclude maf filter for demographic history).

# Make sure edit_sampleID.sh is done to vcf beforehand if needed
# Also need filter.hets.R file in directory for the HWE filter to work
# And need list of scaffolds (for removing small scaffold step) and list of X and Y scaffolds. 

# Set up these before running
#variant_prefix=NAR_allvariantcalls_pp.ID
snps_prefix=killerwhale2_snps.ID
bcftools_path=/home/cailak/programs/bcftools-1.9
gatk_jar_path=/home/cailak/programs/gatk-4.1.9.0
ref_genome_path=/home/cailak/killer_whales/ref_genome
ref_genome_prefix=KW_GCA_937001465.1_mOrcOrc1.1.scafname
min_scaf_list=mOrcOrc_scafmin_100kb #list of scaffolds to keep
#xchr_scaf_list=final_scaffolds_Xlinked.txt #list of scaffolds in X chromosome
#ychr_scaf_list=final_scaffolds_Ylinked.txt #list of scaffolds in Y chromosome



######## Step 1: Remove indels and SNPs that don't have the "PASS" filter. Start with vcf.gz file (with corresponding tbi)
# First count snps for starting count
#vcftools --gzvcf $variant_prefix.vcf.gz --out $variant_prefix.snpcount

# Remove indels (looks like the platypus remove indel filter didn't remove all of them).
#vcftools --gzvcf $variant_prefix.vcf.gz --remove-indels --recode --recode-INFO-all --out $snps_prefix

#mv $snps_prefix.recode.vcf $snps_prefix.vcf

# Zip & index if not zipped already
#bgzip $snps_prefix.vcf
#tabix -p vcf $snps_prefix.vcf.gz

# Count snps
vcftools --gzvcf $snps_prefix.vcf.gz --out $snps_prefix.snpcount

# Filter out non-PASS sites
$bcftools_path/bcftools view -f PASS $snps_prefix.vcf.gz -o $snps_prefix.PASS.vcf

# Zip & index
bgzip $snps_prefix.PASS.vcf
tabix -p vcf $snps_prefix.PASS.vcf.gz

# Count snps
vcftools --gzvcf $snps_prefix.PASS.vcf.gz --out $snps_prefix.PASS.snpcount



######## Step 2: Next filtering round
# Need reference genome .dict file in same directory as reference fasta
#java -jar /home/cailak/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar CreateSequenceDictionary -R $ref_genome_path/$ref_genome_prefix.fasta -O $ref_genome_path/$ref_genome_prefix.dict 

# Filter out QUAL < 50
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar VariantFiltration -R $ref_genome_path/$ref_genome_prefix.fasta -V $snps_prefix.PASS.vcf.gz -O $snps_prefix.PASS.QUAL.vcf.gz --filter-name "QUALlt50" --filter-expression "QUAL < 50" 

# Count snps
vcftools --gzvcf $snps_prefix.PASS.QUAL.vcf.gz --out $snps_prefix.PASS.QUAL.snpcount

# Filter out MQ < 40
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar SelectVariants -R $ref_genome_path/$ref_genome_prefix.fasta -V $snps_prefix.PASS.QUAL.vcf.gz -O $snps_prefix.PASS.QUAL.MQ.vcf.gz -select "MQ >= 40.0"

# Count snps
vcftools --gzvcf $snps_prefix.PASS.QUAL.MQ.vcf.gz --out $snps_prefix.PASS.QUAL.MQ.snpcount

# Filter out QD < 4
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar SelectVariants -R $ref_genome_path/$ref_genome_prefix.fasta -V $snps_prefix.PASS.QUAL.MQ.vcf.gz -O $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz -select "QD >= 4.0"

# Count snps
vcftools --gzvcf $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz --out $snps_prefix.PASS.QUAL.MQ.QD.snpcount

# Rename with "filter1" and clean up files
mv $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz $snps_prefix.filter1.vcf.gz
mv $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz.tbi $snps_prefix.filter1.vcf.gz.tbi

rm $snps_prefix.PASS.QUAL.vcf.gz $snps_prefix.PASS.QUAL.vcf.gz.tbi $snps_prefix.PASS.QUAL.MQ.vcf.gz $snps_prefix.PASS.QUAL.MQ.vcf.gz.tbi $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz.tbi



######## Step 3: Missingness & bi-allelic filter. These two steps should automatically create the log files for snp counts

# Filter out snps with high missingness (in vcftools 1=no missing, 0=all missing)
vcftools --gzvcf $snps_prefix.filter1.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out $snps_prefix.filter1.miss

# Remove non-biallelic sites
vcftools --vcf $snps_prefix.filter1.miss.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel

mv $snps_prefix.filter1.miss.biallel.recode.vcf $snps_prefix.filter1.miss.biallel.vcf

# Zip & index
bgzip $snps_prefix.filter1.miss.biallel.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.vcf.gz



######## Step 4: Remove scaffolds less than 100kb in length. Needs input list of scaffold names to KEEP. I chose 100kb based on my dataset, but should be checked/adjusted for other data (example, if this filter loses a very large chunk of the genome, may need to use a smaller scaffold filter like 10kb).

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $min_scaf_list | paste -s -d, - > scaf_list_line

# Set up list
list=`cat scaf_list_line`

# Filter vcf for these scaffolds
$bcftools_path/bcftools filter --regions $list $snps_prefix.filter1.miss.biallel.vcf.gz > $snps_prefix.filter1.miss.biallel.min100kb.vcf

# Can make a list of scaffolds to double check
grep -v "^#" $snps_prefix.filter1.miss.biallel.min100kb.vcf | cut -f1 | sort | uniq > filtered_contig_list_check.txt

# Zip & index
bgzip $snps_prefix.filter1.miss.biallel.min100kb.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.vcf.gz

# Count snps
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.vcf.gz --out $snps_prefix.filter1.miss.biallel.min100kb.snpcount



######## Step 5: Filter out sex-linked SNPS to create autosomal dataset

# Do the filter to remove snps on X and Y
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.vcf.gz --not-chr OW443365.1 --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes
mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf

# Zip and index
bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf.gz

# Can make a list of scaffolds to double check
#grep -v "^#" $snps.X.vcf | cut -f1 | sort | uniq > X_scaffold_check.txt



######## Step 6: Additional HWE, MAF filter
# Remove SNPs out of Hardy-Weinberg Equilibrium
# First need to create the out.hwe file
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf.gz --hardy

# Run filter.hets.R script to create list of snps with het > 0.6
R --vanilla < filter.hets.R
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf.gz --exclude-positions snps_het06_filter --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf
bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz

# Minor-allele frequency filter
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz --maf 0.05 --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf

bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz



######## Step 7: LD pruning

# Create list of snps to prune
# --double-id b/c multiple "_" in sample name
plink --vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --allow-extra-chr --indep-pairwise 50 5 0.8 --set-missing-var-ids @:#\$1,\$2 --out pruned.50window.r0.8 --double-id

# Preparing the prune-in file to filter vcf
sed 's/:/\t/g' pruned.50window.r0.8.prune.in > temp.in
sed 's/...$//' temp.in > temp2.in
mv temp2.in list.pruned.50window.r0.8.prune.in
rm temp.in

# Prune the snps
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --positions list.pruned.50window.r0.8.prune.in --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf

bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf.gz