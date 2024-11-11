#!/bin/bash

#SBATCH -J vcf
#SBATCH -o vcf.out
#SBATCH -e vcf.err
#SBATCH -p compute
#SBATCH -c 8
#SBATCH -t 24:00:00

# Create output directories
OUT="${TOY_DATASET}/vcf"
TMP="${TOY_DATASET}/vcf/tmp"

mkdir -p ${OUT} ${TMP}

# Load required modules
module load BCFtools

# Loop over the URLs in the file and download the VCF files
while read url; do
    curl -C - -o ${TMP}/$(basename $url) $url
    tabix -f ${TMP}/$(basename $url)
done < resources/vcf_urls.txt

# Concatenate the VCF files
bcftools concat --naive ${TMP}/ALL.*.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
bcftools view --threads ${SLURM_CPUS_PER_TASK} -Oz -o ${TMP}/all.variants.vcf.gz

tabix -f ${TMP}/all.variants.vcf.gz

# Create a subset of 100 samples and 1000 variants
bcftools query -l ${TMP}/all.variants.vcf.gz | shuf -n 100 > ${TMP}/subset.samples.txt
bcftools query -f '%ID\n' ${TMP}/all.variants.vcf.gz | shuf -n 1000 > ${TMP}/subset.variants.txt

# Create the final VCF file
bcftools view -S ${TMP}/subset.samples.txt ${TMP}/all.variants.vcf.gz | \
bcftools view -i ID==@${TMP}/subset.variants.txt | \
bcftools view --threads ${SLURM_CPUS_PER_TASK} -Oz -o ${OUT}/1kg.variants.vcf.gz

tabix -f ${OUT}/1kg.variants.vcf.gz
