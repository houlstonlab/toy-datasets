#!/bin/bash

#SBATCH -J vcf-references
#SBATCH -o vcf-references.out
#SBATCH -e vcf-references.err
#SBATCH -p compute
#SBATCH -c 8
#SBATCH -t 24:00:00

# Create output directories
OUT="${TOY_DATASET}/vcf-references"
TMP="${TOY_DATASET}/vcf-references/tmp"

mkdir -p ${OUT} ${TMP}

GNOMAD="${SCRATCH}/mahmed03/cache/gnomad_v4_exomes"
CLINVAR="${REFS}/NCBI/ClinVar/GRCh38"

# Load the required modules
module load BCFtools

# Concatenate the GnomAD VCF files
bcftools concat --naive \
    ${GNOMAD}/gnomad.exomes.v4.1.sites.*.vcf.bgz \
    --threads ${SLURM_CPUS_PER_TASK} \
    -Oz -o ${TMP}/gnomad.v4.vcf.gz
tabix -f ${TMP}/gnomad.v4.vcf.gz

# Extract random variants from the GnomAD VCF
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${TMP}/gnomad.v4.vcf.gz > ${TMP}/gnomad.variants.txt
cat ${TMP}/gnomad.variants.txt | shuf -n 10000 > ${TMP}/random.variants.txt
cat  ${TMP}/random.variants.txt | sed 's/^chr//' > ${TMP}/random.variants.short.txt

# Subset the GnomAD VCF to the random variants
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${TMP}/gnomad.v4.vcf.gz | \
bcftools view -i ID==@${TMP}/random.variants.txt | \
bcftools view --threads ${SLURM_CPUS_PER_TASK} -Oz -o ${OUT}/gnomad.v4.vcf.gz
tabix -f ${OUT}/gnomad.v4.vcf.gz

# # Subset the ClinVar VCF to the random variants
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${CLINVAR}/clinvar_20200520.vcf.gz | \
bcftools view -i ID==@${TMP}/random.variants.short.txt | \
bcftools view --threads ${SLURM_CPUS_PER_TASK} -Oz -o ${OUT}/clinvar.20200520.vcf.gz
tabix -f ${OUT}/clinvar.20200520.vcf.gz
