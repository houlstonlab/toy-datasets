#!/bin/bash

#SBATCH -J wes-sarek
#SBATCH -o wes-sarek.out
#SBATCH -e wes-sarek.err
#SBATCH -p master-worker
#SBATCH -t 120:00:00

# Create output directories
OUT="${TOY_DATASET}/wes-sarek"
TMP="${TOY_DATASET}/wes-sarek/tmp"

mkdir -p ${OUT} ${TMP} ${TMP}/samples

# # Download first FASTQ reads
# tail -n +2 resources/wes_sarek_input.csv | while IFS=, read -r patient sex status sample lane fastq_1 fastq_2
# do
#   # Download the fastq_1 and fastq_2 files
#   wget -N -P ${TMP}/samples/ "$fastq_1"
#   wget -N -P ${TMP}/samples/ "$fastq_2"

#   sleep 5
# done

# # create sarek input file with local samples
# cat resources/wes_sarek_input.csv |  head -n 1 > ${TMP}/input.csv
# cat resources/wes_sarek_input.csv | tail -n +2 | awk -F, 'BEGIN {OFS=","} {split($6, a, "/"); split($7, b, "/"); $6="samples/"a[length(a)]; $7="samples/"b[length(b)]; print}' >> ${TMP}/input.csv

# run sarek pipeline
# cp resources/wild.config ${TMP}
cd ${TMP}

# # Install nextflow (Run once)
# module load java/jdk15.0.1
# python3 -m venv venv
# source venv/bin/activate
# pip install --upgrade pip
# pip install nextflow
# pip install nf-core

# # Run nextflow
# nextflow run nf-core/sarek -r 3.2.3 \
# 	-profile singularity \
# 	-c wild.config \
# 	-resume \
# 	--input wes_sarek_local_input.csv \
# 	--outdir . \
# 	--wes true \
# 	--vep_loftee true \
# 	--vep_spliceai true \
# 	--tools vep,haplotypecaller \
# 	--vep_cache_version 111 \
# 	--genome 'GATK.GRCh38' \
# 	--igenomes_base ${REFS}/iGenomes/ \
# 	--vep_cache ${SCRATCH}/mahmed03/cache/vep/ \
# 	--spliceai_snv ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz \
# 	--spliceai_snv_tbi ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz.tbi \
# 	--spliceai_indel ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
# 	--spliceai_indel_tbi ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz.tbi \
#     --cadd_indels ${REFS}/CADD/v1.6/prescored/GRCh38_v1.6/no_anno/gnomad.genomes.r3.0.indel.tsv.gz \
#     --cadd_indels_tbi ${REFS}/CADD/v1.6/prescored/GRCh38_v1.6/no_anno/gnomad.genomes.r3.0.indel.tsv.gz.tbi \
#     --cadd_wg_snvs ${REFS}/CADD/v1.6/prescored/GRCh38_v1.6/no_anno/whole_genome_SNVs.tsv.gz \
#     --cadd_wg_snvs_tbi ${REFS}/CADD/v1.6/prescored/GRCh38_v1.6/no_anno/whole_genome_SNVs.tsv.gz.tbi \
# 	--cadd_cache \
# 	--joint_germline

# Choose random samples and variants from a VCF file
module load BCFtools

VCF="annotation/haplotypecaller/recalibrated_joint_variant_calling/joint_germline_recalibrated_VEP.ann.vcf.gz"

# # Select samples and variants
# bcftools query -l ${VCF} | shuf -n 10 > subset.samples.txt
# bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,SpliceAI_pred_DS_AG:Float -f '%ID\t%IMPACT\t%CLIN_SIG\t%SpliceAI_pred_DS_AG\n' ${VCF} > annotations.tsv
# cat annotations.tsv | awk '$1 != "." && ($2 == "HIGH" ||  $3 == "pathogenic" || $4 > 0.5)' > subset.variants.txt

# # Subset the VCF file
# bcftools view -g het -S subset.samples.txt ${VCF} | \
# bcftools view -i ID==@subset.variants.txt | \
# bcftools view --threads 4 -Oz -o ${OUT}/pheno.variants.vcf.gz
# tabix -f ${OUT}/pheno.variants.vcf.gz

# # Create a list of cases
# bcftools query -l ${VCF} > ${OUT}/pheno.cases.txt

# echo '##INFO=<ID=CADD_PHRED,Number=1,Type=Integer,Description="CADD score">' > annots.hdr
# echo '##INFO=<ID=CADD_RAW,Number=1,Type=Integer,Description="CADD score">' >> annots.hdr
echo '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|CADD_RAW|CADD_PHRED">' > annots.hdr
# Create the file, zip, and index in the same step
(
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ|30|30\n' ${OUT}/pheno.variants.vcf.gz
) | bgzip -c > cadd.tab.gz
tabix -s 1 -b 2 -e 2 cadd.tab.gz

bcftools annotate \
    -x INFO/CSQ \
    -a cadd.tab.gz -h annots.hdr \
    -c CHROM,POS,REF,ALT,+CSQ ${OUT}/pheno.variants.vcf.gz \
    -Oz -o ${OUT}/pheno.variants.annotated.vcf.gz
tabix -f ${OUT}/pheno.variants.annotated.vcf.gz

# Create a PED file
# Rscript ${OUT}/pheno.cases.txt ${OUT}/pheno.families.ped