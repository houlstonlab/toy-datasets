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

# Download first FASTQ reads
tail -n +2 resources/wes_sarek_input.csv | while IFS=, read -r patient sex status sample lane fastq_1 fastq_2
do
  # Download the fastq_1 and fastq_2 files
  wget -N -P ${TMP}/samples/ "$fastq_1"
  wget -N -P ${TMP}/samples/ "$fastq_2"

  sleep 5
done

# create sarek input file with local samples
cat resources/wes_sarek_input.csv |  head -n 1 > ${TMP}/input.csv
cat resources/wes_sarek_input.csv | tail -n +2 | awk -F, 'BEGIN {OFS=","} {split($6, a, "/"); split($7, b, "/"); $6="samples/"a[length(a)]; $7="samples/"b[length(b)]; print}' >> ${TMP}/input.csv

# run sarek pipeline
cp resources/wild.config ${TMP}
cd ${TMP}

# Install nextflow (Run once)
module load java/jdk15.0.1
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install nextflow
pip install nf-core

# Run nextflow
nextflow run nf-core/sarek -r 3.2.3 \
	-profile singularity \
	-c wild.config \
	-resume \
	--input wes_sarek_local_input.csv \
	--outdir . \
	--wes true \
	--vep_loftee true \
	--vep_spliceai true \
	--tools vep,haplotypecaller \
	--vep_cache_version 111 \
	--genome 'GATK.GRCh38' \
	--igenomes_base ${REFS}/iGenomes/ \
	--vep_cache ${SCRATCH}/mahmed03/cache/vep/ \
	--spliceai_snv ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz \
	--spliceai_snv_tbi ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz.tbi \
	--spliceai_indel ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
	--spliceai_indel_tbi ${SCRATCH}/mahmed03/cache/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz.tbi \
	--joint_germline

# copy annotated VCF