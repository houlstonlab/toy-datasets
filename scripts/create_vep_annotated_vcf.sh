#!/bin/bash

# URL of the VCF file
url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/"
output_file="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.recalibrated_variants.vcf.gz"
output_dir="vep-annotated-vcf"

# Use curl to download the file
curl -o $output_dir/$output_file $url/$output_file
curl -o $output_dir/$output_file.tbi $url/$output_file.tbi
