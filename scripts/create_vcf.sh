#!/bin/bash

# URL of the VCF file
url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
output_file="ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
output_dir="vcf"

# Use curl to download the file
curl -o $output_dir/$output_file $url/$output_file
curl -o $output_dir/$output_file.tbi $url/$output_file.tbi
