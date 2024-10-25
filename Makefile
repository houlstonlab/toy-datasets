all: setup \
	vcf \
	vep-annotated-vcf \
	list-all

setup:
	@echo "Setting up directories"
	mkdir -p vcf
	mkdir -p vep-annotated-vcf

vcf: 
	@echo "Creating a VCF file"
	sh scripts/create_vcf.sh

vep-annotated-vcf:
	@echo "Creating a VEP annotated VCF file"
	sh scripts/create_vep_annotated_vcf.sh

list-all:
	@echo "Listing all files"
	ls -d {vcf,vep-annotated-vcf}/* | sort > available-files.txt
