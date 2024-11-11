all: setup \
	vcf \
	vcf-references \
	wes-sarek \
	list-all

setup:
	@echo "Setting up directories"
	@[ -z "${TOY_DATASETS}" ] && mkdir -p toy-datasets/ && export TOY_DATASETS=toy-datasets && echo "TOY_DATASETS is set to toy-datasets/" || echo "TOY_DATASETS is defined as ${TOY_DATASETS}"

vcf: scripts/create_vcf.sh resources/vcf_urls.txt
	@echo "Creating a VCF file"
	sbatch scripts/create_vcf.sh

vcf-references: scripts/create_vcf_references.sh
	@echo "Creating vcf files for reference datasets"
	sbatch scripts/create_vcf_references.sh

list-all:
	@echo "Listing all files"
	ls -d ${TOY_DATASETS}/* | sort > available-files.txt
