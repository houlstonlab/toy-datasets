all: setup \
	vcf \
	vcf-references \
	wes-sarek \
	cohort-counts \
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

wes-sarek: scripts/create_wes_sarek.sh resources/wes_sarek_input.csv
	@echo "Creating WES Sarek files"
	sbatch scripts/create_wes_sarek.sh

cohort-counts: scripts/create_cohort_counts.R
	@echo "Creating cohort counts"
	Rscript scripts/create_cohort_counts.R ${TOY_DATASETS} cohort-counts

list-all:
	@echo "Listing all files"
	ls -d ${TOY_DATASETS}/* | sort > available-files.txt
