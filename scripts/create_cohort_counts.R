#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
out_dir <- paste0(args[1], args[2])
if (!dir.exists(out_dir)) dir.create(out_dir)

# Creat a tibble with the cohort information
cohorts_info <- tibble::tibble(
  cohort = c('pheno', 'general'),
  type   = c('cases', 'controls'),
  size   = c(100, 1000)
)

readr::write_csv(
  cohorts_info, 
  file.path(out_dir, 'cohorts_info.csv')
)

# Create a tibble with the gene-level counts for the general cohort
general.aggregate <- tibble::tibble(
  gene = paste0('gene', 1:6),
  nvar = c(5, 5, 10, 10, 5, 10),
  ac = c(5, 5, 110, 110, 5, 110),
  an = rep(2000, 6),
  af = c(.0025, .0025, .0055, .0055, .0025, .0055),
  nhom = c(0, 0, 1, 2, 0, 1)
)

readr::write_tsv(
  general.aggregate, 
  file.path(out_dir, 'general.categoryA.aggregate.tsv')
)

# Reverse the order of the genes
general.aggregate$gene <- general.aggregate$gene[nrow(general.aggregate):1]

readr::write_tsv(
  general.aggregate, 
  file.path(out_dir, 'general.categoryB.aggregate.tsv')
)

# Create a tibble with the gene-level counts for the pheno cohort
pheno.aggregate <- tibble::tibble(
  gene = paste0('gene', 1:5),
  het = 5:1,
  hom = c(1, 0, 0, 0, 0),
  ch = c(0, 1, 1, 0, 0)
)

readr::write_tsv(
  pheno.aggregate, 
  file.path(out_dir, 'pheno.categoryA.aggregate.tsv')
)

# Reverse the order of the genes
pheno.aggregate$gene <- pheno.aggregate$gene[nrow(pheno.aggregate):1]

readr::write_tsv(
  pheno.aggregate, 
  file.path(out_dir, 'pheno.categoryB.aggregate.tsv')
)
