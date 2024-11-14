#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

# setwd('/data/rds/DGE/DUDGE/MOPOPGEN/toy-datasets/wes-sarek/tmp/')
ids <- readr::read_lines(input)

x = pedtools::nuclearPed(nch = 2)
x$ID <- ids[1: 4]
d1 <- as.data.frame(x)
d1$aff <- ifelse(x$ID %in% pedtools::leaves(x), 2, 1)
d1$famid <- 'FAM_01'

x = pedtools::linearPed(n = 2)

x$ID <- ids[5: 9]
d2 <- as.data.frame(x)
d2$aff <- ifelse(x$ID %in% c(pedtools::father(x, pedtools::leaves(x)), pedtools::leaves(x)), 2, 1)
d2$famid <- 'FAM_02'

d <- dplyr::bind_rows(d1, d2)
readr::write_tsv(d, output)
