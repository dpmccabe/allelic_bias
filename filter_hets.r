#!/usr/bin/env Rscript

rm(list = ls())

library(io)
library(plyr)
library(dplyr)
library(argparser)
library(doMC)
library(data.table)

ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores / 2))

ref_snpv <- readRDS("/xchip/scarter/dshih/data/ucsc/hg19/snp146Common-filtered.rds")

setwd("/xchip/scarter/dmccabe/germline/tcga_luad/out")
cov_files <- list.files("./", pattern = "hets$", recursive = T)

aaply(cov_files, 1, function(f) {
  message(f)
  outfile <- paste0(f, ".filtered")
  
  if (!file.exists(outfile)) {
    input <- fread(f, stringsAsFactors = F)
    input[, mut := sprintf("%s:%d%s>%s", CONTIG, POSITION, REF_NUCLEOTIDE, ALT_NUCLEOTIDE)]
    setkey(input, mut)

    input <- input[mut %in% ref_snpv]
    
    message("Writing ", outfile)
    input <- input[, .(CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE)]
    colnames(input) <- c("contig", "position", "t_ref_count", "t_alt_count", "ref_allele", "alt_allele")
    write.table(input[, .(contig, position, t_ref_count, t_alt_count, ref_allele, alt_allele)], file = outfile,
                sep = "\t", row.names = F, quote = F)
  }
  
  return(0)
}, .parallel = wpar)
