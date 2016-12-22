#!/usr/bin/env Rscript

rm(list = ls())

library(plyr)
library(dplyr)
library(doMC)
library(data.table)

ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores / 2))

setwd("/xchip/scarter/dmccabe/germline/tcga_luad/out")
files <- list.files("./", pattern = ".hets.filtered$", recursive = T)

aaply(files, 1, function(f) {
  callstats <- fread(f)
  callstats <- callstats[contig %in% as.character(1:22)]
  callstats[, contig := as.numeric(contig)]
  
  setcolorder(callstats, c("contig", "position", "ref_allele", "alt_allele", "t_ref_count", "t_alt_count"))
  colnames(callstats) <- c("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele1", "i_t_ref_count", "i_t_alt_count")
  out_pathfile <- paste0(f, ".ghet")
  
  if (!file.exists(out_pathfile)) {
    message("Writing ", out_pathfile)
    write.table(callstats, file = out_pathfile, sep = "\t", row.names = F, quote = F)
  } else {
    message(out_pathfile, " exists")
  }
  
  return(0)
}, .parallel = T)

warnings()
