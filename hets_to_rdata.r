rm(list = ls())

.libPaths("/xchip/scarter/dmccabe/R")

library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(doMC)

ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores / 2))

fqfns <- list.files("/xchip/scarter/dmccabe/germline/tcga_luad/out", "*.hets$", full.names = T)

hfg <- alply(fqfns, 1, function(fqfn) {
  message(fqfn)
  d <- read_tsv(fqfn) %>% transmute(sample = str_extract(basename(fqfn), "^[^.]+"),
                                    r = REF_COUNT, 
                                    a = ALT_COUNT,
                                    n = r + a,
                                    locus = paste0(CONTIG, ":", POSITION))
  return(d)
}, .parallel = wpar)

hets <- do.call(bind_rows, hfg)

save(hets, file = "/xchip/scarter/dmccabe/germline/tcga_luad/hets.rdata")
