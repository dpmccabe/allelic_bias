rm(list = ls())

.libPaths("/xchip/scarter/dmccabe/R")

library(dplyr)
library(lme4)
library(boot)

require(doMC)
ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores))

# load("/xchip/scarter/dmccabe/acs_bias/data/sim_hets_small.rdata")
load("/xchip/scarter/dmccabe/germline/tcga_luad/hets.rdata")

locus_counts <- hets %>% group_by(locus) %>% dplyr::summarize(n_samples = n(), N = sum(n))
N_samples <- length(unique(hets$sample))
hets <- semi_join(hets, filter(locus_counts, n_samples >= N_samples * 0.10))

message("Fitting null")
(mod0 <- glmer(cbind(a, r) ~ 0 + (1 | locus), data = hets, family = "binomial"))

message("Fitting full")
(mod <- glmer(cbind(a, r) ~ 0 + (1 | sample) + (1 | locus), data = hets, family = "binomial"))
# save(mod, file = "/xchip/scarter/dmccabe/acs_bias/hier/crossed/sim_mod.rdata")
save(mod, file = "/xchip/scarter/dmccabe/acs_bias/hier/crossed/tcga_luad.rdata")
