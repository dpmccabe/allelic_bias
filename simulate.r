rm(list = ls())

.libPaths("/xchip/scarter/dmccabe/R")

library(dplyr)
library(truncnorm)
library(stringr)

setwd("/xchip/scarter/dmccabe/acs_bias/hier")

set.seed(1)

N <- 200
M <- 1000

samples <- paste0("s", str_pad(1:N, 3, pad = "0"))
loci <- paste0("l", str_pad(1:M, 5, pad = "0"))

s_f <- rnorm(N, 0, 0.2)
s_l <- rnorm(M, 0, 0.4)
depths <- round(runif(N * M, 30, 300))

hets <- dplyr::as_data_frame(expand.grid(sample = samples, locus = loci, stringsAsFactors = F))
hets <- hets %>% mutate(n = depths, a = NA, r = NA, s_f = rep(s_f, times = M),
                        s_l = rep(s_l, each = N), f = plogis(s_f + s_l)) %>%
  mutate(a = rbinom(N * M, size = n, prob = f), r = n - a) %>% arrange(sample, locus)

hets <- sample_frac(hets, 0.75)

save(hets, file = "/xchip/scarter/dmccabe/acs_bias/data/sim_hets_small.rdata")
