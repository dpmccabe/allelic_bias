rm(list = ls())

library(dplyr)
library(lme4)
library(readr)
library(ggplot2)
library(gridExtra)
require(doMC)
ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores / 2))

setwd("/xchip/scarter/dmccabe/acs_bias/hier")

compare_lambdas <- function(mod0, mod1, pon_fqfn) {
  print(anova(mod1, mod0))
  
  pon <- read_tsv(pon_fqfn, comment = "#")
  pon <- pon %>% mutate(lambda = ALPHA / BETA, locus = paste(.$CONTIG, .$POSITION, sep = ":"))
  
  coef0 <- ranef(mod0)
  loci_skews0 <- data.frame(s_l_hat = unname(coef0$locus), locus = rownames(coef0$locus)) %>%
    dplyr::as_data_frame()
  coef_compare0 <- left_join(dplyr::select(pon, locus, lambda), loci_skews0) %>%
    filter(complete.cases(.$s_l_hat))
  coef_compare0 <- coef_compare0 %>%
    mutate(f_hat = plogis(s_l_hat), theta = 1 / (1 + lambda), lambda_ranef = (1 - f_hat) / f_hat,
           residual = lambda_ranef - lambda)
  
  coef1 <- ranef(mod1)
  loci_skews1 <- data.frame(s_l_hat = unname(coef1$locus), locus = rownames(coef1$locus)) %>%
    dplyr::as_data_frame()
  coef_compare1 <- left_join(dplyr::select(pon, locus, lambda), loci_skews1) %>%
    filter(complete.cases(.$s_l_hat))
  coef_compare1 <- coef_compare1 %>%
    mutate(f_hat = plogis(s_l_hat), theta = 1 / (1 + lambda), lambda_ranef = (1 - f_hat) / f_hat,
           residual = lambda_ranef - lambda)
  
  ggplot(coef_compare0, aes(lambda, lambda_ranef)) + geom_point() +
    geom_abline(slope = 1, color = "red") + xlim(0, 8) + ylim(0, 8) + xlab(expression(lambda[j]^GATK)) + ylab(expression(lambda[j]^RANEF)) +
    theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) 
  
  lambda_f <- data_frame(theta = plogis(coef1$sample$`(Intercept)`)) %>%
    mutate(lambda = (1 - theta) / theta, skew_type = "sample")
  lambda_l <- data_frame(theta = plogis(coef1$locus$`(Intercept)`)) %>%
    mutate(lambda = (1 - theta) / theta, skew_type = "locus")
  lambdas <- bind_rows(lambda_f, lambda_l)
  
  bias_f <- data_frame(theta = plogis(coef1$sample$`(Intercept)`) - 0.5, skew_type = "sample")
  bias_l <- data_frame(theta = plogis(coef1$locus$`(Intercept)`) - 0.5, skew_type = "locus")
  bias <- bind_rows(bias_f, bias_l)
  
  ggplot(bias, aes(skew_type, theta, color = skew_type, fill = skew_type)) + geom_violin() +
    guides(color = F, fill = F) +
    theme(axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))) + xlab("Bias type") + ylab("")
}

# sim
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/sim_mod0.rdata")
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/sim_mod.rdata")
load("/xchip/scarter/dmccabe/acs_bias/data/sim_hets_small.rdata")

s_l_hat0 <- ranef(mod0)
s_l_hat <- ranef(mod)

s_l_compare0 <- dplyr::data_frame(locus = rownames(s_l_hat0$locus), s_l_hat0 = plogis(s_l_hat0$locus$`(Intercept)`))
s_l_compare <- dplyr::data_frame(locus = rownames(s_l_hat$locus), s_l_hat = plogis(s_l_hat$locus$`(Intercept)`))
s_l <- hets %>% arrange(locus) %>% dplyr::select(locus, s_l) %>% distinct() %>% mutate(s_l = plogis(s_l))

all_s_l <- s_l %>% left_join(s_l_compare0) %>% left_join(s_l_compare)

all_s_l <- all_s_l %>% mutate(lambda_l_hat0 = (1 - s_l_hat0) / s_l_hat0,
                              lambda_l_hat = (1 - s_l_hat) / s_l_hat,
                              lambda_l = (1 - s_l) / s_l)

ggplot(all_s_l, aes(s_l, s_l_hat0)) + geom_point(size = 0.5) + geom_abline(slope = 1, color = "red")
ggplot(all_s_l, aes(s_l, s_l_hat)) + geom_point(size = 0.5) + geom_abline(slope = 1, color = "red")

ggplot(all_s_l, aes(lambda_l, lambda_l_hat0)) + geom_point(size = 0.5) + geom_abline(slope = 1, color = "red") +
  xlab("Actual") + ylab("Estimated") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) 
ggplot(all_s_l, aes(lambda_l, lambda_l_hat)) + geom_point(size = 0.5) + geom_abline(slope = 1, color = "red") +
  xlab("Actual") + ylab("Estimated") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) 

s_f_compare <- dplyr::data_frame(sample = rownames(s_l_hat$sample), s_f_hat = plogis(s_l_hat$sample$`(Intercept)`))
s_f <- hets %>% arrange(sample) %>% dplyr::select(sample, s_f) %>% distinct() %>% mutate(s_f = plogis(s_f))

all_s_f <- s_f %>% left_join(s_f_compare)

all_s_f <- all_s_f %>% mutate(lambda_f_hat = (1 - s_f_hat) / s_f_hat,
                              lambda_f = (1 - s_f) / s_f)

ggplot(all_s_f, aes(s_f, s_f_hat)) + geom_point(size = 0.5) + geom_abline(slope = 1, color = "red")
ggplot(all_s_f, aes(lambda_f, lambda_f_hat)) + geom_point(size = 0.5) + geom_abline(slope = 1, color = "red") +
  xlab("Actual") + ylab("Estimated") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) 

# thca
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/thca_mod0.rdata")
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/thca_mod.rdata")
pon_fqfn <- "/dsde/working/slee/pon/thca-pon-bayes-s30/pulldown/thca-pon-bayes-s30.tsv"
compare_lambdas(mod0, mod, "/dsde/working/slee/pon/thca-pon-bayes-s30/pulldown/thca-pon-bayes-s30.tsv")

# bm ice
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/bm_ice_mod0.rdata")
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/bm_ice_mod.rdata")
compare_lambdas(mod0, mod, "/xchip/scarter/dmccabe/apon/bm_ice/bm_ice.apon.tsv")

# bm agilent
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/bm_ag_mod0.rdata")
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/bm_ag_mod.rdata")
compare_lambdas(mod0, mod, "/xchip/scarter/dmccabe/apon/bm_agilent/bm_agilent.apon.tsv")

# tcga luad
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/tcga_luad0.rdata")
load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/tcga_luad.rdata")
tcga_luad_plots <- compare_lambdas(mod0, mod, "/xchip/scarter/dmccabe/apon/tcga_luad/tcga_luad.apon.tsv")
save(tcga_luad_plots, file = "crossed/plots/tcga_luad_plots.rdata")
