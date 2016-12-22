rm(list = ls())

.libPaths("/xchip/scarter/dmccabe/R")

library(dplyr)

require(doMC)
ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores / 2))

load("/xchip/scarter/dmccabe/germline/tcga_luad/hets.rdata")

locus_counts <- hets %>% group_by(locus) %>% dplyr::summarize(n_samples = n(), N = sum(n))
N_samples <- length(unique(hets$sample))
hets <- semi_join(hets, filter(locus_counts, n_samples >= N_samples * 0.10))

sample_counts <- (hets %>% group_by(locus) %>% dplyr::summarize(n()))
qplot(sample_counts$`n()`) + theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) +
  xlab(expression(n_samples))

locus_counts <- hets %>% group_by(sample) %>% dplyr::summarize(n())

qplot(locus_counts$`n()`) + theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) +
  xlab(expression(n_loci))

ggplot(sample_n(hets, 10000), aes(r, a)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, color = "red") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("ref") + ylab("alt") +
  xlim(0, 1000) + ylim(0, 1000)

ggplot(sample_n(hets, 10000), aes(r, a)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, color = "red") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("ref") + ylab("alt") +
  xlim(0, 1000) + ylim(0, 1000)

locus_grouped_hets <- ddply(hets, .(locus), function(d) {
  return(data_frame(locus = d$locus[1], sum_a = sum(d$a), sum_r = sum(d$r), n_samples = nrow(d)))
}, .parallel = T) %>% as_data_frame()

sample_grouped_hets <- ddply(hets, .(sample), function(d) {
  return(data_frame(locus = d$locus[1], sum_a = sum(d$a), sum_r = sum(d$r), n_loci = nrow(d)))
}, .parallel = T) %>% as_data_frame()

ggplot(locus_grouped_hets, aes(sum_r, sum_a)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, color = "red") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("ref") + ylab("alt") +
  xlim(0, 4e5) + ylim(0, 4e5)

ggplot(locus_grouped_hets, aes(n_samples, sum_a / (sum_a + sum_r))) + geom_point(alpha = 0.5) +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("n_samples") + ylab(expression(a / (a+r)))

ggplot(sample_grouped_hets, aes(sum_r, sum_a)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, color = "red") +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("ref") + ylab("alt") +
  xlim(0, 1.5e6) + ylim(0, 1.5e6)

ggplot(sample_grouped_hets, aes(n_loci, sum_a / (sum_a + sum_r))) + geom_point(alpha = 0.5) +
  theme(axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5))) + xlab("n_loci") + ylab(expression(a / (a+r)))
