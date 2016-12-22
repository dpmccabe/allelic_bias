rm(list = ls())

library(MASS)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(lme4)
library(ggplot2)
library(gridExtra)
library(stringr)
library(randomForest)
library(PSCBS)

require(doMC)
ncores <- detectCores()
if (wpar <- (ncores > 2)) registerDoMC(cores = floor(ncores / 2))

load("/xchip/scarter/dmccabe/acs_bias/hier/crossed/tcga_luad0.rdata")

coef1 <- ranef(mod0)
loci_skews <- data.frame(s_l_hat = unname(coef1$locus), locus = rownames(coef1$locus)) %>% dplyr::as_data_frame() %>%
  mutate(f_hat = plogis(s_l_hat), lambda = (1 - f_hat) / f_hat, locus = as.character(locus)) %>%
  separate(locus, c("contig", "position"), convert = T, remove = F)

baits <- read_tsv("/home/unix/dshih/data/targets/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.tsv")

if (file.exists("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_skew_baits_df.rdata")) {
  load("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_skew_baits_df.rdata")
} else {
  message("Building skew_baits_df")
  skew_baits_df <- ddply(loci_skews, .(locus), function(r) {
    # message(r$locus)
    baits <- filter(baits, contig == r$contig, start <= r$position, end >= r$position)
    n_baits <- nrow(baits)
    
    if (n_baits == 1) {
      return(bind_cols(r %>% mutate(n_baits = n_baits), baits))
    } else {
      return(r %>% mutate(n_baits = n_baits))
    }
  }, .parallel = wpar)
  
  save(skew_baits_df, file = "/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_skew_baits_df.rdata")
}

if (file.exists("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_skew_baits_covs.rdata")) {
  load("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_skew_baits_covs.rdata")
} else {
  message("Building skew_baits_covs")
  
  skew_baits <- dplyr::as_data_frame(skew_baits_df) %>% filter(complete.cases(.)) %>%
    dplyr::select(locus, contig, position, s_l_hat, start, end, sequence)
  
  interval_length <- skew_baits[1,]$end - skew_baits[1,]$start + 1
  
  skew_baits <- skew_baits %>%
    mutate(bait_pos = position - start + 1,
           bait_pos_scaled = (bait_pos - 1) / (interval_length - 1),
           bait_pos_from_mid = abs(bait_pos - (interval_length - 1) / 2),
           bait_pos_from_mid_scaled = (bait_pos_from_mid - min(bait_pos_from_mid)) / (max(bait_pos_from_mid) - min(bait_pos_from_mid)),
           gc = str_count(.$sequence, "G|C") / (end - start))
  
  alleles <- str_split(skew_baits$sequence, "", simplify = T)

  seq_rle <- alply(alleles, 1, function(r) return(rle(r)))
  max_rle <- laply(seq_rle, function(r) return(max(r$lengths)))
  sd_rle <- laply(seq_rle, function(r) return(sd(r$lengths)))
  mean_rle <- laply(seq_rle, function(r) return(mean(r$lengths)))
  n_switches <- laply(seq_rle, function(r) return(length(r$lengths)))
  
  skew_baits <- cbind(skew_baits, max_rle, sd_rle, mean_rle, n_switches)
  
  skew_baits_covs <- ddply(skew_baits, .(bait_pos), function(d) {
    # d <- filter(skew_baits, bait_pos == 40)
    d_bait_pos <- d$bait_pos[1]
    alleles <- str_split(d$sequence, "", simplify = T)
    
    covs_l <- alply(1:10, 1, function(i) {
      left_i <- max(1, (d_bait_pos - i))
      right_i <- min(interval_length, (d_bait_pos + i))
      i_length <- right_i - left_i + 1
      
      neigh <- alleles[, left_i:right_i]
      
      prop_c <- rowMeans(neigh == "C")
      prop_g <- rowMeans(neigh == "G")
      prop_a <- rowMeans(neigh == "A")
      prop_t <- rowMeans(neigh == "T")
      prop_cg <- prop_c + prop_g
      prop_ca <- prop_c + prop_a
      prop_ct <- prop_c + prop_t
      prop_at <- prop_a + prop_t
      prop_ag <- prop_a + prop_g
      prop_gt <- prop_g + prop_t
      
      n_uniq <- apply(neigh, 1, function(r) return(length(unique(r))))
      
      seq_rle <- alply(neigh, 1, function(r) return(rle(r)))
      max_rle <- laply(seq_rle, function(r) return(max(r$lengths)))
      mean_rle <- laply(seq_rle, function(r) return(mean(r$lengths)))
      n_switches <- laply(seq_rle, function(r) return(length(r$lengths)))
      
      covs <- cbind(prop_c, prop_g, prop_a, prop_t, prop_cg, prop_ca, prop_ct, prop_at, prop_ag, prop_gt, n_uniq, max_rle, mean_rle, n_switches)
      colnames(covs) <- paste(colnames(covs), "_", i, sep = "")
      
      return(covs)
    })
    
    covs <- do.call(cbind, covs_l)
    
    ref <- alleles[, d_bait_pos]
  
    if (d_bait_pos > 1) {
      left_ref <- alleles[, d_bait_pos - 1]
    } else {
      left_ref <- "-"
    }
    
    if (d_bait_pos < interval_length) {
      right_ref <- alleles[, d_bait_pos + 1]
    } else {
      right_ref <- "-"
    }
    
    return(cbind(d, left_ref, ref, right_ref, covs))
  }, .parallel = T) %>% as_data_frame()
  
  save(skew_baits_covs, file = "/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_skew_baits_covs.rdata")
}

if (file.exists("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_training_set.rdata")) {
  load("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_training_set.rdata")
  load("/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_test_set.rdata")
} else {
  set.seed(20161216)
  skew_baits_covs <- skew_baits_covs %>% arrange(contig, position)
  test_i <- sample(1:nrow(skew_baits_covs), 0.2 * nrow(skew_baits_covs))
  skew_baits_covs$tset = "training"
  skew_baits_covs[test_i, ]$tset = "test"
  skew_baits_covs <- mutate(skew_baits_covs, tset = as.factor(tset))
  
  training_set <- filter(skew_baits_covs, tset == "training")
  test_set <- filter(skew_baits_covs, tset == "test")
  
  plot_set <- filter(skew_baits_covs, contig == 11)
  plot_cbs_out <- segmentByCBS(plot_set$s_l_hat, chromosome = plot_set$contig, x = plot_set$position, undo = 0, joinSegments = F, verbose = T, alpha = 0.1, min.width = 3)
  plot(plot_cbs_out)
  
  cbs_out <- segmentByCBS(training_set$s_l_hat, chromosome = training_set$contig, x = training_set$position, undo = 0, joinSegments = F, verbose = T)
  cbs_segs <- cbs_out$output %>% as_data_frame() %>% dplyr::select(-sampleName) %>% filter(complete.cases(.))
  
  save(cbs_segs, file = "/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_cbs_segs.rdata")
  
  training_segs <- ddply(training_set, .(locus), function(d) {
    matched_seg <- filter(cbs_segs, chromosome == d$contig, d$position >= start, d$position <= end)
    return(d %>% mutate(seg_mean = matched_seg$mean))
  }, .parallel = T) %>% as_data_frame()
  
  save(training_segs, file = "/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_training_set.rdata")
  
  mean_s_l_hat <- mean(training_set$s_l_hat)
  
  test_segs <- ddply(test_set, .(locus), function(d) {
    matched_seg <- filter(cbs_segs, chromosome == d$contig, d$position >= start, d$position <= end)
    
    if (nrow(matched_seg) == 1) {
      return(d %>% mutate(seg_mean = matched_seg$mean))
    } else {
      return(d %>% mutate(seg_mean = mean_s_l_hat))
    }
  }, .parallel = T) %>% as_data_frame()
  
  save(test_segs, file = "/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_test_set.rdata")
}

colnames(training_segs)
(rf_out <- randomForest(s_l_hat ~ . - tset - contig - position - locus - start - end - sequence - bait_pos_scaled - bait_pos_from_mid - bait_pos_from_mid_scaled, data = training_segs))

save(rf_out, file = "/xchip/scarter/dmccabe/bait_bias/data/tcga_luad_rf_out2.rdata")
stop()
plot(rf_out$y, rf_out$predicted)

rf_imp <- importance(rf_out, scale = T)
t(t(rf_imp[order(-rf_imp),]))

test_pred <- predict(rf_out, test_segs)
plot(test_pred, test_segs$s_l_hat)

mean((test_pred - test_segs$s_l_hat)^2)

ggplot(filter(skew_baits_covs, contig == 1), aes(position, s_l_hat)) + geom_point() +
  theme(axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))) + xlab("Chr11") + ylab(expression(kappa))

ggplot(filter(skew_baits_covs, contig == 1, position >= 1e8, position <= 1e8 * 1.2), aes(position, s_l_hat)) + geom_point() +
  theme(axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))) + xlab("Chr11") + ylab(expression(kappa))
