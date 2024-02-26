library(XPASS)
library(data.table)
#library(RhpcBLASctl)
#blas_set_num_threads(30)
source("local_xpass.R")
print("TOOOOOOOOOOOOOOOOOP")

args = commandArgs(trailingOnly=TRUE)

chrom <- args[1]
target_eth <- args[2]
aux_eth <- args[3]
d <- args[4]
len_target <- args[5]
len_aux <- args[6]
big_index <- as.numeric(args[7])

author <- target_eth

change_names <- data.frame("new" = c("eur", "tot", "hsp", "afr", "eas"),
                           "old" = c("european", "total", "hispanic", "african", "eastasian"), stringsAsFactors = F)
short_target_eth <- change_names$new[change_names$old == target_eth]
short_aux_eth <- change_names$new[change_names$old == aux_eth]


fit_bbj <-localXPASS(file_z1 = paste0(d, "/target.ss"), file_z2 = paste0(d, "/aux.ss"),
                file_ref1 = paste0("../../refs/1000genomes/", short_target_eth, ".", chrom),
                file_ref2 = paste0("../../refs/1000genomes/", short_aux_eth, ".", chrom),
                file_cov1 = paste0("../get_pca/", short_target_eth, ".proper_pca"),
                file_cov2 = paste0("../get_pca/", short_aux_eth, ".proper_pca"),
                compPosMean = T,
                file_out = paste0(d, "/BMI_bbj_ukb_ref_TGP"))
saveRDS(fit_bbj, paste0(d, "/out.RDS"))
ss <- read.table(paste0("../raw_ss/", target_eth, "/chr_ss/", target_eth, "_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
ss <- ss[ss$RSID %in% fit_bbj$mu$SNP,]
ss$BETA <- fit_bbj$mu$mu_XPASS1

write.table(ss, paste0("../mod_sets/", author, "/", author, ".", chrom, ".xpass.", big_index, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")
print("first write")
print(dim(ss))
print(paste("big_index", big_index))


if(as.numeric(len_target) == 0){
  print("clump 1")
  use_snps <- read.table(paste0(d, "/out_target.clumped"), stringsAsFactors=F, header=T)

  fit_bbj <-localXPASS(file_z1 = paste0(d, "/target.ss"), file_z2 = paste0(d, "/aux.ss"),
                file_ref1 = paste0("../../refs/1000genomes/", short_target_eth, ".", chrom),
                file_ref2 = paste0("../../refs/1000genomes/", short_aux_eth, ".", chrom),
                file_cov1 = paste0("../get_pca/", short_target_eth, ".proper_pca"),
                file_cov2 = paste0("../get_pca/", short_aux_eth, ".proper_pca"),
                snps_fe1 = use_snps$SNP,
                compPosMean = T,
                file_out = paste0(d, "/BMI_bbj_ukb_ref_TGP"))
  saveRDS(fit_bbj, paste0(d, "/out.RDS"))
  ss <- read.table(paste0("../raw_ss/", target_eth, "/chr_ss/", target_eth, "_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
  ss <- ss[ss$RSID %in% fit_bbj$mu$SNP,]
  ss$BETA <- fit_bbj$mu$mu_XPASS1

  write.table(ss, paste0("../mod_sets/", author, "/", author, ".", chrom, ".xpass.", big_index+1, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")
  print("second write")
  print(dim(ss))
  print(paste("big_index", big_index+1))
}



if(as.numeric(len_target) == 0 & as.numeric(len_aux) == 0){
  print("clump 2")
  use_snps <- read.table(paste0(d, "/out_target.clumped"), stringsAsFactors=F, header=T)
  use_aux_snps <- read.table(paste0(d, "/out_aux.clumped"), stringsAsFactors=F, header=T)

  fit_bbj <-localXPASS(file_z1 = paste0(d, "/target.ss"), file_z2 = paste0(d, "/aux.ss"),
                file_ref1 = paste0("../../refs/1000genomes/", short_target_eth, ".", chrom),
                file_ref2 = paste0("../../refs/1000genomes/", short_aux_eth, ".", chrom),
                file_cov1 = paste0("../get_pca/", short_target_eth, ".proper_pca"),
                file_cov2 = paste0("../get_pca/", short_aux_eth, ".proper_pca"),
                snps_fe1 = use_snps$SNP, snps_fe2 = use_aux_snps$SNP,
                compPosMean = T,
                file_out = paste0(d, "/BMI_bbj_ukb_ref_TGP"))
  saveRDS(fit_bbj, paste0(d, "/out.RDS"))
  ss <- read.table(paste0("../raw_ss/", target_eth, "/chr_ss/", target_eth, "_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
  ss <- ss[ss$RSID %in% fit_bbj$mu$SNP,]
  ss$BETA <- fit_bbj$mu$mu_XPASS1

  write.table(ss, paste0("../mod_sets/", author, "/", author, ".", chrom, ".xpass.", big_index+2, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")
  print("third write")
  print(dim(ss))
  print(paste("big_index", big_index+2))
}


