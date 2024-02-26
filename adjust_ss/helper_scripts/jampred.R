library(bigsnpr)
.libPaths( c( "/home/kulmsc/R/x86_64-pc-linux-gnu-library/3.6",  .libPaths() ) )
library(R2BGLiMS)

args <- commandArgs(trailingOnly=TRUE)

chr <- args[1]
author <- args[2]
d <- args[3]

if(!file.exists(paste0(d, "/for_jampred.rds"))){
  snp_readBed(paste0(d, "/for_jampred.bed"))
}
obj.bigSNP <- snp_attach(paste0(d, "/for_jampred.rds"))
ped <- big_copy(snp_fastImputeSimple(obj.bigSNP$genotypes, "mean0"), type = "integer")
ped <- ped$bm()
options(bigmemory.allow.dimnames=TRUE)
colnames(ped) <- obj.bigSNP$map$marker.ID
rm(obj.bigSNP)

ss <- read.table(paste0("temp_files/ss.", tolower(author), ".", chr), stringsAsFactors=F, header=T)
ss <- ss[ss$RSID %in% colnames(ped),]
marg_beta <- ss$BETA
names(marg_beta) <- ss$RSID
marg_se <- ss$SE
names(marg_se) <- ss$RSID

meta_stats <- read.table("~/athena/doc_score/raw_ss/meta_stats", sep = ",", stringsAsFactors=F, header = T)
meta_line <- meta_stats[meta_stats[,1] == tools::toTitleCase(author),]

lambdas <- read.table("all_specs/jampred_param_specs", header = T)
try_ind <- unname(which(marg_beta != 0))

for(i in 1:nrow(lambdas)){
  print(str(lambdas))
  print(lambdas[i,1])

  jampred.res.bin <- JAMPred(
   marginal.betas = marg_beta[try_ind],
   n.training = as.numeric(meta_line$sampe_size),
   marginal.logor.ses = marg_se[try_ind], # Only necessary for a binary trait
   p.cases.training = meta_line$cases/meta_line$sampe_size, # Only necessary for a binary trait
   ref.geno = ped[1:nrow(ped),try_ind],
   total.snps.genome.wide = meta_line$snps, # Total SNPs across all chromosomes
   n.mil = 0.2,
   n.cores = 1,
   beta.binom.b.lambda = lambdas[i,1],
   debug = TRUE,
   save.path = paste0("/home/kulmsc/athena/doc_score/adjust_ss/", d),
   seed = 1 # For re-producibility. If not set a random seed is used
  )

  saveRDS(jampred.res.bin, paste0(d, "/jampred.RDS"))
  newss <- ss[ss$RSID %in% names(jampred.res.bin$step2.posterior.mean.snp.weights),]
  newss <- newss[order(newss$RSID)[rank(names(jampred.res.bin$step2.posterior.mean.snp.weights))],]
  newss$BETA <- jampred.res.bin$step2.posterior.mean.snp.weights

  write.table(newss, paste0("../mod_sets/", author, "/", tolower(author), ".", chr,  ".jampred.", i, ".ss"), row.names = F, col.names = T, quote = F, sep = '\t')
}
