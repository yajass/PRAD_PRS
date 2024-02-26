args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
short_eth <- args[2]
author <- args[3]
chr <- args[4]
low_author <- tolower(author)

print(author)

ss <- read.table(paste0("temp_files/ss.", author, ".", chr), stringsAsFactors=F, header=T)
res <- read.table(paste0(d, "/polyfun_output.agg.txt.gz"), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% res$SNP,]
res <- res[res$SNP %in% ss$RSID,]

res <- res[order(res$SNP)[rank(ss$RSID)],]

ss$BETA <- res$BETA_MEAN

i <- 1
write.table(ss, paste0("~/athena/SPORE/mod_sets/", author, "/", author, ".", chr,  ".polypred.", i ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)

