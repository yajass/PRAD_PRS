args = commandArgs(trailingOnly=TRUE)
chrom <- args[1]
eth <- args[2]
d <- args[3]
ess <- args[4]
type_ss <- args[5]

ss <- read.table(paste0("../raw_ss/", eth, "/chr_ss/", eth, "_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
new_ss <- ss[,c(3, 1, 8, 4, 5)]
new_ss[,2] <- ess
new_ss[,3] <- abs(qnorm(ss$P)/2) * sign(ss$BETA)
new_ss <- new_ss[!is.na(new_ss[,3]),]

colnames(new_ss) <- c("SNP", "N", "Z", "A1", "A2")
write.table(new_ss, paste0(d, "/", type_ss, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")
