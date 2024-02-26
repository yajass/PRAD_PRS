args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
short_eth <- args[2]
author <- args[3]
chr <- args[4]
low_author <- tolower(author)

print(author)

ss <- read.table(paste0("temp_files/ss.", author, ".", chr), stringsAsFactors=F, header=T)
#ss$MAF <- runif(nrow(ss), 0.002, 0.5)
maf <- read.table(paste0("temp_files/", author, ".", chr, ".frq"), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% maf$SNP,]
maf <- maf[maf$SNP %in% ss$RSID,]

maf <- maf[order(maf$SNP)[rank(ss$RSID)],]

ss$MAF <- maf$MAF

write.table(ss, paste0(d, "/ss.", author, ".", chr, ".maf"), row.names = F, col.names = T, quote = F, sep = "\t")
