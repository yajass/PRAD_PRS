
args <- commandArgs(trailingOnly=TRUE)

chr <- args[1]
author <- args[2]
d <- args[3]
i <- args[4]
lowauthor <- tolower(author)

dbslmm <- read.table(paste0(d,"/ss.", tolower(author), ".dbslmm.txt"), stringsAsFactors=F)
ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", chr), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% dbslmm[,1],]
dbslmm <- dbslmm[order(dbslmm[,1])[rank(ss$RSID)],]
ss$BETA <- dbslmm[,3]

write.table(ss, paste0("../mod_sets/", author, "/", lowauthor, ".", chr, ".dbslmm.", i, ".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
