
args <- commandArgs(trailingOnly=TRUE)

chr <- args[1]
author <- args[2]
d <- args[3]
i <- args[4]
lowauthor <- tolower(author)

sbay <- read.table(paste0(d,"/sbay_out.snpRes"), stringsAsFactors=F, header = T)
ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", chr), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% sbay[,2],]
sbay <- sbay[order(sbay[,2])[rank(ss$RSID)],]
ss$BETA <- sbay$A1Effect

write.table(ss, paste0("../mod_sets/", author, "/", lowauthor, ".", chr, ".sbayesr.", i, ".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
