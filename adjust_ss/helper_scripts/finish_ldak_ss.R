args <- commandArgs(trailingOnly=TRUE)
lowauthor <- tolower(args[1])

ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", args[2]), stringsAsFactors=F, header=T)
mega <- read.table(paste0(args[3], "/megabayesr.effects.best"), stringsAsFactors=F, header=T)

ss$new <- paste0(ss[,1], ":", ss[,2])
ss <- ss[ss$new %in% mega[,1],]
ss <- ss[order(ss$new)[rank(mega[,1])],]
ss[,7] <- mega[,5]

write.table(ss, paste0("~/athena/doc_score/mod_sets/", args[1], "/", lowauthor, ".", args[2], ".ldak.", args[4], ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")

#FINISH

