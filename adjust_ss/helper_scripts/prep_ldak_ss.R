args <- commandArgs(trailingOnly=TRUE)
lowauthor <- tolower(args[1])

ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", args[2]), stringsAsFactors=F, header=T)
ss1 <- ss[,c(3,4,5,7,8)]
ss1[,4] <- sign(ss1[,4])
ss <- ss[ss1[,4] != 0,]
ss1 <- ss1[ss1[,4] != 0,]

colnames(ss1) <- c("Predictor", "A1", "A2", "Direction", "P")


write.table(ss1, paste0("temp_files/ss.", lowauthor, ".", args[2],".ldak"), row.names = F, col.names = T, sep = '\t', quote = F)

ss1[,1] <- paste0(ss$CHR, ":", ss$BP)
write.table(ss1, paste0("temp_files/ss.", lowauthor, ".", args[2],".pos_ldak"), row.names = F, col.names = T, sep = '\t', quote = F)
