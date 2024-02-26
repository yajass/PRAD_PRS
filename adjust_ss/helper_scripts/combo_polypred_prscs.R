
args = commandArgs(trailingOnly=TRUE)
chrom <- args[1]
ancestry <- args[2]
start_ind <- as.numeric(args[3])
prscs_ind <- args[4]

polyfun <- read.table(paste0("../mod_sets/total/total.", chrom, ".polypred.1.ss"), stringsAsFactors=F, header=T)
prscs <- read.table(paste0("../mod_sets/", ancestry, "/", ancestry, ".", chrom, ".prscs.", prscs_ind, ".ss"), stringsAsFactors=F, header=T)

prscs <- prscs[prscs$RSID %in% polyfun$RSID,]
polyfun <- polyfun[polyfun$RSID %in% prscs$RSID,]
prscs <- prscs[order(prscs$RSID),]
polyfun <- polyfun[order(polyfun$RSID),]

polyfun_beta <- polyfun$BETA
prscs_beta <- prscs$BETA

prscs$BETA <- (1/3)*prscs_beta + (2/3)*polyfun_beta
write.table(prscs, paste0("../mod_sets/", ancestry, "/", ancestry, ".", chrom, ".polyfun.", 1+start_ind, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")

prscs$BETA <- (1/2)*prscs_beta + (2/2)*polyfun_beta
write.table(prscs, paste0("../mod_sets/", ancestry, "/", ancestry, ".", chrom, ".polyfun.", 2+start_ind, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")

prscs$BETA <- (2/3)*prscs_beta + (1/3)*polyfun_beta
write.table(prscs, paste0("../mod_sets/", ancestry, "/", ancestry, ".", chrom, ".polyfun.", 3+start_ind, ".ss"), row.names = F, col.names = T, quote = F, sep = "\t")


