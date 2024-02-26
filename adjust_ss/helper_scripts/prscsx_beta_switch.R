#arguments are: d, i, author, chr 
args <- commandArgs(trailingOnly=TRUE)
d <- args[1]
short_eth <- args[2]
author <- args[3]
chr <- args[4]
i <- args[5]
phi <- args[6]
low_author <- tolower(author)


file_name <- list.files(d, pattern = paste0("x_", short_eth))
print(file_name)
print(paste0("1e-0", nchar(phi)-2))
if(phi == 1){
  file_name <- grep("1e+00", file_name, value = T)
} else {
  file_name <- grep(paste0("1e-0", nchar(phi)-2), file_name, value = T)
}
print(d)
print(file_name)
print(paste0(d, "/", file_name))
beta <- read.table(paste0(d, "/", file_name), stringsAsFactors = F)
ss <- read.table(paste0("temp_files/ss.", low_author, ".", chr), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% beta[,2],]
beta <- beta[order(beta[,2])[rank(ss$RSID)],]
ss$BETA <- beta[,6]

write.table(ss, paste0("~/athena/SPORE/mod_sets/", author, "/", low_author, ".", chr,  ".prscsx.", i ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)

