
df <- readRDS("~/athena/pgs_working_group/ancestry_match/data/matched-pop.rds")
eid <- readRDS("~/athena/pgs_working_group/ancestry_match/data/ukbb_eid.RDS")
df$eid <- eid

qc <- read.table("~/athena/ukbiobank/custom_qc/fam_files/qc_fam", stringsAsFactors=F)

df <- df[!is.na(df$pop),]
df <- df[!is.na(df$continent),]
df <- df[df$eid %in% qc[,1],]

set.seed(34)
df <- df[sample(1:nrow(df)),]

euro_eid <- df$eid[df$continent == "Europe"][1:round(sum(df$continent == "Europe")*0.2)]
africa_eid <- df$eid[df$continent == "Africa"][1:round(sum(df$continent == "Africa")*0.2)]
asia_eid <- df$eid[df$continent == "Asia"][1:round(sum(df$continent == "Asia")*0.2)]
total_eid <- df$eid[round(nrow(df)*0.2)]

euro_eid <- euro_eid[1:10000]
total_eid <- total_eid[1:10000]

write.table(euro_eid, "helper_files/european_eid", row.names = F, col.names = F, quote = F)
write.table(africa_eid, "helper_files/african_eid", row.names = F, col.names = F, quote = F)
write.table(asia_eid, "helper_files/eastasian_eid", row.names = F, col.names = F, quote = F)
write.table(total_eid, "helper_files/total_eid", row.names = F, col.names = F, quote = F)
