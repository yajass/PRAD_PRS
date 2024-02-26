args <- commandArgs(trailingOnly=TRUE)
eth <- args[1]
d <- args[2]
chr <- args[3]

all_files <- list.files(paste0("../prep_impact/", eth, "_final_results/"))
all_h2 <- rep(0, length(all_files))
for(i in 1:length(all_files)){
  impact_out <- read.table(paste0("../prep_impact/", eth, "_final_results/", all_files[i]), stringsAsFactors=F, header=T)

  check_val <- impact_out$Prop._h2[nrow(impact_out)] - impact_out$Prop._h2_std_error[nrow(impact_out)]
  if(check_val > 0){
    all_h2[i] <- impact_out$Prop._h2[nrow(impact_out)]
  }
}

if(eth == "eur"){
best_ind <- as.numeric(strsplit(all_files[which.max(all_h2)], ".", fixed = T)[[1]][2])
} else {
best_ind <- as.numeric(strsplit(all_files[which.max(all_h2)], ".", fixed = T)[[1]][3])
}
best_ind <- best_ind + 4

system(paste0("zcat ../prep_impact/ready_files/", eth, ".", chr, ".annot.gz | tail -n+2 | cut -f3 > ", d ,"/other_annot_rsids"))
system(paste0("zcat ~/athena/exome_score/funct_anno/impact/IMPACT707_EUR_chr", chr, ".annot.gz | cut -f1-4,", best_ind," | fgrep -w -f ", d, "/other_annot_rsids > ", d, "/curr_annot"))

annot <- read.table(paste0(d, "/curr_annot"), stringsAsFactors=F)
write.table(annot[annot[,5] > quantile(annot[,5], 0.9),3], paste0(d, "/rsid.0.9.txt"), row.names = F, col.names = F, quote = F)
write.table(annot[annot[,5] > quantile(annot[,5], 0.95),3], paste0(d, "/rsid.0.95.txt"), row.names = F, col.names = F, quote = F)
write.table(annot[annot[,5] > quantile(annot[,5], 0.99),3], paste0(d, "/rsid.0.99.txt"), row.names = F, col.names = F, quote = F)

#system("rm temp_files/total_annot")
#for(i in 1:22){
# system(paste0("zcat ../prep_impact/ready_files/", eth, ".", i, ".annot.gz | tail -n+2 | cut -f3 > temp_files/other_annot_rsids"))
# system(paste0("zcat ~/athena/exome_score/funct_anno/impact/IMPACT707_EUR_chr", i, ".annot.gz | cut -f1-4,", best_ind," | fgrep -w -f temp_files/other_annot_rsids > temp_files/curr_annot"))
# system("cat temp_files/curr_annot >> temp_files/total_annot")
#}

