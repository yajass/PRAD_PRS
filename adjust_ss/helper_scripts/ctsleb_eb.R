library(CTSLEB)
library(data.table)
library(dplyr)


temp.dir = "ctsleb/run/temp.dir/"
data.dir = "ctsleb/data/"

target_pop="AFR"
helper_pop="EUR"
chrom = 22



sum_EUR = read.table(paste0("../raw_ss/european/chr_ss/european_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
sum_AFR = read.table(paste0("../raw_ss/african/chr_ss/african_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
sum_TOT = read.table(paste0("../raw_ss/total/chr_ss/total_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
sum_EAS = read.table(paste0("../raw_ss/eastasian/chr_ss/eastasian_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)

colnames(sum_EUR)[3] <- "SNP"
colnames(sum_AFR)[3] <- "SNP"
colnames(sum_TOT)[3] <- "SNP"
colnames(sum_EAS)[3] <- "SNP"

sum_other_list = list(sum_EUR,sum_AFR,sum_TOT,sum_EAS)
other_ans_names = c("EUR","AFR","TOT","EAS")

#this is used with EB
sum_com_many <- AlignSumMulti(sum_tar = sum_AFR,
                    sum_other_list = sum_other_list,
                    other_ans_names = other_ans_names)

sum_com <- AlignSum(sum_tar = sum_other_list[[which(other_ans_names == target_pop)]],
                    sum_other = sum_other_list[[which(other_ans_names == helper_pop)]])

sum_com_split <- SplitSum(sum_com)
sum_other_ref = sum_com_split[[1]]
sum_tar_ref = sum_com_split[[2]]


#start here by combining together the snp_lists

snp_list <- list()
for(i in 1:length(big_list[[1]])){
  snp_list[[i]] <- do.call("rbind", lapply(big_list, function(x) x[[i]]))
}


prs_mat <- readRDS(paste0("ctsleb_scores/score.", target_pop, ".", helper_pop, ".", chrom, ".RDS"))
unique_infor <- paste0("ctsleb_scores/uinfor.", target_pop, ".", helper_pop, ".", chrom, ".RDS")
score_file <- paste0("ctsleb_scores/adj_ss.", target_pop, ".", helper_pop, ".", chrom, ".RDS")

prs_tun = prs_mat[1:10000,]
max_ind <- 8



snp_set_ind = colnames(prs_tun)[max_ind]
SNP_set = GetSNPSet(snp_set_ind, score_file, unique_infor)
#unique_infor_post = EBpost(unique_infor,SNP_set)
unique_infor_post = EBpostMulti(unique_infor,SNP_set,
                        sum_com_many,other_ans_names)


#post_beta_mat = cbind(unique_infor_post$BETA_EB_target,unique_infor_post$BETA_EB_other)
#colnames(post_beta_mat) = c("EB_target","EB_eur")
eb_post_col_names = c("BETA_EB_target",paste0("BETA_EB_",other_ans_names[1]))
post_beta_mat = unique_infor_post %>% 
  select(all_of(eb_post_col_names))


#plink_file_eb = PreparePlinkFileEB(snp_list, unique_infor_post, post_beta_mat)
plink_file_eb = PreparePlinkFileEB(snp_list,
                            unique_infor_post,
                            post_beta_mat)


score_file = plink_file_eb[[1]]
eb_score_file <- data.frame(score_file)
saveRDS(paste0("ctsleb_scores/eb_ss.", target_pop, ".", helper_pop, ".", chrom, ".RDS"))
write.table(score_file,file = paste0(temp.dir,"score_file_eb"),row.names = F,col.names = F,quote=F)

p_value_file = plink_file_eb[[2]]


p_value_file_temp = p_value_file
for(k1 in 1:length(pthres)){
  idx <- which(unique_infor$P_other<=pthres[k1])
  p_value_file_temp$P[idx] = 0
  write.table(p_value_file_temp,file = paste0(temp.dir,"p_value_file"),col.names = F,row.names = F,quote=F)
  n_col = ncol(score_file)

  res = system(paste0("plink2 ",
    "--q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file ",
    "--score-col-nums 3-",n_col," ",
    "--score ",temp.dir,"score_file_eb  ",
    "--bgen /home/kulmsc/athena/ukbiobank/imputed/ukbb.", chrom, ".bgen ref-first ",
    "--sample /home/kulmsc/athena/ukbiobank/imputed/ukbb.", chrom, ".sample ",
    "--out ",temp.dir,"eb_prs_p_other_",k1))
}



prs_list = list()
temp = 1
names = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    prs_temp = fread(paste0(temp.dir,"eb_prs_p_other_",k1,".p_tar_",k2,".sscore"))
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]*2*nrow(score_file)
   
    colnames(prs_list[[temp]]) = paste0(names,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp = temp + 1
  }
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
eb_prs_mat <- data.frame(prs_mat)
saveRDS(eb_prs_mat, paste0("ctsleb_scores/eb_mat.", target_pop, ".", helper_pop, ".", chrom, ".RDS"))

##############################

#now we have to save the clump_prs_mat and eb_prs_mat
#similarly with clump_score_file and eb_score_file

