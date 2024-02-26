library(CTSLEB)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

temp.dir = args[1]
target_pop = args[2]
helper_pop = args[3]
chrom = args[4]


#switch names
conv <- c("european", "african", "total", "eastasian")
ronv <- c("EUR", "AFR", "TOT", "EAS")
target_pop <- ronv[conv == target_pop]
helper_pop <- ronv[conv == helper_pop]


#read in summary statistics ------------------------------
sum_EUR = read.table(paste0("/home/kulmsc/athena/SPORE/raw_ss/european/chr_ss/european_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
sum_AFR = read.table(paste0("/home/kulmsc/athena/SPORE/raw_ss/african/chr_ss/african_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
sum_TOT = read.table(paste0("/home/kulmsc/athena/SPORE/raw_ss/total/chr_ss/total_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)
sum_EAS = read.table(paste0("/home/kulmsc/athena/SPORE/raw_ss/eastasian/chr_ss/eastasian_", chrom, ".ss.gz"), stringsAsFactors=F, header=T)

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


r2_vec = c(0.3, 0.5, 0.7)
wc_base_vec = c(100)

write.table(sum_other_ref,paste0(temp.dir,"sum_other_ref"),col.names = T,row.names = F,quote=F)
write.table(sum_tar_ref,paste0(temp.dir,"sum_tar_ref"),col.names = T,row.names = F,quote=F)


#Do the 2D clumping -------------------------------------
snp_list = list()
temp = 1
for(r_ind in 1:length(r2_vec)){
  #create the window size given the clumping r2
  wc_vec = round(wc_base_vec/r2_vec[r_ind])
  for(w_ind in 1:length(wc_vec)){
    pthr = 1
    r2thr = r2_vec[r_ind]
    kbpthr = wc_vec[w_ind]
    #for the first group, we perform clumping using EUR popultion as the reference   
    system(paste0("plink ",
    "--bfile /home/kulmsc/athena/refs/1000genomes/", tolower(helper_pop), ".", chrom, " ",
    "--clump ",temp.dir,"sum_other_ref ",
    "--clump-p1 ",pthr," ",
    "--clump-r2 ",r2thr," ",
    "--clump-kb ",kbpthr," ",
    "--out ", temp.dir, helper_pop, "_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    #for the second group, we perform clumping using AFR population as the reference
    system(paste0("plink ",
    "--bfile /home/kulmsc/athena/refs/1000genomes/", tolower(target_pop), ".", chrom, " ",
    "--clump ",temp.dir,"sum_tar_ref ",
    "--clump-p1 ",pthr," ",
    "--clump-r2 ",r2thr," ",
    "--clump-kb ",kbpthr," ",
    "--out ", temp.dir, target_pop, "_ref_CT_rind_",r_ind,"_wcind_",w_ind))
    
    #combine the SNPs from the two clumping groups
    LD_EUR= fread(paste0(temp.dir, helper_pop, "_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD_tar = fread(paste0(temp.dir, target_pop, "_ref_CT_rind_",r_ind,"_wcind_",w_ind,".clumped"))[,3,drop=F]
    LD  = rbind(LD_EUR,LD_tar)
    snp_list[[temp]] = LD
    names(snp_list[[temp]]) = paste0("clump_r2_",r2thr,"_ws_",kbpthr)
    temp = temp + 1
  }
}


saveRDS(snp_list, paste0(temp.dir,"snp_list.", chrom, ".RDS"))

plink_file = PreparePlinkFile(snp_list,sum_com)



score_file = plink_file[[1]]
clump_score_file <- data.frame(score_file)
write.table(score_file,file = paste0(temp.dir,"score_file"),row.names = F,col.names = F,quote=F)
saveRDS(score_file, paste0("ctsleb_scores/adj_ss.", target_pop, ".", helper_pop, ".", chrom, ".RDS"))
write.table(score_file$SNP, paste0(temp.dir, "use_snps"), row.names=F, col.names=F, quote=F)
system(paste0("bgenix -g ~/athena/ukbiobank/imputed/ukbb.", chrom, ".bgen -incl-rsids ", temp.dir, "use_snps > ", temp.dir, "temp.bgen"))

p_value_file = plink_file[[2]]
unique_infor = plink_file[[3]]
saveRDS(unique_infor, paste0("ctsleb_scores/uinfor.", target_pop, ".", helper_pop, ".", chrom, ".RDS"))
pthres <- c(5E-08,5E-05,5E-02,0.5)
options(warn=-1)
q_range = CreateQRange(pthres)
options(warn=0)
q_range$filename <- paste0("p_tar_", 1:nrow(q_range))

write.table(q_range,file = paste0(temp.dir,"q_range_file"),row.names = F,col.names = F,quote=F)



#score the clump data -------------------------------------
p_value_file_temp = p_value_file
for(k1 in 1:length(pthres)){
  idx <- which(unique_infor$P_other<=pthres[k1])
  p_value_file_temp$P[idx] = 0
  write.table(p_value_file_temp,file = paste0(temp.dir,"p_value_file"),col.names = F,row.names = F,quote=F)
  n_col = ncol(score_file)

  res = system(paste0("plink2 ",
     "--q-score-range ",temp.dir,"q_range_file ",temp.dir,"p_value_file ",
     "--score-col-nums 3-",n_col," ",
     "--score ",temp.dir,"score_file  ",
     "--bgen ", temp.dir, "temp.bgen ref-first ",
     "--sample /home/kulmsc/athena/ukbiobank/imputed/ukbb.", chrom, ".sample ",
     "--out ",temp.dir,"prs_p_other_",k1))
}






prs_list = list()
temp = 1
usenames = colnames(score_file)[3:ncol(score_file)]
for(k1 in 1:length(pthres)){
  for(k2 in 1:length(pthres)){
    
    prs_temp = fread(paste0(temp.dir,"prs_p_other_",k1,".p_tar_",k2,".sscore"))
    prs_list[[temp]] = prs_temp[,5:ncol(prs_temp)]*2*nrow(score_file)
   
    colnames(prs_list[[temp]]) = paste0(usenames,"_","p_other_",pthres[k1],"_p_tar_",pthres[k2])
    temp = temp + 1
  }
}
prs_mat = as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
clump_prs_mat <- data.frame(prs_mat)
saveRDS(clump_prs_mat, paste0("ctsleb_scores/score.", target_pop, ".", helper_pop, ".", chrom, ".RDS"))
