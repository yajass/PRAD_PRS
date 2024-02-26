chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


if [ $author == "total" ];then

  ess_eur=`cat temp_files/ss.european.${chr} | head -2 | tail -n +2 | cut -f9`
  ess_afr=`cat temp_files/ss.african.${chr} | head -2 | tail -n +2 | cut -f9`
  ess_eas=`cat temp_files/ss.eastasian.${chr} | head -2 | tail -n +2 | cut -f9`
  ess_tot=`cat temp_files/ss.total.${chr} | head -2 | tail -n +2 | cut -f9`


  cp helper_scripts/prscs_header ${d}/ss.eur
  cat temp_files/ss.european.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.eur

  cp helper_scripts/prscs_header ${d}/ss.afr
  cat temp_files/ss.african.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.afr

  cp helper_scripts/prscs_header ${d}/ss.eas
  cat temp_files/ss.eastasian.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.eas

  cp helper_scripts/prscs_header ${d}/ss.tot
  cat temp_files/ss.total.${chr} | tail -n+2 | cut -f3,4,5,7,8 >> ${d}/ss.tot



  i=1 #cannot create prscsx for total because there is not any reference created
      #or I could just use european ?
  cat all_specs/prscs_param_specs | while read phi;do

    let i=i+1

    python ~/Programs/PRScsx/PRScsx.py --ref_dir=/home/kulmsc/athena/refs --bim_prefix=geno_files/${low_author}.${chr} --sst_file=${d}/ss.eur,${d}/ss.afr,${d}/ss.eas --phi=${phi} --n_gwas=${ess_eur},${ess_afr},${ess_eas} --pop=eur,afr,eas --chrom=$chr --out_dir=${d} --out_name=x_eur,x_afr,x_eas

    Rscript helper_scripts/prscsx_beta_switch.R $d afr african $chr $i $phi
    Rscript helper_scripts/prscsx_beta_switch.R $d eas eastasian $chr $i $phi
    Rscript helper_scripts/prscsx_beta_switch.R $d hsp hispanic $chr $i $phi
    Rscript helper_scripts/prscsx_beta_switch.R $d eur european $chr $i $phi

    let i=i+1


done


fi
