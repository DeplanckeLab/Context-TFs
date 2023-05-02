# # Functions for RF cluster mechanism analysis incroporating GM12878 data

# the following table, takes a peak file, a variant file, and the peak_IDs that want to be extracted as input and merges peak and variant file
# in addition further information on the SNP is given, like distance to peak center
## USED
get_peak_var_tab=function(peaks,vars,PIDs){
  gvars =do.call(rbind,lapply(PIDs,function(x){
    vars[vars$Peak == x,]}))
  gpeaks = do.call(rbind,lapply(PIDs,function(x){
    peaks[peaks$Peak == x,]}))
  vcmf = cbind(gvars,gpeaks[match(gvars$Peak,gpeaks$Peak),])
  vcmf$Peak_width = vcmf$Pos_Right - vcmf$Pos_Left
  vcmf$Peak_center = vcmf$Pos_Left + round( (vcmf$Peak_width/2),0)
  vcmf$Dist_SNP_to_peak_center = abs(vcmf$Pos-vcmf$Peak_center)
  vcmf$SNP_centrality = vcmf$Dist_SNP_to_peak_center / vcmf$Peak_width
  vcmf$Chr = paste0("chr",vcmf$Chr)
  return(vcmf)
}
## USED
# making a GRanges object from df, given specif column-ids
make_GR_from_df = function(df,cids){
  makeGRangesFromDataFrame(data.frame(chr = df[,cids[1]],start = df[,cids[2]],end=df[,cids[3]]))
}

## adding the info of dependent enhancer to a lead peak var tab, this requires a list of DEP Peak_ids for each LEad_Peak_id
## USED
add_dep_info = function(LEAD_rs_tab,lead_dep_list,vars,peaks){
  do.call(rbind,lapply(1:length(LEAD_rs_tab[,1]),function(x){
    deps=lead_dep_list[[match(LEAD_rs_tab$Peak[x],names(lead_dep_list))]]
    ddf=do.call(rbind,lapply(deps,function(j){
      get_peak_var_tab(peaks,vars,j)
    }))
    colnames(ddf) = paste0("dep_",colnames(ddf))
    jdf=cbind(LEAD_rs_tab[rep(x,length(ddf[,1])),],ddf)
    return(jdf)}))
}







# ## Adding in ASA and ASB info for lead and dependen SNPs

# +
##USED

get_leadSNPs_in_GM_info = function(rs_tab,ASB_mat){
  keep_ids =match(rs_tab$RsID,ASB_mat$ID)
  keep=rs_tab
  keep$GM_RsID = ASB_mat$ID[keep_ids]
  keep$GM_REF = ASB_mat$REF[keep_ids]
  keep$GM_ALT = ASB_mat$ALT[keep_ids]
  keep$GM_Pos = ASB_mat$POS[keep_ids]
  keep$GM_phasing = ASB_mat$GENOTYPE[keep_ids]
  keep$Dist_GM_SNP_to_peak_center = abs(keep$GM_Pos - keep$Peak_center)
  #keepf = keep[!is.na(keep$GM_RsID),]
  return(keep)
}
## USED
# ASAB
get_dom_ASA_in_GM = function(snp_tab,encode_map,encode_cids,cnames=c("c_ref","c_alt","c_act","c_pas")){
  cactid = c(1,2)[as.integer(snp_tab$Beta>0 )+1]
  tab = snp_tab
  tab$ref = encode_map[match(snp_tab$RsID,encode_map$ID),encode_cids[1]]
  tab$alt = encode_map[match(snp_tab$RsID,encode_map$ID),encode_cids[2]]
  ad=encode_map[match(snp_tab$RsID,encode_map$ID),encode_cids]
  tab$act = sapply(1:length(tab[,1]),function(j){
    ad[j,cactid[j]]})
  tab$pas = sapply(1:length(tab[,1]),function(j){
    ad[j,setdiff(c(1,2),cactid[j])]})
  colnames(tab)[(length(snp_tab[1,])+1):length(tab[1,])]=cnames
  return(tab)  
}



# +
# adding in GM SNPs in dependent enhancer


get_SNPs_in_deppeak_info = function(rs_tab,ASB_mat){
  GR_GM_SNPS = makeGRangesFromDataFrame(data.frame(chr = paste0("chr",ASB_mat$X.CHROM),start = ASB_mat$POS,end=ASB_mat$POS))
  GR_rs = makeGRangesFromDataFrame(data.frame(chr=rs_tab$dep_Chr,start = rs_tab$dep_Peak_center - 350,end =  rs_tab$dep_Peak_center + 349))
  ovrl = as.data.frame(findOverlaps(GR_rs,GR_GM_SNPs,maxgap = 0))
  max_snps = max(sapply(1:length(rs_tab[,1]),function(i){
    length(which(ovrl[,1]==i))
  }))
  dep_GM_info = lapply(1:length(rs_tab[,1]),function(x){
    if (length(which(ovrl[,1]==x))==0){
      df=data.frame(lead_peak_id = c(rs_tab$Peak[x],rep(NA,(max_snps-1))),dep_peakid = c(rs_tab$dep_Peak[x],rep(NA,(max_snps-1))),
                        depGMSNP_RsID = NA, depGMSNP_REF = NA,depGMSNP_ALT=NA,depGMSNP_Pos =NA,depGMSNP_phasing = NA,depGMSNP_dist_to_dep_peak_center = NA)
      return(df)}
    if (length(which(ovrl[,1]==x))>0){
      ids = which(ovrl[,1]==x)
      d = length(ids)
      df = data.frame(lead_peak_id = c(rs_tab$Peak[ovrl[ids,1]],rep(NA,(max_snps-d))),dep_peakid = c(rs_tab$dep_Peak[ovrl[ids,1]],rep(NA,(max_snps-d))),
                      depGMSNP_RsID = c(ASB_mat$ID[ovrl[ids,2]],rep(NA,(max_snps-d))), depGMSNP_REF = c(ASB_mat$REF[ovrl[ids,2]],rep(NA,(max_snps-d))),
                      depGMSNP_ALT=c(ASB_mat$ALT[ovrl[ids,2]],rep(NA,(max_snps-d))),depGMSNP_Pos =c(ASB_mat$POS[ovrl[ids,2]],rep(NA,(max_snps-d))),
                      depGMSNP_phasing = c(ASB_mat$GENOTYPE[ovrl[ids,2]],rep(NA,(max_snps-d))),
                      depGMSNP_dist_to_dep_peak_center = c(abs(ASB_mat$POS[ovrl[ids,2]] - rs_tab$dep_Peak_center[ovrl[ids,1]]),rep(NA,(max_snps-d))))
      return(df)}
  })
  names(dep_GM_info) = paste0("depPEAD_eq_",rs_tab$dep_Peak)
  return(dep_GM_info)
}
# -

# get GM-SNP info in dep peaks together with ATAC coverage

## USED
get_SNPs_in_deppeak_info_plus = function(rs_tab,ASB_mat){
  GR_GM_SNPS = makeGRangesFromDataFrame(data.frame(chr = paste0("chr",ASB_mat$X.CHROM),start = ASB_mat$POS,end=ASB_mat$POS))
  GR_rs = makeGRangesFromDataFrame(data.frame(chr=rs_tab$dep_Chr,start = rs_tab$dep_Peak_center - 350,end =  rs_tab$dep_Peak_center + 349))
  ovrl = as.data.frame(findOverlaps(GR_rs,GR_GM_SNPs,maxgap = 0))
  max_snps = max(sapply(1:length(rs_tab[,1]),function(i){
    length(which(ovrl[,1]==i))
  }))
  dep_GM_info = lapply(1:length(rs_tab[,1]),function(x){
    if (length(which(ovrl[,1]==x))==0){
      df=data.frame(lead_peak_id = c(rs_tab$Peak[x],rep(NA,(max_snps-1))),dep_peakid = c(rs_tab$dep_Peak[x],rep(NA,(max_snps-1))),
                    depGMSNP_RsID = NA, depGMSNP_REF = NA,depGMSNP_ALT=NA,depGMSNP_Pos =NA,depGMSNP_phasing = NA,depGMSNP_dist_to_dep_peak_center = NA,
                    depGMSNP_refcount= NA,
                    depGMSNP_altcount= NA)
      return(df)}
    if (length(which(ovrl[,1]==x))>0){
      ids = which(ovrl[,1]==x)
      d = length(ids)
      df = data.frame(lead_peak_id = c(rs_tab$Peak[ovrl[ids,1]],rep(NA,(max_snps-d))),dep_peakid = c(rs_tab$dep_Peak[ovrl[ids,1]],rep(NA,(max_snps-d))),
                      depGMSNP_RsID = c(ASB_mat$ID[ovrl[ids,2]],rep(NA,(max_snps-d))), depGMSNP_REF = c(ASB_mat$REF[ovrl[ids,2]],rep(NA,(max_snps-d))),
                      depGMSNP_ALT=c(ASB_mat$ALT[ovrl[ids,2]],rep(NA,(max_snps-d))),depGMSNP_Pos =c(ASB_mat$POS[ovrl[ids,2]],rep(NA,(max_snps-d))),
                      depGMSNP_phasing = c(ASB_mat$GENOTYPE[ovrl[ids,2]],rep(NA,(max_snps-d))),
                      depGMSNP_dist_to_dep_peak_center = c(abs(ASB_mat$POS[ovrl[ids,2]] - rs_tab$dep_Peak_center[ovrl[ids,1]]),rep(NA,(max_snps-d))),
                      depGMSNP_refcount= c(ASB_mat[ovrl[ids,2],7],rep(NA,(max_snps-d))),
                      depGMSNP_altcount= c(ASB_mat[ovrl[ids,2],8],rep(NA,(max_snps-d))))
      return(df)}
  })
  names(dep_GM_info) = paste0("depPEAD_eq_",rs_tab$dep_Peak)
  return(dep_GM_info)
}

# # add ASA for depGM

# USED
get_dep_ASA_in_GM = function(snp_tab,encode_map,encode_cids,cnames=c("c_ref","c_alt","c_act","c_pas")){
  altact = as.integer(snp_tab$Beta>0)
  refact = as.integer(snp_tab$Beta<0)
  geno_gm_same = snp_tab$depGMSNP_phasing == snp_tab$GM_phasing
  cactid = c(1,2)[as.integer(altact==T &geno_gm_same==T | refact==T & geno_gm_same==F)+1]
  tab = snp_tab
  tab$ref = encode_map[match(snp_tab$depGMSNP_RsID,encode_map$ID),encode_cids[1]]
  tab$alt = encode_map[match(snp_tab$depGMSNP_RsID,encode_map$ID),encode_cids[2]]
  ad=encode_map[match(snp_tab$depGMSNP_RsID,encode_map$ID),encode_cids]
  tab$act = sapply(1:length(tab[,1]),function(j){
    if (is.na(cactid[j])){
      return(NA)}
    if (!is.na(cactid[j])){
      return(ad[j,cactid[j]])}
    })
  tab$pas = sapply(1:length(tab[,1]),function(j){
    if (is.na(cactid[j])){
      return(NA)}
    if (!is.na(cactid[j])){
      return(ad[j,setdiff(c(1,2),cactid[j])])}
    })
  colnames(tab)[(length(snp_tab[1,])+1):length(tab[1,])]=cnames
  return(tab)  
}

