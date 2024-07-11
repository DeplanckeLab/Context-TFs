## This file contains the following functions: 

# 1) initiator score
# 2) context score
# 3) IP/ATAC-seq coverage
# 4) association between TF score and epigenomic features (K27 and BRD4) (use randomly sampled enhancers of LCL data set)
# 5) GC and p_caQTL matching in lead enhancers
# 6) distance matching between enhancer pairs to CM enhancers
# 7) Z-score and SNP disruption frequency across TFs for lead enhancer classes
# 8) Z-scoreacross TFs for neighbor/dep enhancer classes




# 1) initiator score
generate_caqtl_hyper = function(Enh_tab,peak_id_sc,dom_act_top_ae,dom_pas_top_ae,TFs_hoco){
    # Enh_tab = master peak table
    # peak_id_sc = peak IDs in order as they occurr in the  matrix containing the peak top TFBS scores
    # dom_act_top_ae = TFBS scores in the more accessible genotype (dimension = n - Enh x 401 TFs)
    # dom_pas_top_ae = TFBS scores in the less accessible genotype (dimension = n - Enh x 401 TFs)
    # TFs_hoco = TF names from hocomoco (see colnames of peak scoring matrix)
    
    dom_act_top = dom_act_top_ae[match(Enh_tab$Peak,peak_id_sc),]
    dom_pas_top = dom_pas_top_ae[match(Enh_tab$Peak,peak_id_sc),]
    # hypergeometric test & fisher.test
    N_TFBSdisr_per_TF_dhyper = do.call(rbind,lapply(1:401,function(i){
        dc=do.call(rbind,tapply((dom_act_top[,i]-dom_pas_top[,i])>0,Enh_tab$causal_label,table)) # how many SNPs in each enhancer class create a new best site in more versus less acc genotype?
        sig=phyper(dc[1,2],m=sum(dc[1,]),n=sum(dc[2,]),k=sum(dc[,2]),lower.tail=FALSE)
        OR = fisher.test(dc[2:1,])$estimate # initiator score
        SNP_hypdf = data.frame(TF=TFs_hoco[i],phyper_SNP = sig,OR_SNP=OR)
        return(SNP_hypdf)
    }))
    return(N_TFBSdisr_per_TF_dhyper)
}


# 2) computation of context score
get_effsize_ceofs_perTF = function(Enh_tab,dom_peak_id_sc,scoreslacc,scoresmacc,TF,label,atac,levels){
    # Enh_tab = master peak table
    # dom_peak_id_sc = peak IDs in order as they occurr in the  matrix containing the peak top TFBS scores
    # scoreslacc,scoresmacc = scores in less and more acc genotype
    # TF =  TF name for which context score is to be generated
    # label = 'causal_label' column of Master peak table used as input
    # atac = column in master table giving ATAC-seq bigwig coverage for each peak ('Lead_bwcum_ATAC_GM12878dedup')
    # levels = levels of label (to assure correct order, i.e. 'no-effect' & 'caQTL')
    # 
    scoreless=scoreslacc[match(Enh_tab$Peak,dom_peak_id_sc),which(TFs_hoco==TF)]
    scoremore=scoresmacc[match(Enh_tab$Peak,dom_peak_id_sc),which(TFs_hoco==TF)]
    mod_df = data.frame(motif_score = scoreless,ATAC = log10(atac),E_type=label,motif_scoremore = scoremore)
    mod_df$E_type = factor(mod_df$E_type,levels=levels)
    fdiff = !((mod_df$motif_scoremore-mod_df$motif_score)>0) ## this removes all enhancers where the SNP creates a new best site, meaning the motif strength might not come from the seqeunce context, but the SNP
    fna = !is.na(mod_df$motif_score) & is.finite(mod_df$motif_score)
    fatac = !is.na(mod_df$ATAC) & is.finite(mod_df$ATAC)
    linmod=lm(motif_score ~ ATAC + E_type ,data=mod_df[fdiff & fna & fatac ,])
    df = data.frame(qtlATAC_coef = summary(linmod)$coef[2,1],qtlATAC_tstat = summary(linmod)$coef[2,3],qtlATAC_pval = summary(linmod)$coef[2,4],
                    qtlEtype_coef = summary(linmod)$coef[3,1],qtlEtype_tstat = summary(linmod)$coef[3,3],qtlEtype_pval = summary(linmod)$coef[3,4],
                    tf=TF)
    return(df)
}


# 3) get IP or ATAC-seq coverage
get_peakcum_bw_score=function(gr,bw){
    # gr = GRanges object containing the peak coordinates (e.g. +/- 350bp or +/- 650bp around peak center
    # bw = bigwig file with ATAC or IP coverage
    ovrl = as.data.frame(findOverlaps(gr,bw))
    s=sapply(1:length(gr),function(i){
        sum(bw$score[ovrl[which(ovrl[,1]==i),2]],na.rm=T)
    }) # cumulative coverage
return(s)}


# 4) TF score to epigenome association while controlling for ATAC
get_Feat_betas_perTF_comb = function(scores,feat1,feat2,feat3,scrtop,tf,co,lab){
    # scores = TF motif scores across enhancers
    # feat1 = ATAC coverage of enhancers
    # feat2 = K27 IP coverage
    # feat3 = BRD4 IP covera
    # scrtop = scores of scrambled peak seqeunces to control for base composition bias
    # co = column id of linear model coefficient summary to extrac ( 1 = coef, 3= t-statistic, 4= p.val,..)
    # label = label to name resulting data frame (1 = coef, 3= t-statistic, 4= p.val)
    
    mod_df = data.frame(motif_score = scores,FEAT1 = log10(feat1),FEAT2 = log10(feat2),FEAT3 = log10(feat3),scr_topsc = scrtop)
    fna = !is.na(mod_df$motif_score) & is.finite(mod_df$motif_score)& !is.na(mod_df$scr_topsc) & is.finite(mod_df$scr_topsc)
    fatac = !is.na(mod_df$FEAT1) & is.finite(mod_df$FEAT1) &!is.na(mod_df$FEAT2) & is.finite(mod_df$FEAT2) &!is.na(mod_df$FEAT3) & is.finite(mod_df$FEAT3)
    # two equivalent ways: either predict BRD4 or K27 as funciton of score and ATAC or predict score as a function of BRD4/K27 and ATAC
    linmod2=lm(FEAT2 ~ motif_score + FEAT1 + scr_topsc  ,data=mod_df[ fna & fatac ,])
    linmod3=lm(FEAT3 ~ motif_score + FEAT1 +scr_topsc  ,data=mod_df[ fna & fatac ,])
    linmod4=lm( motif_score ~ FEAT2 +FEAT1 + scr_topsc ,data=mod_df[ fna & fatac ,])
    linmod5=lm( motif_score ~ FEAT3 +FEAT1 + scr_topsc ,data=mod_df[ fna & fatac ,])
    df=data.frame(beta_K27_atcontr = summary(linmod2)$coef[2,co],
                  beta_mK27_atcontr = summary(linmod4)$coef[2,co],
                  beta_BRD4_atcontr = summary(linmod3)$coef[2,co],
                  beta_mBRD4_atcontr = summary(linmod5)$coef[2,co],TF=tf)
    colnames(df)[1:4] = paste0(lab,c('_K27_atcontr','_motY_K27_atcontr','_BRD4_atcontr','_motY_BRD4_atcontr'))
    return(df)
}

#5) match CM-leads and LOCAL-leads using caQTL leads

get_sample_LEAD_caQTLm= function(df,seedi = 1,MASTER_tab){
    #df = data frame of MASTER peak table 
    # seedi = set seed for reproducibility
    # MASTER_tab = the complete MASTER peak table that contains annotation data
    ## STEP 1 match P_caQTL 
    set.seed(seedi) # orig4
    
    pcaq1 = table(df$P_caQTL==1,df$Enh_type) # 
    rat_non1_vcm = pcaq1[1,'CM']/sum(pcaq1[,'CM'])
    sample_loc_req = round(round(pcaq1[2,'LOCAL'] /(1- rat_non1_vcm),0)*rat_non1_vcm,0)
    ## all cm lead enhancers with P_caQTL  smaller than 1 
    vcm_pcaqn1 = df$P_caQTL[df$Enh_type=='CM' & df$P_caQTL<1]
    ## same for local
    loc_mattab = df[df$Enh_type=='LOCAL' &df$P_caQTL<1,]
    loc_pcaqn1 = loc_mattab$P_caQTL
    keep_loc_pcaids = c()
    ## generate targets from CM pcaqtl<1 sample
    target = sample(vcm_pcaqn1,round(sample_loc_req,0),replace = T)
    for (i in 1:length(target)){
      w=which.min(abs(target[i]-loc_pcaqn1))
      loc_pcaqn1[w]=NA
      keep_loc_pcaids[i]=w
    }
    peak_ids_loc_pamat= loc_mattab$Peak[keep_loc_pcaids]
    dfn=df
    dfn$keep_PcaQTL_matched_ids = df$P_caQTL==1 | df$Enh_type %in% c('CM','no-effect') | df$Peak %in% peak_ids_loc_pamat
    
    # # STEP 2 GC CONT
    Lead_dfk = dfn[dfn$keep_PcaQTL_matched_ids | dfn$Enh_type=='no-effect',]
    #print(table(Lead_dfk$Enh_type))
    
    vcm_kid = which(Lead_dfk$Enh_type=='CM')
    loc_idsa = which(Lead_dfk$Enh_type=='LOCAL')
    nc_idsa = which(Lead_dfk$Enh_type=='no-effect')

    ## to sample
    
    vcm_gc = Lead_dfk$Lead_GC_scored[vcm_kid]

    add_vcm_gc_samp =  vcm_gc[sample(1:length(vcm_gc),(sum(Lead_dfk$Enh_type=='LOCAL')-300-length(vcm_kid)))]
    

    vcm_gcs= c(vcm_gc,add_vcm_gc_samp)
    #print(length(vcm_gcs))

    # locals
    loc_gcs = Lead_dfk$Lead_GC_scored[loc_idsa]

    loc_keep_idb = c()
    for (i in 1:length(vcm_gcs)){
      w=which.min(abs(vcm_gcs[i]-loc_gcs))
      loc_gcs[w]=NA
      loc_keep_idb[i]=loc_idsa[w]
    }
    # no-effects
    nc_gcs =Lead_dfk$Lead_GC_scored[nc_idsa]

    nc_keep_idb = c()
    for (i in 1:length(vcm_gcs)){
      w=which.min(abs(vcm_gcs[i]-nc_gcs))
      nc_gcs[w]=NA
      nc_keep_idb[i]=nc_idsa[w]
    }
    
    final_Lead_df = Lead_dfk[c(vcm_kid,loc_keep_idb,nc_keep_idb),]
    
    ## now keeping a dependent for each of the matched SNPs
    ## fisrt remove all leads that do not have an enhancer as neighbor or dep
    depanno = lapply(final_Lead_df$Peak,function(x){unique(MASTER_tab$dep_anno_tss500[which(MASTER_tab$Peak==x)])})
    f2_dannoE = sapply(depanno,function(x){
            ('Intron' %in% x | 'Distal' %in% x) })
    
    final_Lead_dff= final_Lead_df[ which(f2_dannoE) ,]     
   
    lpidk = unique(final_Lead_dff$Peak) # unique peak ids
    
    Dep_dfsglraw = do.call(rbind,lapply(1:length(lpidk),function(x){
      kt=Gaf_LEDE_IPtab[which(Gaf_LEDE_IPtab$Peak==lpidk[x] & Gaf_LEDE_IPtab$dep_anno_tss500 %in% c('Intron','Distal')),]
        # if no SNP in GM12878
        if (length(which(!is.na(kt$depGMSNP_RsID)))==0){ 
        # choose at random a dependent
        ret=kt[sample(dim(kt)[1],1),]
      }
        # if snp in GM12878
      if(length(which(!is.na(kt$depGMSNP_RsID)))>0){
        kkt = kt[!is.na(kt$depGMSNP_RsID),]
        if (length(which(kkt$ATACj_dep_pas>0))==0){
          ret=kkt[sample(dim(kkt)[1],1),]
        }
        if (length(which(kkt$ATACj_dep_pas>0))>0){
          kkkt = kkt[which(kkt$ATACj_dep_pas>0),]
          if (dim(kkkt)[1]==1){
            ret = kkkt
          }
          if (dim(kkkt)[1]>1){
            ret = kkkt[which.min(kkkt$depGMSNP_dist_to_dep_peak_center),]
          }
        }
      }
      
      return(ret)
    }))
    
    ## confrim that both pairs are indeadd enhancers
    f2_lannoE = Dep_dfsglraw$lead_anno_tss500 %in% c('Intron','Distal')
    f2_dannoE = Dep_dfsglraw$dep_anno_tss500 %in% c('Intron','Distal')
    

    Dep_dfsglEE = Dep_dfsglraw[ f2_dannoE ,]
    
    return(Dep_dfsglEE)
}
    
# 6) Distance match CM and independent enhancer pairs

get_dep_pairs = function(df,seed=1){
    # df = master peak data frame from LEAD enhancer GC and P_caQTL matching
    # seed = set seed for reproducibility
    # --------------
    #keep_df = df # all CMs without differentiating coordinated from uncoordinated
    keep_df = df[df$LD_coor_direction=='same' | df$Enh_type %in% c('LOCAL','no-effect'),] # only coordinated CMs 
    
    set.seed(seed) #orig   4
    vcm_kid = which(keep_df$Enh_type=='CM')
    loc_idsa = which(keep_df$Enh_type=='LOCAL')
    nc_idsa = which(keep_df$Enh_type=='no-effect')

    len_k = round(0.6*length(loc_idsa),0)
    len_k2 = round(0.75*length(loc_idsa),0)

    vcm_dists1 = abs(keep_df$Peak_center - keep_df$dep_Peak_center)[vcm_kid]
    
    add_vcm_dists_samp =  vcm_dists1[sample(1:length(vcm_dists1),(len_k-length(vcm_kid)))]
    vcm_dists = c(vcm_dists1,add_vcm_dists_samp)
   
    
    add_vcm_dists_samp2 =  vcm_dists1[sample(1:length(vcm_dists1),(len_k2-length(vcm_kid)))]

    vcm_dists2 = c(vcm_dists1,add_vcm_dists_samp2)

    
    # locals
    loc_dists = abs(keep_df$Peak_center - keep_df$dep_Peak_center)[loc_idsa]

    loc_keep_idb = c()
    for (i in 1:length(vcm_dists)){
      w=which.min(abs(vcm_dists[i]-loc_dists))
      loc_dists[w]=NA
      loc_keep_idb[i]=loc_idsa[w]
    }
    # no-effects
    nc_dists = abs(keep_df$Peak_center - keep_df$dep_Peak_center)[nc_idsa]

    nc_keep_idb = c()
    for (i in 1:length(vcm_dists2)){
      w=which.min(abs(vcm_dists2[i]-nc_dists))
      nc_dists[w]=NA
      nc_keep_idb[i]=nc_idsa[w]
    }

    final_keep_df = keep_df[c(vcm_kid,loc_keep_idb,nc_keep_idb),]
    return(final_keep_df)
}


# 7) enhancer context Z-score and SNP BS disruption of TFs in lead enhancers (run across all samples of LOCAL-to-CM matching)
generate_cont_caqtl_df = function(Enh_tab,peak_id_sc,dom_act_topZ_a,dom_pas_topZ_a,TFs_hoco,lab){
    # Enh_tab = master_peak table subsetted to a specific enhancer class (e.g. CM-Leads)
    # peak_id_sc = peak_ids as they occurr in the matrix obtaining top scores
    #dom_act_topZ_a,dom_pas_topZ_a = Z-transformed matrices of top scores across all TFs and enhancers (dim = n enh x TF motifs) for less and more accessible genotype (for mean and st-dev use randomly sampled and scored enhancers for each TF)
    #TFs _hoco = TF motif names from hocomoco (col names of scores matrix)
    # lab =  e.g. for which enhancer type the computation is done
    # ------
    
    dom_act_top_Zs =dom_act_topZ_a[match(Enh_tab$Peak,peak_id_sc),] # get scoring matrix in order of peak table
    dom_pas_top_Zs =dom_pas_topZ_a[match(Enh_tab$Peak,peak_id_sc),] # get scoring matrix in order of peak table
    ## snp disruption for a given TF
    N_TFBSdisr_per_TF = sapply(1:401,function(i){
        sum((dom_act_top_Zs[,i]-dom_pas_top_Zs[,i])>0,na.rm=T)
    })
    Nexp = length(Enh_tab[,1])
  
    dis_df=data.frame(TF=TFs_hoco,Ndis = N_TFBSdisr_per_TF/Nexp,label = lab)
    ## average score in less acc context
    avg_perTF_panE_pas_topZ = do.call(rbind,lapply(1:401,function(i){
        w=which(!((dom_act_top_Zs[,i]-dom_pas_top_Zs[,i])>0)) # remove enhancers were TF sits at SNP to assure  binding to enhancer context
        d=mean(dom_pas_top_Zs[w,i],na.rm=T)
        return(data.frame(avg_Z = d,label= lab))}))
    rownames(avg_perTF_panE_pas_topZ)=TFs_hoco
  
    comb_df = data.frame(frac_top_BS_incr=dis_df[,2],label = dis_df$label,
                       avg_cau_Zpas_disrexcl=avg_perTF_panE_pas_topZ[,1],
                       TF=TFs_hoco)
  
  
  return(comb_df)
  
}

# 8) same as 7 for neighbor/dependent elements 
generate_cont_caqtl_basic_df = function(Enh_tab,peak_id_sc,peak_Zscore_ref,TFs_hoco,lab){
    #peak_Zscore_ref = Z-scores for all TFs across reference seq of neighbor or dep enh (dim = n enh x TF motifs)
  
    
    pidk=match(Enh_tab$dep_Peak,peak_id_sc)
      ## average score in less acc context
      avg_perTF_panE_pas_topZ = do.call(rbind,lapply(1:401,function(i){
        d=mean(peak_Zscore_ref[pidk,i],na.rm=T)
        return(data.frame(avg_Z = d,label= lab))}))
      rownames(avg_perTF_panE_pas_topZ)=TFs_hoco
  
      comb_df = data.frame(label = avg_perTF_panE_pas_topZ$label,
                       avg_cau_Zpas_disrexcl=avg_perTF_panE_pas_topZ[,1],
                       TF=TFs_hoco)
  
  
      return(comb_df)  
}
    
