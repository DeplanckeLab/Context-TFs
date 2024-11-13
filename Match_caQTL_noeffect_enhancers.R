
match_nc_to_cau=function(MASTER_tab,seedi=2){
    # MASTER_tab contains all the peak information of caQTL and no-effect peaks
    
    ### match cg 
    set.seed(seedi)
    cauidcg = which(MASTER_tab$causal_label=='caQTL')
    ncidcg = which(MASTER_tab$causal_label=='no-effect')

    ## repeat more caQTL enhancers to have a buffer when matching SNP distances below
    #print((length(ncidcg)-length(cauidcg)-150))
    cau_cgo =c(MASTER_tab$Lead_GC_scored[cauidcg],MASTER_tab$Lead_GC_scored[cauidcg][sample(1:length(cauidcg),(length(ncidcg)-length(cauidcg)-100))])
    ## randomize order of caQTL GC content
    cau_cg = cau_cgo[sample(1:length(cau_cgo),length(cau_cgo))]
    
    # noncausals
    nc_cg =MASTER_tab$Lead_GC_scored[ncidcg]

    # length
    len_cg = length(cau_cg)
    
    # sample
    nc_keep_idscg = c()
    for (i in 1:len_cg){
      w=which.min(abs(cau_cg[i]-nc_cg))
      nc_cg[w]=NA
      nc_keep_idscg[i]=ncidcg[w]
    }

    MASTER_tabCGm = MASTER_tab[c(cauidcg,nc_keep_idscg),]
    
    
    
    ### 2nd match SNP to peak center dist
    set.seed(seedi)
    cauid = which(MASTER_tabCGm$causal_label=='causal')
    ncid = which(MASTER_tabCGm$causal_label=='noncausal')


    # take all caus snp to peak center distances and match one to one
    cau_scdist =MASTER_tabCGm$Dist_SNP_to_peak_center[cauid]
    # noncausals
    nc_scdist =MASTER_tabCGm$Dist_SNP_to_peak_center[ncid]

    # length 
    len_scdk = length(cau_scdist)

    # sample
    nc_keep_idscdist = c()
    for (i in 1:len_scdk){
      w=which.min(abs(cau_scdist[i]-nc_scdist))
      nc_scdist[w]=NA
      nc_keep_idscdist[i]=ncid[w]
    }

    Cau_enh_tabscdCGm = Cau_enh_tabCGm[c(cauid,nc_keep_idscdist),]
    
    return(Cau_enh_tabscdCGm)
}