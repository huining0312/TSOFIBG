#' simu_gv
#'
#' estimate genotypic value from collected phenotype data; the generating genotypic value will be treated as true phenotype data at the following steps
#'
#' @param pheno collected phenotype data; please include species ID and the column names
#' @param kin kinship matrix
#' @noRd
#'

simu_gv <- function(pheno,kin){
  GV_gen = data.frame(matrix(NA,nrow = nrow(pheno), ncol = ncol(pheno)))
  rownames(GV_gen) = rownames(kin)
  colnames(GV_gen) = colnames(pheno)
  GV_gen$ID = pheno$ID
  for(trait in colnames(pheno)[-1]){
    #gblup model
    trainset = pheno$ID
    p = pheno[pheno$ID%in%trainset,trait] #sorted
    kTR = kin[rownames(kin)%in%trainset,rownames(kin)%in%trainset]
    fit_gblup=suppressWarnings(BGLR::BGLR(y=p,ETA=list(MRK=list(K=kTR,model="RKHS")),nIter=5000,burnIn=10^3,verbose=FALSE))

    g_hat <- fit_gblup$ETA$MRK$u
    GV_gen[,trait] = g_hat
  }
  return(GV_gen)
}

#' gen_table
#'
#' generate result table
#'
#' @param gblup_res methods: GBLUP; genotypic values generated from collected phenotype data and the optimal training set
#' @param rkhs_res methods: RKHS; genotypic values generated from collected phenotype data and the optimal training set
#' @param rf_Res methods: random forest; genotypic values generated from collected phenotype data and the optimal training set
#' @param mcrank_res methods: McRank; genotypic values generated from collected phenotype data and the optimal training set
#' @param k the first k rank
#' @noRd
#'

gen_table = function(gblup_res, rkhs_res, rf_Res, mcrank_res, k){
  # trait
  df_final = sapply(names(gblup_res), function(x) NULL)
  for( col in 1:length(gblup_res)){
    gblup_t = gblup_res[[col]]
    rkhs_t = rkhs_res[[col]]
    rf_t = rf_Res[[col]]
    mcrank_t = mcrank_res[[col]]

    # training set size
    df_list_k = sapply(colnames(gblup_res[[col]]), function(x) NULL)

    for(size in 1:ncol(gblup_t)){

      gblup_sort = gblup_t[order(gblup_t[,size], decreasing = T),size]
      names(gblup_sort) = substring(names(gblup_sort),2)
      gblup_sort_df = data.frame(gblup_sort)
      gblup_sort_df$score = 1:nrow(gblup_sort_df)

      rkhs_sort = rkhs_t[order(rkhs_t[,size], decreasing = T),size]
      names(rkhs_sort) = substring(names(rkhs_sort),2)

      rf_sort = rf_t[order(rf_t[,size], decreasing = T),size]
      names(rf_sort) = rownames(rf_t)[order(rf_t[,size], decreasing = T)]

      mcrank_sort = mcrank_t[order(mcrank_t[,size], decreasing = T),size]
      names(mcrank_sort) = rownames(mcrank_t)[order(mcrank_t[,size], decreasing = T)]

      final_score_df = gblup_sort_df
      final_score_df$rkhs_score = match(rownames(final_score_df),names(rkhs_sort))
      final_score_df$rf_score = match(rownames(final_score_df),names(rf_sort))
      final_score_df$mcrank_score = match(rownames(final_score_df),names(mcrank_sort))

      final_score_df$sumScore = apply(final_score_df[,2:5], 1, sum)
      final_score_df = final_score_df[order(final_score_df$sumScore),]

      df_k = list()

      for(rank in k){
        df_k[paste("k=",rank,sep = "")] = list(rownames(final_score_df)[1:rank])
        df_list_k[colnames(gblup_res[[col]])[size]] = list(df_k)
      }

    }

    df_final[col] = list(df_list_k)
  }

  return(df_final)
}



#' EstimationM
#'
#' compare bayesian RKHS (GBLUP model), RKHS, random forest and McRank
#'
#' @param snp_matrix the SNP marker matrix that has been organized (finishing quality control, imputation and standardized)
#' @param opt_trainSet the number of group for the size num1
#' @param pheno collected phenotype data; please include species ID and the column names; if data have population structure, please include subpopulations information in the dataframe
#' @param delta training set size = delta *n; default:1/5,1/4,1/3,1/2,2/3,3/4,1
#' @param n the number of genotypes in the dataset
#' @param k the first k individuals (after sorting)
#' @export
#'
EstimationM <- function(snp_matrix,opt_trainSet,pheno,delta=c(1/5,1/4,1/3,1/2,2/3,3/4,1),n=nrow(kin),k=c(1,5,10)){
  kin = (as.matrix(snp_matrix)%*%t(as.matrix(snp_matrix)))/ncol(as.matrix(snp_matrix))

  TRsize = round(delta*n)
  GV_gen = simu_gv(pheno,kin)
  res_gblup_rkhs = gblup_rkhs_process(opt_trainSet=sort_data_mid,GV_gen,subpopTag=T,TRsize=TRsize)
  gblup_result = res_gblup_rkhs[[1]]
  rkhs_result = res_gblup_rkhs[[2]]
  picklist = res_gblup_rkhs[[3]]
  #res_ndcg = GBLUP_RKHS_NDCG(GV_gen,gblup_result,rkhs_result,TRsize=TRsize)

  # --- python --- #
  #reticulate::virtualenv_list()
  reticulate::py_install("pandas")
  reticulate::py_install("scikit-learn")
  sklearn <- NULL
  sklearn <- reticulate::import("sklearn")
  reticulate::py_run_string('import pandas')

  RF <- sklearn$ensemble$RandomForestRegressor
  RFC <- sklearn$ensemble$RandomForestClassifier
  reticulate::py_run_file(system.file("python", "McRank.py", package = "TSOFIBG"))

  # return a list
  reticulate::py_run_file(system.file("python", "mainFunRF.py", package = "TSOFIBG"))
  reticulate::py_run_file(system.file("python", "mainFunMcRank.py", package = "TSOFIBG"))
  randomforest_res <- reticulate::py$main_fun_RF(snp_matrix_up,TRsize,GV_gen,picklist,RF)
  mcrank_res <- reticulate::py$main_fun_Mcrank(snp_matrix_up,TRsize,GV_gen,picklist,RFC)

  # display table
  #final_res = gen_table(gblup_result, rkhs_result, randomforest_res, mcrank_res, k)

  #return(final_res)
  return(list("gblup"=gblup_result,"rkhs"=rkhs_result,"rf"=randomforest_res,"mcrank"=mcrank_res))
  #return(list("GBLUP_NDCG"=res_ndcg[[1]],"RKHS_NDCG"=res_ndcg[[2]],"RF"=randomforest_res,"McRank"=mcrank_res))
}
