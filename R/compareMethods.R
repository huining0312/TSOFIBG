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

#' EstimationM
#'
#' compare bayesian RKHS (GBLUP model), RKHS, random forest and McRank
#'
#' @param snp_matrix the SNP marker matrix that has been organized (finishing quality control, imputation and standardized)
#' @param opt_trainSet the number of group for the size num1
#' @param pheno collected phenotype data; please include species ID and the column names; if data have population structure, please include subpopulations information in the dataframe
#' @param delta training set size = delta *n; default:1/5,1/4,1/3,1/2,2/3,3/4,1
#' @param n the number of genotypes in the dataset
#' @export
#'
EstimationM <- function(snp_matrix,opt_trainSet,pheno,delta=c(1/5,1/4,1/3,1/2,2/3,3/4,1),n=nrow(kin)){
  kin = (as.matrix(snp_matrix)%*%t(as.matrix(snp_matrix)))/ncol(as.matrix(snp_matrix))

  TRsize = round(delta*n)
  GV_gen = simu_gv(pheno,kin)
  res_gblup_rkhs = gblup_rkhs_process(opt_trainSet=sort_data_mid,GV_gen,subpopTag=T,TRsize=TRsize)
  gblup_result = res_gblup_rkhs[[1]]
  rkhs_result = res_gblup_rkhs[[2]]
  picklist = res_gblup_rkhs[[3]]
  res_ndcg = GBLUP_RKHS_NDCG(GV_gen,gblup_result,rkhs_result,TRsize=TRsize)

  # --- python --- #
  #reticulate::virtualenv_list()
  reticulate::py_install("pandas")
  reticulate::py_install("scikit-learn")
  pd <- NULL
  sklearn <- NULL
  pd <- reticulate::import("pandas")
  sklearn <- reticulate::import("sklearn")

  RF <- sklearn$ensemble$RandomForestRegressor
  RFC <- sklearn$ensemble$RandomForestClassifier
  reticulate::py_run_file(system.file("python", "McRank.py", package = "TSOFIBG"))
  reticulate::py_run_file(system.file("python", "NDCG_fun.py", package = "TSOFIBG"))
  reticulate::py_run_file(system.file("python", "mainFunRF.py", package = "TSOFIBG"))
  reticulate::py_run_file(system.file("python", "mainFunMcRank.py", package = "TSOFIBG"))
  randomforest_res <- reticulate::py$main_fun_RF(snp_matrix_up,TRsize,GV_gen,picklist,RF)
  mcrank_res <- reticulate::py$main_fun_Mcrank(snp_matrix_up,TRsize,GV_gen,picklist,RFC)
  #return(picklist)
  return(list("GBLUP_NDCG"=res_ndcg[[1]],"RKHS_NDCG"=res_ndcg[[2]],"RF"=randomforest_res,"McRank"=mcrank_res))
}
