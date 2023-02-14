#' GBLUP_RKHS_NDCG
#'
#' Calculate RE_NDCG and RE_mean_NDCG of bayesian RKHS and RKHS
#'
#' @param GV_gen genotypic value estimated by collected phenotype data (estimate based on GBLUP model); please add species ID at the first column; treat as true phenotype data
#' @param gblup_result predicted genotypic values obtained from bayesian RKHS method
#' @param rkhs_result predicted genotypic values obtained from RKHS method
#' @param TRsize training set size
#' @param NDCG_K the ranking user desire; default is k=1. k=5 and k=10; please input a number vector
#' @param mean_NDCG_K please input the desired average NDCG value; default is mean_NDCG for k =10
#' @noRd
#'
GBLUP_RKHS_NDCG <- function(GV_gen,gblup_result,rkhs_result,TRsize,NDCG_K=c(1,5,10),mean_NDCG_K=c(10)){
  # --- calculate NDCG --- #
  # --- generate container --- #

  build_df = function(TRsize,NDCG_K,mean_NDCG_K,col_name){
    row_num = length(c(NDCG_K,mean_NDCG_K))+1
    t = matrix(NA,nrow =row_num, ncol = length(col_name))
    dimnames(t) = list(c("TRsize",paste("k",c(NDCG_K),sep=""),paste("mean_NDCG@k",meanK,sep = "")),c(col_name))
    return(as.data.frame(t))
  }

  gblup_ndcg_res = lapply(1:(ncol(GV_gen)-1),function(x) build_df(TRsize,NDCG_K,mean_NDCG_K,col_name=colnames(GV_gen)[-1]))

  rkhs_ndcg_res = lapply(1:(ncol(GV_gen)-1),function(x) build_df(TRsize,NDCG_K,mean_NDCG_K,col_name=colnames(GV_gen)[-1]))

  #names(gblup_ndcg_res) = names(rkhs_ndcg_res)= colnames(GV_gen)[-1]


  for(i in 1:length(TRsize)){
    for (trait in names(gblup_result)){

      g_gblup = gblup_result[[trait]][,i]
      g_rkhs = rkhs_result[[trait]][,i]
      g_true = GV_gen[,trait]
      g_true_sort = sort(GV_gen[,trait], decreasing = T)

      # --- gblup --- #
      #k=1
      ndcg_1 = DCG_k_pre(g_gblup,g_true,1)/DCG_k_true(g_true_sort,1)
      #k=5
      ndcg_5 = DCG_k_pre(unname(g_gblup),g_true,5)/DCG_k_true(g_true_sort,5)
      #k=10
      ndcg_10 = DCG_k_pre(g_gblup,g_true,10)/DCG_k_true(g_true_sort,10)
      #meank
      ndcg_mean = mean_NDCGk(g_gblup,g_true,10,g_true_sort)

      gblup_ndcg_res[[trait]][i,2:5] = round(c(ndcg_1,ndcg_5,ndcg_10,ndcg_mean),4)

      # --- rkhs --- #
      #k=1
      ndcg_1 = DCG_k_pre(g_rkhs,g_true,1)/DCG_k_true(g_true_sort,1)
      #k=5
      ndcg_5 = DCG_k_pre(g_rkhs,g_true,5)/DCG_k_true(g_true_sort,5)
      #k=10
      ndcg_10 = DCG_k_pre(g_rkhs,g_true,10)/DCG_k_true(g_true_sort,10)
      #meank
      ndcg_mean = mean_NDCGk(g_rkhs,g_true,10,g_true_sort)

      rkhs_ndcg_res[[trait]][i,2:5] = round(c(ndcg_1,ndcg_5,ndcg_10,ndcg_mean),4)

    }
  }
  return(list("GBLUP_NDCG"=gblup_ndcg_res,"RKHS_NDCG"=rkhs_ndcg_res))
}
