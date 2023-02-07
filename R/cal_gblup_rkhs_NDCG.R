#' GBLUP_RKHS_NDCG
#'
#' Calculate RE_NDCG and RE_mean_NDCG of bayesian RKHS and RKHS
#'
#' @param GV_gen genotypic value estimated by collected phenotype data (estimate based on GBLUP model); please add species ID at the first column; treat as true phenotype data
#' @param gblup_result predicted genotypic values obtained from bayesian RKHS method
#' @param rkhs_result predicted genotypic values obtained from RKHS method
#' @param TRsize training set size
#' @noRd
#'
GBLUP_RKHS_NDCG <- function(GV_gen,gblup_result,rkhs_result,TRsize){
  # --- calculate NDCG --- #
  gblup_ndcg_res = lapply(1:(ncol(GV_gen)-1),function(x) data.frame("TRsize" = TRsize, "k1" = rep(0,length(TRsize)),
                                                                    "k5" = rep(0,length(TRsize)),"k10" = rep(0,length(TRsize)),
                                                                    "mean@k10" = rep(0,length(TRsize))))

  rkhs_ndcg_res = lapply(1:(ncol(GV_gen)-1),function(x) data.frame("TRsize" = TRsize, "k1" = rep(0,length(TRsize)),
                                                                   "k5" = rep(0,length(TRsize)),"k10" = rep(0,length(TRsize)),
                                                                   "mean@k10" = rep(0,length(TRsize))))
  names(gblup_ndcg_res) = names(rkhs_ndcg_res)= colnames(GV_gen)[-1]


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
