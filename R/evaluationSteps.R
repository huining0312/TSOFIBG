#' paral_trainSetSize
#'
#' Function deal with GBLUP model in evaluation steps.
#'
#' @param kinship  kinship matrix
#' @param optimal_trainSet dataframe; generated from previous steps; provide information of genotypes' selection Index
#' @param subpop dataframe; provide subpopulation information of genotypes
#' @param trainSize training set size that users desire
#' @param p_true we assume phenotype data collected from field is unknown; therefore we generate phenotype value by using MCMC.
#' @param g_true true genotypic values are inestimable; therefore we generate by using MCMC.
#' @noRd
#'
paral_trainSetSize <- function(kinship,optimal_trainSet,subpop,trainSize,p_true,g_true){
  if(!subpop){
    strategy = sort(optimal_trainSet[1:trainSize,3])
  }else{
    gpName = names(table(optimal_trainSet[,4]))
    for(i in gpName){
      nam <- paste(i)
      assign(nam, optimal_trainSet%>%filter(optimal_trainSet[,4]==i)%>%.$posIndex)
    }

    subgpNum = c()
    for(sp in gpName){
      subgpNum = c(subgpNum, length(get(sp)))
    }
    ratio_pp = subgpNum/sum(subgpNum)
    PPpickNum = floor(ratio_pp*trainSize)
    PickNumpp = sampleSizefun(pickNum = PPpickNum,train_size = trainSize,subgpNum = subgpNum)

    pick_pp = c()
    for(gpnum in 1:length(gpName)){
      nam <- paste("pp",gpName[gpnum],sep = "_")
      assign(nam, eval(parse(text=gpName[gpnum]))[c(1:PickNumpp[gpnum])])
      pick_pp = c(pick_pp, eval(parse(text = nam)))
    }
    strategy = sort(pick_pp)
  }

  kTR = kinship[strategy,strategy]
  p = p_true[strategy]

  fit=suppressWarnings(BGLR::BGLR(y=p,ETA=list(MRK=list(K=kTR,model="RKHS")),nIter=5000,burnIn=10^3,verbose=FALSE))

  g_hat=fit$ETA$MRK$u
  g_bp = kinship[,strategy]%*%solve(kTR)%*%g_hat
  row.names(g_bp) <- rownames(kinship)
  test_mat = cbind(g_true,g_bp)
  test_mat_sort = test_mat[order(-test_mat[,1]),]
  colnames(test_mat_sort)[2]<-"g_bp"
  sort_g_true = test_mat_sort[,1]
  #k=1
  ndcg_1 = DCG_k_pre(g_bp,g_true,1)/DCG_k_true(sort_g_true,1)
  #k=5
  ndcg_5 = DCG_k_pre(g_bp,g_true,5)/DCG_k_true(sort_g_true,5)
  #k=10
  ndcg_10 = DCG_k_pre(g_bp,g_true,10)/DCG_k_true(sort_g_true,10)
  #k=20
  ndcg_mean = mean_NDCGk(g_bp,g_true,10,sort_g_true)

  NDCGVals = c(ndcg_1,ndcg_5,ndcg_10,ndcg_mean)
  return(NDCGVals)
}

#' evalsteps
#'
#' Main function deal with GBLUP model in evaluation steps.(consider different heritability settings)
#'
#' @param optimal_trainSet  dataframe; generated from previous steps; provide information of genotypes' selection Index
#' @param kinship kinship matrix
#' @param nsim Number of repetitions
#' @param subpop Boolean value; set T/TRUE if consider population structure
#' @param desireH we assume phenotype data collected from field is unknown; therefore we generate phenotype value by using MCMC.
#' @param mu true genotypic values are inestimable; therefore we generate by using MCMC.
#' @param sg variance of genotypic values; default is 25
#' @param n the total number of candidates
#' @param desireDelta The proportion of the candidate dataset that users would like to use as the training set.
#' @noRd
#'
evalsteps <- function(optimal_trainSet,kinship,nsim=2000,subpop,desireH,mu=100,sg=25,n,desireDelta){
  trainSize = round(desireDelta*n)
  innerList1 = matrix(0,length(trainSize),4)
  innerList2 = matrix(0,length(trainSize),4)
  EIResult = lapply(1:length(desireH),function(x) innerList1)
  RE_EIResult = lapply(1:length(desireH),function(x) innerList2)
  names(EIResult) = names(RE_EIResult) <- paste("h",desireH,sep = "_")

  gT = mvrnorm(nsim,rep(0,n),sg*kinship)
  for(i in 1:length(desireH)){
    h = desireH[i]
    se=sg*(1-h)/h
    pT = mu+gT+mvrnorm(nsim,rep(0,n),se*diag(1,nrow = n,ncol = n))
    temp <-
      foreach(s=1:nsim, .combine = "+") %:%
        foreach(t=1:length(trainSize), .combine="rbind") %dopar% {
          paral_trainSetSize(kinship,optimal_trainSet,subpop,trainSize[t],p_true = pT[s,],g_true=gT[s,])
        }
    EIResult[[i]] = round((temp/nsim),4)
    RE_EIResult[[i]] = apply(EIResult[[i]],2,function(x) x/x[length(trainSize)])
    colnames(EIResult[[i]]) = c("k1","k5","k10","meank10")
    colnames(RE_EIResult[[i]]) = c("RE_k1","RE_k5","RE_k10","RE_meank10")
    rownames(EIResult[[i]]) = rownames(RE_EIResult[[i]]) = paste("N",trainSize,sep = "")
  }
  return(list("EIResult"=EIResult,"RE_EIResult"=RE_EIResult))
}
