#' AugEI_fun
#'
#' Function deal with GBLUP model and calculate AugEI value for each genotype
#' This can go on for multiple lines.
#'
#' @param CV_gpnumber  number of folds for cross validation
#' @param candi number vector; length is equal to total candidates
#' @param gp1 the number of group for the size num1
#' @param n the total number of candidates
#' @param num1 num1+num2=n
#' @param num2 num1+num2=n
#' @param kinship kinship matrix
#' @param p_sim phenotype data for genotypes (calculate by GBLUP model)


AugEI_fun = function(CV_gpnumber,candi,gp1,n,num1,num2,kinship,p_sim){
  EIfinal = data.frame(meanEI=rep(0,n))
  EImat = matrix(NA, nrow = n, ncol = CV_gpnumber)
  rownames(EImat) = rownames(EIfinal) = rownames(kinship)
  #
  part = list()
  pool = candi

  for ( i in 1:CV_gpnumber){
    if(i < gp1){
      numSP = num1
      part = append(part,list(sample(pool,numSP)))
      pool = setdiff(pool,part[[i]])
    }else{
      numSP = num2
      part = append(part,list(sample(pool,numSP)))
      pool = setdiff(pool,part[[i]])
    }

  }

  for ( i in 1:length(part)){
    ini = part[[i]]
    te=setdiff(1:n,ini)

    p_new = p_sim
    # subgroup
    k11=kinship[ini,ini]
    k12=kinship[ini,te]
    k21=kinship[te,ini]
    k22=kinship[te,te]

    p = p_new[ini]
    ### BGLR gblup model ###
    fit=suppressWarnings(BGLR::BGLR(y=p,ETA=list(MRK=list(K=k11,model="RKHS")),nIter=5000,burnIn=10^3,verbose=FALSE))
    sg1=fit$ETA$MRK$varU
    se1=fit$varE
    g1=fit$ETA$MRK$u

    ### conditional mean and conditional variance ###
    mu.g2=k21%*%solve(k11)%*%g1
    var.g2=diag(sg1*(k22-k21%*%solve(k11)%*%k12))

    ### expected improvement ###
    uf=g1-diag(sqrt(sg1)*k11)
    muf=max(uf)
    mx=g1[which(uf==muf)]
    z=(mu.g2-mx)/sqrt(var.g2)
    ei=((mu.g2-mx)*pnorm(z)+sqrt(var.g2)*dnorm(z))*(1-sqrt(se1)/sqrt(var.g2+se1))

    for(index in 1:nrow(ei)){
      EImat[which(rownames(EImat)==rownames(ei)[index]),i] <- ei[index]
    }
  }
  for (k in 1:nrow(EImat)){
    EIfinal[k,1] = EIfinal[k,1] + mean(EImat[k,], na.rm=T)
  }
  return(EIfinal)
}


CV_process = function(kinship,nsim,CV_gpnumber,h,n,sg,mu){
  # --- parameters --- #
  candi = 1:n
  se=sg*(1-h)/h ### noise variance ###
  gp1 = floor(n/round(n/CV_gpnumber))
  num1 = round(n/CV_gpnumber)
  num2 = floor((n-(round(n/CV_gpnumber)*(floor(n/round(n/CV_gpnumber))-1))))/2
  g=MASS::mvrnorm(nsim,rep(0,n),sg*kinship)
  p_new = mu+g+MASS::mvrnorm(nsim,rep(0,n),se*diag(1,nrow = n,ncol = n))

  result <- foreach(i=1:nsim, .combine=cbind) %dopar% {
    AugEI_fun(CV_gpnumber,candi,gp1,n=n,num1,num2,kinship=kinship,p_new[i,])
  }

  final_dat = data.frame(augEI=apply(result, 1, mean))
  rank = order(final_dat[,1], decreasing = T)
  sort_result = data.frame(speciesName=rownames(final_dat)[rank],aug.EI=final_dat[rank,],posIndex=rep(0,n))
  for(i in 1:nrow(sort_result)){
    sort_result[sort_result[,1]==rownames(kinship)[i],3] = i
  }
  return(sort_result)
}


#' gen_sel_Index
#'
#' Main function.This function will generate the optimal training set for the given dataset.
#' Users can input the desiring sample size and heritability to obtain the NDCG and relative NDCG results based on the settings and the optimal training set.
#' The index of the suggested training set will also be returned.
#'
#' @param kinship kinship matrix
#' @param nOpsim  Number of repetitions for cross-validation
#' @param nEvalsim Number of repetitions for evaluate the performance of optimal training set
#' @param CV_gpnumber number of folds for cross validation
#' @param h heritability parameters for obtaining the optimal training set; default is 0.5
#' @param n the total number of candidates
#' @param sg variance of genotypic values; default is 25
#' @param mu general mean; default is 100
#' @param subpopTag Information of population structure. Input dataset should include genotypes' name and subpopulation tags. If the dataset don not consider the population structure, the parameters is allowed to miss.
#' @param desireH The heratibility parameters for user's desiring traits. This would be used in the simulation step.
#' @param desireDelta The proportion of the candidate dataset that users would like to use as the training set.
#' @export
#'
gen_sel_Index = function(kinship,nOpsim=1000,nEvalsim=2000,CV_gpnumber=5,h=0.5,n,sg=25,mu=100,subpopTag,desireH = c(0.5),desireDelta=c(1/5,1/3,2/3)){
  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)

  sort_result = CV_process(kinship,nOpsim,CV_gpnumber,h,n,sg,mu)

  if(missing(subpopTag)){
    simuRes = evalsteps(optimal_trainSet=sort_result,kinship=kinship,nsim=nEvalsim,subpop = F,desireH=desireH,mu=mu,sg=sg,desireDelta=desireDelta,n=n)
  }else{
    cat(paste("evaluate step processing..."),sep = "\n")
    sort_result$tag = rep(NA,n)
    for(i in 1:n){
      sort_result[i,4] <- as.character(subpopTag[subpopTag[,1] == sort_result[i,1],2])
    }
    simuRes = evalsteps(optimal_trainSet=sort_result,kinship=kinship,nsim=nEvalsim,subpop = T,desireH=desireH,mu=mu,sg=sg,n=n,desireDelta = desireDelta)
  }
  return(list(selIndex=sort_result,simuNDCG=simuRes[[1]]))
  doParallel::stopImplicitCluster()
}
