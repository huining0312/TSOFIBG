#' gblup_rkhs_process
#'
#' Function generate genotypic value by using bayesian RKHS (GBLUP model) and RKHS
#'
#' @param opt_trainSet optimal training set
#' @param GV_gen genotypic value estimated by collected phenotype data (estimate based on GBLUP model); please add species ID at the first column
#' @param TRsize training set size
#' @param subpopTag boolean; whether the dataset has string population structure
#' @param trait_num the number of traits
#' @noRd

gblup_rkhs_process <- function(opt_trainSet,GV_gen,TRsize,subpopTag,trait_num){

  gpName = names(table(opt_trainSet$tag))
  for(i in gpName){
    nam <- paste(i)
    assign(nam, opt_trainSet%>%filter(tag==i)%>%.$name)
  }
  subgpNum = c()
  for(i in 1:length(gpName)){
    subgpNum = c(subgpNum,length(get(gpName[i])))
  }
  # proportion
  ratio_pp = subgpNum/sum(subgpNum)
  pickList = list()

  TRsize = round(delta*n)

  if(missing(trait_num)){
    trait_num = ncol(GV_gen)-1
  }else{
    trait_num = trait_num
  }

  gblup_result = lapply(1:trait_num,function(x) matrix(0,nrow(kin),length(TRsize)))
  names(gblup_result) = colnames(GV_gen)[colnames(GV_gen)!="ID"]
  rkhs_result = lapply(1:trait_num,function(x) matrix(0,nrow(kin),length(TRsize)))
  names(rkhs_result) = colnames(GV_gen)[colnames(GV_gen)!="ID"]


  for(trait in names(gblup_result)){
    dimnames(gblup_result[[trait]]) = list(colnames(kin),paste("N",TRsize,sep = ""))
    dimnames(rkhs_result[[trait]]) = list(colnames(kin),paste("N",TRsize,sep = ""))

    for(s in 1:length(TRsize)){
      size = TRsize[s]
      cat(paste(size),sep = "\n")
      if(!subpopTag){
        pick_pp = opt_trainSet$name[1:size]
      }else{
        #gblup model
        PPpickNum = floor(ratio_pp*size)
        PickNumpp = sampleSizefun(pickNum = PPpickNum,train_size = size,subgpNum = subgpNum)
        pick_pp = c()
        for(gpnum in 1:length(gpName)){
          if(PickNumpp[gpnum]!=0){
            nam <- paste("pp",gpName[gpnum],sep = "_")
            assign(nam, eval(parse(text=gpName[gpnum]))[c(1:PickNumpp[gpnum])])
          }
          pick_pp = c(pick_pp, eval(parse(text = nam)))
        }
        pickList[[s]] = pick_pp
      }
      trainset = sort(pick_pp)
      p = GV_gen[GV_gen$ID%in%trainset,trait] #sorted
      kTR = kin[rownames(kin)%in%trainset,rownames(kin)%in%trainset]
      fit_gblup=suppressWarnings(BGLR::BGLR(y=p,ETA=list(MRK=list(K=kTR,model="RKHS")),nIter=5000,burnIn=10^3,verbose=FALSE))

      g_hat <- fit_gblup$ETA$MRK$u
      g_gblup = kin[,rownames(kin)%in%trainset]%*%solve(kTR)%*%g_hat

      gblup_result[[trait]][,s] = g_gblup
      # RKHS
      df = data.frame(y=p,gid=rownames(kTR))
      ans = rrBLUP::kin.blup(df,geno="gid",pheno="y",K=kTR)
      g_rkhs = ans$g
      g_bp_rkhs = kin[,rownames(kin)%in%trainset]%*%solve(kTR)%*%g_rkhs
      rkhs_result[[trait]][,s] = g_bp_rkhs
    }

  }

  return(list("gblup_result"=gblup_result,"rkhs_result"=rkhs_result, "genoList"=pickList))
}
