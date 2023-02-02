#' sampleSizefun
#'
#' Function that deal with determine training set when considering population structure in generating optimal training set.
#'
#' @param pickNum the number of genotypes provided by each subpopulation group; number vector
#' @param train_size training set size
#' @param subgpNum the number of genotypes of each subpopulation group
#' @noRd
#'
sampleSizefun = function(pickNum,train_size,subgpNum){
  while(sum(pickNum)<train_size){
    gap = sort(sample(1:length(subgpNum),train_size-sum(pickNum),replace = T)) #random pick up gps
    pickNum[gap] = pickNum[gap]+1}
  repeatInside = FALSE
  if(sum(pickNum>subgpNum)!=0) {repeatInside = TRUE}
  while(repeatInside){
    gpSite = which(pickNum>subgpNum)
    gap2 = sample(setdiff(1:length(subgpNum),gpSite),sum((pickNum[gpSite]-subgpNum[gpSite])),replace = T)
    pickNum[as.numeric(names(table(gap2)))] = pickNum[as.numeric(names(table(gap2)))]+table(gap2)
    pickNum[gpSite] = subgpNum[gpSite]
    if(sum(pickNum>subgpNum)==0) repeatInside = FALSE
  }
  return(pickNum)
}
