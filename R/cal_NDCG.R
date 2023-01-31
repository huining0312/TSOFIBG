#' DCG_k_pre
#'
#' calculate DCG score for predicted ranking
#'
#' @param g_bp genotypic value predicted by GBLUP model
#' @param g_true true genotypic value (unsorted)
#' @param k the ranking of the genotypes; default is 1, 5, 10 and the mean of top 10
#' @noRd
#'
DCG_k_pre <- function(g_bp,g_true,k){
  if(missing(k)) k <- length(g_bp)
  stopifnot(k <= length(g_bp))
  ans = c()
  for (i in 1:k){
    target = order(-g_bp)[i]
    ans = c(ans,g_true[target]/log(i + 1, 2))
  }
  return(sum(ans))
}

#' DCG_k_true
#'
#' calculate DCG score for ideal ranking
#'
#' @param true_gen true genotypic value sorted in the decreasing order
#' @param k the ranking of the genotypes; default is 1, 5, 10 and the mean of top 10
#' @noRd
#'
DCG_k_true <- function(true_gen,k){
  if(missing(k)) k <- length(true_gen)
  stopifnot(k <= length(true_gen))
  return(sum(true_gen[seq(k)]/log(seq(k) + 1, 2)))
}


#' mean_NDCGk
#'
#' calculate the mean of DCG score for top k ranking
#'
#' @param g_bp genotypic value predicted by GBLUP model
#' @param g_true true genotypic value(unsorted)
#' @param k the ranking of the genotypes; default is 1, 5, 10 and the mean of top 10
#' @param true_gen true genotypic value sorted in the decreasing order
#' @noRd
#'
mean_NDCGk <- function(g_bp,g_true,k,true_gen){
  temp = c()
  for (i in 1:k){
    temp = c(temp,DCG_k_pre(g_bp,g_true,i)/DCG_k_true(true_gen,i))
  }
  ans = mean(temp)
  return(ans)
}

