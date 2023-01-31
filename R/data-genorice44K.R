#' A genomic relationship data (kinship matrix) of 44K rice (data with strong population structure).
#'
#' The kinship matrix is calculated by the function: K = X%*%t(X)/ncol(X), where X stands for standardized SNP marker matrix.
#' The original dataset included 413 genotypes, but we used only 301 genotypes.
#' Report ...
#'
#' @format
#' A matrix with 301 rows and columns:
#' The row names and column names referred to the name of genotypes.
#' @source <https://doi.org/10.1038/ncomms1467>
"geno_rice44K"
