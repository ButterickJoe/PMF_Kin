

#' Title
#'
#' @param r rate
#' @param k occurance
#'
#' @return Poisson distribution with parameters r,k

pois_fert_pdf <- function(r,k){
  pdf <- r^k*exp(-r)/factorial(k)
  return(pdf)
}

#' Title
#'
#' @param F_mat fertility matrix with age-specific rates on top-row
#' @param Q Maximum number of offspring
#'
#' @return Newborn probability distribution: Prob(number of newborns = 0,1,2,...,Q)

fert_dists <- function(F_mat, Q){

  n <- ncol(F_mat)
  fert_pdfs_list <- list()
  fert_pdfs_list <- lapply(1:n, function(x){
    pvec <- rep(0, Q)
    for(k in 1:length(pvec)){
      pvec[k] <- pois_fert_pdf(F_mat[1,x] , (k-1) )
    }
    fert_pdfs_list[[x]] <- pvec
  })
  return(fert_pdfs_list)
}
