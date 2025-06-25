


#' Title
#'
#' @param age age of individual at reproduction
#' @param Q upper bound for number of offspring
#' @param F_mat fertility matrix with age-specific rates on top-row
#'
#' @return The matrix F as in Eq 12 of manuscript. Projects a pmf of offrping,

Q_matrix <- function(age, Q, F_mat){

  fert_dist <- fert_dists(F_mat, Q)

  phi_ <- fert_dist[[(age+1)]]
  if(length(phi_)!=Q ){"error"}
  list_phi_ <- list()
  for(i in 1 : Q){
    list_phi_[[i]] <- phi_
  }

  Q_mat <- matrix(0, Q, Q)
  for(i in 2:(ncol(Q_mat)) ){
    Q_mat[,i] <- convoluion_nth((i-1), list_phi_[1:(i-1)])
    #Q_mat[,i] <- fourier_fft_nth(i, list_phi_[1:i])
  }
  Q_mat[,1] <- c(1, rep(0, Q-1))

  return(Q_mat)
}
