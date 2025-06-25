

#' Title
#'
#' @param U_mat survival matirx with age-specific probs of moving age x to x+1 on subdiagonal
#' @param F_mat fertility matirx with age-specific fertility rate on top-row
#' @param x possible age of mothering

#' @return probability that a randomly selected newborn has a mother of age x

mothers_age <- function(U_mat, F_mat, x){
  A_proj <- U_mat + F_mat
  w_vec <- SD(A_proj)
  lam <- lambda(A_proj)
  pb <- F_mat[1, x]
  pop_freq <- w_vec[x]
  w1 <- w_vec[1]
  prob <- pb*pop_freq/(lam*w1)
  return(prob)
}
