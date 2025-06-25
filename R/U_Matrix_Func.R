

#' Title all survive prob
#'
#' @param k number of kin to survive
#' @param s number of years kin survive (from birth)
#' @param U_mat mortality matrix
#'
#' @return probability that "k" individuals survive for "s" years

all_k_suv <- function(k,  n_start , n_stop, U_mat){

  n_start <- n_start + 1
  n_stop <- n_stop + 1
  n <- nrow(U_mat)

  if(n_stop < n){
    u <- diag(U_mat[-1,-ncol(U_mat)])
    u <- u[n_start : n_stop]%>%prod()
  }
  else{
    u <- diag(U_mat[-1,-ncol(U_mat)])
    u <- c(u[start : length(u)] , U_mat[n, n] )%>%prod()
  }

  f <- u^k
  return(f)
}

#' Title prob that j out of k survive (Binomial expansion)
#'
#' @param k number of kin
#' @param j number of kin survive (i.e., k - j die)
#' @param s number of years
#' @param U_mat mortality matrix
#'
#' @return probability that exactly "j out of k" individuals survive, over "s" years

k_to_j <- function(k, j, n_start , n_stop, U_mat){

  n <- nrow(U_mat)

  n_start <- n_start + 1
  n_stop <- n_stop + 1

  if(n_stop < n){
    u <- diag(U_mat[-1,-ncol(U_mat)])
    u <- u[n_start:n_stop]%>%prod()
  }
  else{
    u <- diag(U_mat[-1,-ncol(U_mat)])
    u <- c(u[n_start : length(u)] , U_mat[n,n] )%>%prod()
  }
  f = choose(k,j)*u^(j)*(1-u)^(k-j)
  return(f)

}

#' Title
#'
#' @param s years of survival
#' @param Q number of kin-number classes (from 0 to Q-1 max)
#' @param U_mat mortality matrix
#'
#' @return Eqs 6,8,9 in manuscript. A matrix which projects the pdf of newborns for "s" years and
#'  reshuffles them into kin-number classes conditional on death

U_kin_death <- function(n_start , n_stop, Q, U_mat){

  if(n_start >= n_stop){
    U_prob <- diag(Q)
  }
  else{
  U_prob <- matrix(0, Q, Q)
  U_prob[1,1] <- 1

  if(ncol(U_prob)>2){
    for(i in 2:(-1+nrow(U_prob))){
      for(j in (i+1):ncol(U_prob)){
        U_prob[i,j] <- k_to_j(j-1 ,i-1 , n_start , n_stop , U_mat)
        }
      }
    for(i in 2:nrow(U_prob)){
      U_prob[1,i] <- k_to_j(i-1, 0, n_start , n_stop , U_mat)
      U_prob[i,i] <- all_k_suv(i-1 , n_start , n_stop, U_mat)
      }
  }
  else{U_prob[2,2] <- all_k_suv(1 , n_start , n_stop, U_mat)
       U_prob[1,2] <- 1- U_prob[2,2]}
  }
  return(U_prob)
}
