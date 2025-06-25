


#' Title
#'
#' @param y age of Focal
#' @param s0 age of mother when Focal is y
#' @param U_mat matrix of mortality rates
#' @param F_mat matrix of fertility rates
#' @param Q number of kin classes 0,...,Q-1

#' @return an age-specific pdf (data frame) with number = 0, 1, ..., Q-1 and resp. probs of kin

unconditional_mother_pdf <- function(y, s0, U_mat, F_mat, Q){

  n_min <- which(mothers_age(U_mat, F_mat)>0)[1]
  n_max <- which(mothers_age(U_mat, F_mat)>0)[length(which(mothers_age(U_mat, F_mat)>0))]
  b1 <- s0 - y
  if(b1 >= (n_min) & b1 <= (n_max)){
  p <- mothers_age(U_mat, F_mat, b1)
  vv <- rep(0, 2)
  vv[1] <- 1-p
  vv[2] <- p
  mother <-  U_kin_death(b1, s0-1, 2, U_mat) %*% vv
  mother_df <- data.frame(number = c(0,1),
                            prob = mother,
                            age = s0)}
  else{mother_df <- data.frame(number = c(0,1),
                               prob = c(1,0),
                               age = s0)}

  return(mother_df)

}



#' Title
#'
#' @param y age of Focal
#' @param s0 age of gran when Focal is y
#' @param U_mat mortality matrix
#' @param F_mat fertility matrix
#' @param Q number of kin number classses
#'
#' @return a pdf in data frame output with numbers and resp probs

unconditional_grandmother_pdf <- function(y, s0, U_mat, F_mat, Q){
  n <- nrow(F_mat)
  fertile_values <- unlist(lapply(1:n, function(i){mothers_age(U_mat, F_mat, i)}))
  fertile_ages <- which(fertile_values>0)
  #fertile_values <- fertile_values[fertile_ages]
  sum <- s0 - y
  grandmother <- rep(0, Q)
  pp <- 0
  for(i in fertile_ages){
    pm <- fertile_values[i]
    for(j in fertile_ages){
      pg <- fertile_values[j]
      if(i + j == sum){
        prob <- pm*pg
        pp <- pp + prob
        vv <- rep(0, Q)
        vv[1] <- 0
        vv[2] <- prob
        grandmother <-  grandmother + U_kin_death(j, i + j + y - 1 , Q, U_mat) %*% vv
      }
    }
  }
  pnogran <- 1 - pp
  grandmother[1] <- pnogran
  grandmother_df <- data.frame(number = seq(0, Q-1),
                          prob = grandmother,
                          age_gran = s0)

  return(grandmother_df)
}

