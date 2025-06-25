



#' Title calculate the kin-number distribution of Focal's cousin descending from younger aunts
#'
#' @param y age of Focal
#' @param s2 age of Focal's cousin from younger aunt, when Focal is y
#' @param U_mat matrix with survival probs
#' @param F_mat matrix of fertility rates
#' @param Q number of kin classes (0,1,2,...,Q)
#'
#' @return data frame: number of kin and probability for kin of age s2

unconditional_grans_age_CYA  <- function(y, s2, U_mat, F_mat, Q){
  n <- nrow(U_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  U_prob <- U_kin_death(0, s2-1 , Q, U_mat)
  cus_list <- list()
  for(s1 in (s2 + actual_ages[1]) : (s2 + actual_ages[length(actual_ages)]) ){
    age_YA_at_Cousin <- s1 - s2
    age_foc_when_Cousin_made <- y - s2

    if(age_YA_at_Cousin %in% actual_ages){
      ya_pdf <- unconditional_grans_age_YA(age_foc_when_Cousin_made, age_YA_at_Cousin, U_mat, F_mat, Q)
      ya_pdf <- ya_pdf$prob
      QQ <- Q_matrix(age_YA_at_Cousin, Q, F_mat)
      newborn_cousin <- QQ %*% ya_pdf
      suv_cousin <- U_prob %*% newborn_cousin
      cus_list[[(1+length(cus_list))]] <- suv_cousin}
  }
  prob <- convoluion_nth( length(cus_list) , cus_list )
  df <- data.frame(number = seq(0, Q-1),
                   prob = prob)
  return(df)
}


#' Title calculate the kin-number distribution of Focal's cousin descending from younger aunts
#'
#' @param y age of Focal
#' @param s2 age of cousin from older aunt, when Focal is y
#' @param U_mat survival matrix
#' @param F_mat Fertilty matrix
#' @param Q number of possible kin
#'
#' @return vector in pmf form
#'
Matrix_func_age_CYA  <- function(y, s2, U_mat, F_mat, Q){
  n <- nrow(U_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0)
  U_prob <- U_kin_death(0, s2-1 , Q, U_mat)
  cus_list <- list()
  for(s1 in (s2 + actual_ages_of_mothering[1]) : (s2 + actual_ages_of_mothering[length(actual_ages_of_mothering)]) ){
    age_YA_at_Cousin <- s1 - s2
    age_foc_when_Cousin_made <- y - s2
    if(age_YA_at_Cousin %in% actual_ages_of_mothering){
      ya_pdf <- Matrix_func_age_YA_CALC(age_foc_when_Cousin_made, age_YA_at_Cousin, U_mat, F_mat, Q)
      QQ <- Q_matrix(age_YA_at_Cousin, Q, F_mat)
      newborn_cousin <- QQ %*% ya_pdf
      suv_cousin <- U_prob %*% newborn_cousin
      cus_list[[(1+length(cus_list))]] <- suv_cousin}
  }
  prob <- convoluion_nth( length(cus_list) , cus_list )
  return(as.vector(prob))
}

#' Title calculate the kin-number distribution of Focal's cousin descending from younger aunts
#'
#' @param y age of Focal
#' @param s2 age of cousin from older aunt, when Focal is y
#' @param U_mat survival matrix
#' @param F_mat Fertilty matrix
#' @param Q number of possible kin
#'
#' @return vector in pmf form
#'
CYA_dist_quick <- function(y, s2, U_mat, F_mat, Q){

  age_pdfs <- fert_dists(F_mat, Q) ### age specific reproduction pmfs (based on Poisson dist.)
  U_suv_CYA <- U_kin_death(0, s2-1, Q, U_mat) ## matrix to project the pmf of newborn cousins to age s2
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){mothers_age(U_mat, F_mat, x)})) # pdf of mother's ages
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) # ages of reproduction under-bar/bar n
  age_combs <- expand.grid(actual_ages_of_mothering, actual_ages_of_mothering) ## possible b1/b2 combinations (repro. ages of mother and gran)
  aunts_age_at_cousin <- (actual_ages_of_mothering + y - s2) ## convolution limits for Y_aunt's reproduction
  ## Main function here:
  matrix_result_CYA <- lapply( actual_ages_of_mothering , function(YA_age_CYA){
    rho_probs <- do.call("cbind", lapply(1:nrow(age_combs),
                                         function(x){b11 <- age_combs[,1][x]; b22 <- age_combs[,2][x]
    probable_ages_of_mothering[b11]*probable_ages_of_mothering[b22]} )
                        ) # cbind a vector pmf of probable combs. of ages of mothering
    rho_probs <- rep(1, Q) %*% rho_probs ## make into matrix for Schur/Hamadarnd product with matrix of respective fert dists.
    dist_mat_Gran <- do.call("cbind", lapply(1:nrow(age_combs), function(x){b11 <- age_combs[,1][x]; b22 <- age_combs[,2][x];
    if( (b11 + y > (YA_age_CYA + s2)) & ((b11 + b22 - (YA_age_CYA + s2) + y) %in% actual_ages_of_mothering) ){
      age_gran_ya <- b11 + b22 - (YA_age_CYA + s2) + y;
      gran_vec <- e_vector(2, Q)
      gran_vec_at_ya <- U_kin_death(b22, age_gran_ya-1, Q, U_mat) %*% gran_vec
      Q_mat_Gran_YA <- Q_matrix(age_gran_ya, Q, F_mat)
      younger_aunt_at_birth <- Q_mat_Gran_YA %*% gran_vec_at_ya
      U_prob_aunt <- U_kin_death(0, YA_age_CYA - 1, Q, U_mat);
      return(U_prob_aunt %*% younger_aunt_at_birth )}
    else{return(c(1, rep(0, Q-1)))}})
                            ) # cbind a matrix with cols representing age-specific repro of gran, given b1,b2
    YA_pdf_project <- (rho_probs*dist_mat_Gran) %*% rep(1, ncol(dist_mat_Gran)) # Schur multiplication of rho mat and fert mat = grans repro pmf for YA
    Q_mat_YA_CYA <- Q_matrix(YA_age_CYA, Q, F_mat)
    newborn_CYA <- Q_mat_YA_CYA %*% YA_pdf_project
    suv_CYA <- U_suv_CYA %*% newborn_CYA
    return(suv_CYA) }  )
  prob <- convoluion_nth( length(matrix_result_CYA) , matrix_result_CYA )
  return(prob)
}

