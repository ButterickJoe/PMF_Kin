



#' Title
#'
#' @param y age of Focal
#' @param s2 age of cousin from older aunt, when Focal is y
#' @param U_mat survival matrix
#' @param F_mat Fertilty matrix
#' @param Q number of possible kin
#'
#' @return data frame: number of kin with corresponding probabilities, for kin at age s2

unconditional_grans_age_COA  <- function(y, s2, U_mat, F_mat, Q){
  n <- nrow(U_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  U_prob <- U_kin_death(0, s2-1 , Q, U_mat)
  cus_list <- list()
  for(s1 in (s2 + actual_ages[1]) : (s2 + actual_ages[length(actual_ages)]) ){
    age_OA_at_Cousin <- s1 - s2
    age_foc_when_Cousin_made <- y - s2
    if(age_OA_at_Cousin %in% actual_ages){
    oa_pdf <- uncondition_grans_age_OA(age_foc_when_Cousin_made, age_OA_at_Cousin, U_mat, F_mat, Q)
    oa_pdf <- oa_pdf$prob
    #oa_pdf <- Matrix_func_age_OA(age_foc_when_Cousin_made, age_OA_at_Cousin, U_mat, F_mat, Q)
    QQ <- Q_matrix(age_OA_at_Cousin, Q, F_mat)
    newborn_cousin <- QQ %*% oa_pdf
    suv_cousin <- U_prob %*% newborn_cousin
    cus_list[[(1+length(cus_list))]] <- suv_cousin}
  }
  prob <- convoluion_nth( length(cus_list) , cus_list )
  df <- data.frame(number = seq(0, Q-1),
                   prob = prob)
  return(df)

}


#' Title
#'
#' @param y age of Focal
#' @param s2 age of cousin from older aunt, when Focal is y
#' @param U_mat survival matrix
#' @param F_mat Fertilty matrix
#' @param Q number of possible kin
#'
#' @return data frame: number of kin with corresponding probabilities, for kin at age s2

Matrix_func_age_COA  <- function(y, s2, U_mat, F_mat, Q){
  n <- nrow(U_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  U_prob <- U_kin_death(0, s2-1 , Q, U_mat)
  cus_list <- list()
  for(s1 in (s2 + actual_ages[1]) : (s2 + actual_ages[length(actual_ages)]) ){
    age_OA_at_Cousin <- s1 - s2
    age_foc_when_Cousin_made <- y - s2
    if(age_OA_at_Cousin %in% actual_ages){
      oa_pdf <- Matrix_func_age_OA_CALC(age_foc_when_Cousin_made, age_OA_at_Cousin, U_mat, F_mat, Q)
      QQ <- Q_matrix(age_OA_at_Cousin, Q, F_mat)
      newborn_cousin <- QQ %*% oa_pdf
      suv_cousin <- U_prob %*% newborn_cousin
      cus_list[[(1+length(cus_list))]] <- suv_cousin}
  }
  prob <- convoluion_nth( length(cus_list) , cus_list )
  return(as.vector(prob))
}


### Fastest way so far

COA_dist_quick <- function(y, s2, U_mat, F_mat, Q){

  age_pdfs <- fert_dists(F_mat, Q) ### age specific reproduction pmfs (based on Poisson dist.)
  U_suv_COA <- U_kin_death(0, s2-1, Q, U_mat) ## matrix to project the pmf of newborn cousins to age s2
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) # ages of reproduction under-bar/bar n
  age_combs <- expand.grid(actual_ages_of_mothering, actual_ages_of_mothering) ## possible b1/b2 combinations (repro. ages of mother and gran)
  aunts_age_at_cousin <- (actual_ages_of_mothering + y + 1 - s2) ## convolution limits for O_aunt's reproduction
  ## Main function here:
  matrix_result_COA <- lapply(actual_ages_of_mothering, function(OA_age_COA){
    rho_probs <- do.call("cbind", lapply(1:nrow(age_combs), function(x){b11 <- age_combs[,1][x]; b22 <- age_combs[,2][x]
    probable_ages_of_mothering[b11]*probable_ages_of_mothering[b22]} )
                        ) # cbind a vector pmf of probable combs. of ages of mothering
    rho_probs <- rep(1, Q) %*% rho_probs ## make into matrix for Schur/Hamadarnd product with matrix of respective fert dists.
    dist_mat_Gran <- do.call("cbind", lapply(1:nrow(age_combs), function(x){b11 <- age_combs[,1][x]; b22 <- age_combs[,2][x];
    if( (b11 + y < (OA_age_COA + s2)) & ((b11 + b22 - (OA_age_COA + s2) + y) %in% actual_ages_of_mothering) ){
      age_gran_oa <- b11 + b22 - (OA_age_COA + s2) + y;
      older_aunt_at_birth <- age_pdfs[[age_gran_oa]];
      U_prob_aunt <- U_kin_death(0, OA_age_COA - 1, Q, U_mat);
      return(U_prob_aunt %*% older_aunt_at_birth )}else{return(c(1, rep(0, Q-1)))}})
                        ) # cbind a matrix with cols representing age-specific repro of gran, given b1,b2
    OA_pdf_project <- (rho_probs*dist_mat_Gran) %*% rep(1, ncol(dist_mat_Gran))
    Q_mat_OA_COA <- Q_matrix(OA_age_COA, Q, F_mat)
    newborn_COA <- Q_mat_OA_COA %*% OA_pdf_project
    suv_COA <- U_suv_COA %*% newborn_COA
    return(suv_COA) }  )
  prob <- convoluion_nth(length(matrix_result_COA), matrix_result_COA)
  return(prob)
}
