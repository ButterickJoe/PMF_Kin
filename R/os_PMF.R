

#' Title Conditional PMF for Older sisters
#'
#' @param y age of Focal at present
#' @param s1 age of Focal's older sister at present
#' @param b1 (age of mother when had focal, b_1 in pdf document)
#' @param U_mat Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal
#' @param F_mat Matrix. fertility rate of age class i in the i-th entry of top-row
#' @param Q Scalar. Life-time maximum kin-number
#'
#' @return Vector. PMF for the number-distribution of older sisters of age s1 when Focal is y, conditional on mother's age b1

os_PMF_conditional <- function(y, s1, b1, U_mat, F_mat, Q){
  if(s1 <= y){stop("Older sister has to be strictly older")}
  age_mother_os <- b1 - s1 + y
  age_pdfs <- fert_dists(F_mat, Q)
  older_sis_at_birth <- age_pdfs[[(age_mother_os+1)]]
  U_prob <- U_kin_death(0, (s1-1), Q, U_mat)
  older_sis_after_suv <- U_prob %*% older_sis_at_birth
  return(older_sis_after_suv)
}



#' Title Unconditional PMF for Older sisters
#'
#' @param y Scalar. age of Focal at present
#' @param s1 Scalar. age of Focal's older sister at present
#' @param U_mat Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal
#' @param F_mat Matrix. fertility rate of age class i in the i-th entry of top-row
#' @param Q Scalar. Life-time maximum kin-number
#'
#' @return Vector. PMF for the number-distribution of older sisters of age s1 when Focal is y
#'
os_PMF <- function(y, s1, U_mat, F_mat, Q){

  if(y >= s1){stop("Older sister has to be strictly older")}
  probable_ages_of_mothering <- sapply(1:ncol(F_mat), function(x) mothers_age(U_mat, F_mat, x))
  # Actual ages of mothering (non-zero probability)
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) - 1
  matrix_result <- 0
  for(mothers_age_Focal in actual_ages_of_mothering){
    mothers_age_OS <- mothers_age_Focal - s1 + y
    pm <- probable_ages_of_mothering[(1+mothers_age_Focal)]
    if(mothers_age_OS %in% actual_ages_of_mothering){
      matrix_result <- matrix_result + pm*os_PMF_conditional(y, s1, mothers_age_Focal, U_mat, F_mat, Q)
    }
    else{matrix_result <- matrix_result + pm*c(1, rep(0, Q-1))}
  }
  return(as.vector(matrix_result))
}
