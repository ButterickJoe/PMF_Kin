




#' Title Conditional dist of Focal's older aunts: conditioned on mother age b1 and gran age b2
#'
#' @param y age of Focal at present
#' @param s1 age of older aunt at present
#' @param b1 age of mother when had Focal
#' @param b2 age of gran when had mother
#' @param U_mat Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal
#' @param F_mat Matrix. fertility rate of age class i in the i-th entry of top-row
#' @param Q Scalar. Life-time maximum kin-number
#'
#' @return Vector. PMF for the number-distribution of older aunts of age s1 when Focal is y, conditional on both mother's age b1 and gran's age b2

oa_PMF_conditional <- function(y, s1, b1, b2, U_mat, F_mat, Q){
  if(s1 <= y + b1){stop("Older aunts have to be strictly older than mother")}
  age_gran_oa <- b1 + b2 - s1 + y
  age_pdfs <- fert_dists(F_mat, Q)
  older_aunt_at_birth <- age_pdfs[[age_gran_oa]]
  U_prob <- U_kin_death(0, s1-1, Q, U_mat)
  older_aunts_after_suv <- U_prob %*% older_aunt_at_birth
  return(older_aunts_after_suv)
}



#' Title Unconditional PMF for older aunts
#'
#' @param y age of Focal at present
#' @param s1 age of older aunt at present
#' @param U_mat Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal
#' @param F_mat Matrix. fertility rate of age class i in the i-th entry of top-row
#' @param Q Scalar. Life-time maximum kin-number
#'
#' @return Vector. PMF for the number-distribution of older aunts of age s1 when Focal is y
#'
oa_PMF  <- function(y, s1, U_mat, F_mat, Q){
  U_prob_aunt <- U_kin_death(0, s1-1, Q, U_mat)
  # Compute age-specific reproduction probability mass functions
  age_pdfs <- fert_dists(F_mat, Q)
  # Probabilities of mother ages
  probable_ages_of_mothering <- sapply(1:ncol(F_mat), function(x) mothers_age(U_mat, F_mat, x))
  # Actual ages of mothering (non-zero probability)
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) - 1
  combinations <- expand.grid(a = actual_ages_of_mothering, b = actual_ages_of_mothering)
  rho_probs <- matrix(diag(outer(probable_ages_of_mothering[(combinations$a+1)], probable_ages_of_mothering[(combinations$b+1)], "*")), nrow = 1)
  rho_probs_mat <- rep(1, Q) %*% rho_probs ## gives the bi-variate pmf of mother/grans ages of reproduction
  ### Conditional on each probability from rho_probs_mat we now derive grans reproduction of aunts -- we then multiply these element-wise
  b1 <- combinations[, 1] # age of mother
  b2 <- combinations[, 2] # age of grandmother
  # age range of gran when having aunt
  age_gran_range <- b1 + b2 - s1 + y
  index_list <- which((b1 + y < s1) & (age_gran_range %in% actual_ages_of_mothering)) ## filter to condition on aunt younger than mom
  ## matrix with columns ages of gran, and rows her pmf of newborn younger aunts
  dist_mat_gran <- matrix(0, nrow = Q, ncol = nrow(combinations))
  dist_mat_gran[1, ] <- 1  # so that unless filled the rows are pmfs of zero newborns (no aunts yet!)
  if (length(index_list) > 0) {

    for (i in index_list) {
      age_gran_oa <- age_gran_range[i]
      older_aunt_at_birth <- age_pdfs[[(age_gran_oa+1)]]
      dist_mat_gran[, i] <- U_prob_aunt %*% older_aunt_at_birth
    }
  }
  matrix_result <- (rho_probs_mat*dist_mat_gran) %*% rep(1, ncol(dist_mat_gran))
  return(matrix_result)
}


