

#' Title Unconditional PDF for nieces through Focal's older sisters
#' @param y age of Focal
#' @param s2 age of older niece when Focal is y
#' @param U_mat mort mat
#' @param F_mat fert mat
#' @param Q kin-numbers
#'
#' @return data frame with numbers and probs of kin aged s2
nos_PMF  <- function(y, s2, U_mat, F_mat, Q){

  probable_ages_of_mothering <- sapply(1:ncol(F_mat), function(x) mothers_age(U_mat, F_mat, x))
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) - 1
  U_prob_niece <- U_kin_death(0, s2-1, Q, U_mat)
  age_dists <- fert_dists(F_mat, Q)
  Q_dists <- lapply(actual_ages_of_mothering, function(x){Q_matrix(x, Q, F_mat)})
  rho_probs <- probable_ages_of_mothering[(actual_ages_of_mothering+1)]
  rho_probs_mat <- rep(1, Q) %*% matrix(rho_probs, nrow = 1)

  result_list <- lapply(actual_ages_of_mothering, function(OS_AT_NOS) {
    b1 <- actual_ages_of_mothering
    mum_age_at_os <- b1 + y - (OS_AT_NOS + s2) # age of mum when having older sis
    index_list <- which((y < OS_AT_NOS + s2) & (mum_age_at_os %in% actual_ages_of_mothering)) ## ages of gran at which she can have aunt, such that aunt is younger then mother
    # Conditional reproductive pmfs of mum (for element-wise product through resp. probs in "rho_probs_mat")
    dist_mat_sis <- matrix(0, nrow = Q, ncol = length(actual_ages_of_mothering))
    dist_mat_sis[1, ] <- 1  # by default we assume no reproduction of newborns unless otherwise
    # If condition_T is met
    if (length(index_list) > 0) { # if otherwise
      U_prob_sis <- U_kin_death(0, OS_AT_NOS - 1, Q, U_mat) # prob sis survives from birth to age when she has niece
      for (i in index_list) {
        mum_age_os <- mum_age_at_os[i]
        older_sis_newborn <- age_dists[[(1+mum_age_os)]]
        dist_mat_sis[, i] <- U_prob_sis %*% older_sis_newborn ## pmf of these aunts when they have cousin
      }
    }
    # Projected PMF for grandmother's reproduction
    OS_pfm_at_niece <- (rho_probs_mat * dist_mat_sis) %*% rep(1, ncol(dist_mat_sis)) ## pmf of aunts at cousin using the probabilistic sum over b1, b2
    Q_OS <- Q_dists[[(OS_AT_NOS+1-actual_ages_of_mothering[1])]] # reproduction of aunt
    newborn_NOS <- Q_OS %*% OS_pfm_at_niece ## --> newborn cousins
    suv_NOS <- U_prob_niece %*% newborn_NOS ## Cousins who suvive to age s2 when Focal is y
    return(suv_NOS)
  })
  # Convolution of distributions over all ages of aunt s1 in actual_ages_of_mothering, under the constraint b1 + y > (YA_age_CYA + s2)
  prob <- convoluion_nth(length(result_list), result_list)
  return(prob)
}
