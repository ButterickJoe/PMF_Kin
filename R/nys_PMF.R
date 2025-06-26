




#' Title Unconditional PDF for nieces through Focal's younger sisters
#' @param y age of Focal
#' @param s2 age of younger niece when Focal is y
#' @param U_mat mort mat
#' @param F_mat fert mat
#' @param Q kin-numbers
#'
#' @return data frame with numbers and probs of kin aged s2

nys_PMF  <- function(y, s2, U_mat, F_mat, Q){

  mother_vec <- e_vector(2, Q)
  probable_ages_of_mothering <- sapply(1:ncol(F_mat), function(x) mothers_age(U_mat, F_mat, x))
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) - 1
  U_prob_niece <- U_kin_death(0, s2-1, Q, U_mat)
  Q_dists <- lapply(actual_ages_of_mothering, function(x){Q_matrix(x, Q, F_mat)})
  rho_probs <- probable_ages_of_mothering[(actual_ages_of_mothering+1)]
  rho_probs_mat <- rep(1, Q) %*% matrix(rho_probs, nrow = 1)

  result_list <- lapply(actual_ages_of_mothering, function(YS_AT_NYS) {
    b1 <- actual_ages_of_mothering
    mum_age_at_ys <- b1 + y - (YS_AT_NYS + s2) # age of mum when having younger sis
    index_list <- which((y > YS_AT_NYS + s2) & (mum_age_at_ys %in% actual_ages_of_mothering)) ## ages of gran at which she can have aunt, such that aunt is younger then mother
    # Conditional reproductive pmfs of gran (for element-wise product through resp. probs in "rho_probs_mat")
    dist_mat_sis <- matrix(0, nrow = Q, ncol = length(actual_ages_of_mothering))
    dist_mat_sis[1, ] <- 1  # by default we assume no reproduction of newborns unless otherwise
    # If condition_T is met
    if (length(index_list) > 0) { # if otherwise
      U_prob_sis <- U_kin_death(0, YS_AT_NYS - 1, Q, U_mat) # prob sis survives from birth to age when she has niece
      for (i in index_list) {
        mum_age_ys <- mum_age_at_ys[i]
        U_prob_mum <- U_kin_death(b1[i], mum_age_ys - 1, Q, U_mat)
        pmf_mum_ys <- U_prob_mum %*% mother_vec # pmf of gran (Bernoulli) when having aunt
        Q_mum <- Q_dists[[(mum_age_ys+1-actual_ages_of_mothering[1])]] ## reproduction pmf at this age
        younger_sis_newborn <- Q_mum %*% pmf_mum_ys ## pmf of newborn aunts
        dist_mat_sis[, i] <- U_prob_sis %*% younger_sis_newborn ## pmf of these aunts when they have cousin
      }
    }
    # Projected PMF for grandmother's reproduction
    YS_pfm_at_niece <- (rho_probs_mat * dist_mat_sis) %*% rep(1, ncol(dist_mat_sis)) ## pmf of aunts at cousin using the probabilistic sum over b1, b2
    Q_YS <- Q_dists[[(YS_AT_NYS+1-actual_ages_of_mothering[1])]] # reproduction of aunt
    newborn_NYS <- Q_YS %*% YS_pfm_at_niece ## --> newborn cousins
    suv_NYS <- U_prob_niece %*% newborn_NYS ## Cousins who suvive to age s2 when Focal is y
    return(suv_NYS)
  })
  # Convolution of distributions over all ages of aunt s1 in actual_ages_of_mothering, under the constraint b1 + y > (YA_age_CYA + s2)
  prob <- convoluion_nth(length(result_list), result_list)
  return(prob)
}
