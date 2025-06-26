

#' Title Unconditional dist of Focal's cousins via younger aunts
#'
#' @param y Scalar. age of Focal at present
#' @param s2 Scalar. age of cousin at present (when Focal is y)
#' @param U_mat Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal
#' @param F_mat Matrix. fertility rate of age class i in the i-th entry of top-row
#' @param Q Scalar. Life-time maximum kin-number
#'
#' @return Vector. PMF for the number-distribution of cousins via younger aunts of age s1 when Focal is y
#'
cya_PMF <- function(y, s2, U_mat, F_mat, Q) {

  gran_vec <- e_vector(2, Q)
  # Matrix projecting PMF of newborn cousins to age s2
  U_suv_CYA <- U_kin_death(0, s2 - 1, Q, U_mat)
  # Probabilities of mother ages
  probable_ages_of_mothering <- sapply(1:ncol(F_mat), function(x) mothers_age(U_mat, F_mat, x))
  # Actual ages of mothering (non-zero probability)
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) - 1
  ## make a list of possible age-specific pmf reproduction matrices for to save computing each time in the lapply below
  Q_dists <- lapply(actual_ages_of_mothering, function(x){Q_matrix(x, Q, F_mat)})
  # All possible combinations of mother and grandmother ages in matrix form (rows repeated)
  combinations <- expand.grid(a = actual_ages_of_mothering, b = actual_ages_of_mothering)
  # Respecitive probabilities that mother (b1) and gran (b2) had Focal's direct ancestor at these ages
  rho_probs <- matrix(diag(outer(probable_ages_of_mothering[(combinations$a+1)], probable_ages_of_mothering[(combinations$b+1)], "*")), nrow = 1)
  ## in matrix form
  rho_probs_mat <- rep(1, Q) %*% rho_probs

  # list of potential youger cousins -- the convolution of which is the result:
  result_list <- lapply(actual_ages_of_mothering, function(YA_age_CYA) {
    b1 <- combinations[, 1] # age of mother when she had Focal
    b2 <- combinations[, 2] # age of gran when she had mother
    gran_age <- b1 + b2 - (YA_age_CYA + s2) + y # age of gran when having aunt
    index_list <- which((b1 + y > (YA_age_CYA + s2)) & (gran_age %in% actual_ages_of_mothering)) ## ages of gran at which she can have aunt, such that aunt is younger then mother
    # Conditional reproductive pmfs of gran (for element-wise product through resp. probs in "rho_probs_mat")
    dist_mat_gran <- matrix(0, nrow = Q, ncol = nrow(combinations))
    dist_mat_gran[1, ] <- 1  # by default we assume no reproduction of newborns unless otherwise
    # If condition_T is met
    if (length(index_list) > 0) { # if otherwise
      U_prob_aunt <- U_kin_death(0, YA_age_CYA - 1, Q, U_mat) # prob gran survives from having mother to age when she has aunt
      for (i in index_list) {
        age_gran_ya <- gran_age[i] # gran's age at younger aunt
        gran_at_ya <- U_kin_death(b2[i], age_gran_ya - 1, Q, U_mat) %*% gran_vec # pmf of gran (Bernoulli) when having aunt
        Q_gran <- Q_dists[[(age_gran_ya+1-actual_ages_of_mothering[1])]] ## reproduction pmf at this age
        younger_aunt <- Q_gran %*% gran_at_ya ## pmf of newborn aunts
        dist_mat_gran[, i] <- U_prob_aunt %*% younger_aunt ## pmf of these aunts when they have cousin
      }
    }
    # Projected PMF for grandmother's reproduction
    YA_pdf_project <- (rho_probs_mat * dist_mat_gran) %*% rep(1, ncol(dist_mat_gran)) ## pmf of aunts at cousin using the probabilistic sum over b1, b2
    Q_YA <- Q_dists[[(YA_age_CYA+1-actual_ages_of_mothering[1])]] # reproduction of aunt
    newborn_CYA <- Q_YA %*% YA_pdf_project ## --> newborn cousins
    suv_CYA <- U_suv_CYA %*% newborn_CYA ## Cousins who suvive to age s2 when Focal is y
    return(suv_CYA)
  })
  # Convolution of distributions over all ages of aunt s1 in actual_ages_of_mothering, under the constraint b1 + y > (YA_age_CYA + s2)
  prob <- convoluion_nth(length(result_list), result_list)
  return(prob)
}
