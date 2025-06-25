

#' Title Unconditional dist of Focal's cousins via older aunts
#'
#' @param y Scalar. age of Focal at present
#' @param s2 Scalar. age of cousin at present (when Focal is y)
#' @param U_mat Matrix. probabilities of surviving from age class i to i+1 in the i-th entry of sub-diagonal
#' @param F_mat Matrix. fertility rate of age class i in the i-th entry of top-row
#' @param Q Scalar. Life-time maximum kin-number
#'
#' @return Vector. PMF for the number-distribution of cousins via older aunts of age s1 when Focal is y
#'
coa_PMF <- function(y, s2, U_mat, F_mat, Q){
  # Compute age-specific reproduction probability mass functions
  age_pdfs <- fert_dists(F_mat, Q)
  U_suv_COA <- U_kin_death(0, s2-1, Q, U_mat) ## matrix to project the pmf of newborn cousins to age s2
  probable_ages_of_mothering <- sapply(1:ncol(F_mat), function(x) mothers_age(U_mat, F_mat, x))
  # Actual ages of mothering (non-zero probability)
  actual_ages_of_mothering <- which(probable_ages_of_mothering != 0) - 1
  # All possible combinations of mother and grandmother ages in matrix form (rows repeated)
  combinations <- expand.grid(a = actual_ages_of_mothering, b = actual_ages_of_mothering)
  rho_probs <- matrix(diag(outer(probable_ages_of_mothering[(combinations$a+1)], probable_ages_of_mothering[(combinations$b+1)], "*")), nrow = 1)
  rho_probs_mat <- rep(1, Q) %*% rho_probs

  # Main computation:
  result_list <- lapply(actual_ages_of_mothering, function(OA_age_COA) {
    b1 <- combinations[, 1] # age combs for mother and grandmother
    b2 <- combinations[, 2]
    # age of gran when having aunt
    gran_age <- b1 + b2 - (OA_age_COA + s2) + y
    index_list <- which((b1 + y < (OA_age_COA + s2)) & (gran_age %in% actual_ages_of_mothering)) ## index_list = indices representing ages at which gran can have aunt, such that aunt is younger then mother
    # Conditional reproductive pmfs of gran (for element-wise product through resp. probs in "rho_probs_mat")
    dist_mat_gran <- matrix(0, nrow = Q, ncol = nrow(combinations))
    dist_mat_gran[1, ] <- 1  #
    # If condition_T is met
    if (length(index_list) > 0) {
      U_prob_aunt <- U_kin_death(0, OA_age_COA - 1, Q, U_mat) # prob gran survives from having mother, up to age when she has aunt
      for (i in index_list) {
        age_gran_oa <- gran_age[i]
        older_aunt_at_birth <- age_pdfs[[(age_gran_oa+1)]]
        dist_mat_gran[, i] <- U_prob_aunt %*% older_aunt_at_birth
      }
    }
    # Projected PMF for grandmother's reproduction
    OA_pdf_project <- (rho_probs_mat * dist_mat_gran) %*% rep(1, ncol(dist_mat_gran))
    # Distribution of newborn cousins to aunt at age "OA_age_COA"
    Q_mat_OA_COA <- Q_matrix(OA_age_COA, Q, F_mat)
    newborn_COA <- Q_mat_OA_COA %*% OA_pdf_project
    # Distribution of cousins after U_suv_CYA
    suv_COA <- U_suv_COA %*% newborn_COA
    return(suv_COA)
  })
  # Convolve distributions over all to get final probability
  prob <- convoluion_nth(length(result_list), result_list)
  return(prob)
}
