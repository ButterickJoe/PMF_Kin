


#' Title conditional distribution of older sisters, assuming we know mother's age at having Focal
#'
#' @param y (at present -- y in pdf document)
#' @param s1 (at present -- s_1 in pdf document)
#' @param b1 (age of mother when had focal, b_1 in pdf document)
#' @param U_mat (matrix of static survival probs -- subdiagonal)
#' @param F_mat (matrix of static fertility rates -- first row)
#' @param Q (the number of kin categories 0, 1, 2, ..., Q-1, i.e., max Q-1 kin over life)
#'
#' @return the pdf of older sisters of Focal age s_1, when focal is aged y, conditioned on mother being age b_1 at focal

unconditional_Focals_age_Dau <- function(y, s1, U_mat, F_mat, Q){

  age_pdfs <- fert_dists(F_mat, Q)
  daughter_at_birth <- age_pdfs[[(y-s1)]]
  U_prob <- U_kin_death(0, (s1-1), Q, U_mat)
  daughter_after_suv <- U_prob %*% daughter_at_birth
  return(daughter_after_suv)
}
