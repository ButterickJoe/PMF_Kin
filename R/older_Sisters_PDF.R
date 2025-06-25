

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

conditional_mothers_age_OS <- function(y, s1, b1, U_mat, F_mat, Q){
  if(s1 <= y){stop("Older sister has to be strictly older")}
  age_mother_os <- b1 - s1 + y
  age_pdfs <- fert_dists(F_mat, Q)
  older_sis_at_birth <- age_pdfs[[(age_mother_os+1)]]
  U_prob <- U_kin_death(0, (s1-1), Q, U_mat)
  older_sis_after_suv <- U_prob %*% older_sis_at_birth
  return(older_sis_after_suv)
}


#' Title unconditional distribution of older sisters, weighted over all ages when mother had Focal
#'
#' @param y (at present -- y in pdf document)
#' @param s1 (at present -- s_1 in pdf document)
#' @param U_mat (matrix of static survival probs -- subdiagonal)
#' @param F_mat (matrix of static fertility rates -- first row)
#' @param Q (the number of kin categories 0, 1, 2, ..., Q-1, i.e., max Q-1 kin over life)
#'
#' @return data frame with columns: numbers 0,..,Q-1, probs of numbers,

unconditional_mothers_age_OS  <- function(y, s1, U_mat, F_mat, Q){
  n <- nrow(U_mat)
  if(y >= s1){stop("Older sister has to be strictly older")}
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  actual_ages <- actual_ages - 1
  dd <- 0
  full_df <- data.frame()
  for(age_mother_focal in actual_ages){
    age_mother_os <- age_mother_focal - s1 + y
    prob_mother_age_Focal <- probable_ages_of_mothering[(1+age_mother_focal)]
    dd <- dd + prob_mother_age_Focal
    if(age_mother_os %in% actual_ages){
      temp_df <- data.frame(number = seq(0, Q-1),
                 prob = conditional_mothers_age_OS(y, s1, age_mother_focal, U_mat, F_mat, Q),
                 age_older_sis = rep(s1, Q),
                 age_mother_older_sis = rep(age_mother_os, Q),
                 reason = "mother reproduces OS",
                 prob_mf = prob_mother_age_Focal)

    }else{
      temp_df <- data.frame(number = seq(0, Q-1),
                            prob = c(1, rep(0, Q-1)),
                            age_older_sis = rep(s1, Q),
                            age_mother_older_sis = rep(age_mother_os, Q),
                            reason = "mother not eligable for OS",
                            prob_mf = prob_mother_age_Focal)

    }
    full_df <- rbind(full_df, temp_df)
    }
    full_df <- full_df %>%
      dplyr::mutate(prob = prob*prob_mf)
    full_df <- full_df %>% dplyr::group_by(number) %>%
    dplyr::summarise(prob = sum(prob)) %>%
    dplyr::ungroup()

  return(full_df)
}


#' Title Matrix implementation of the above "unconditional_mothers_age_OS"
#'
#' @param y (at present -- y in pdf document)
#' @param s1 (at present -- s_1 in pdf document)
#' @param U_mat (matrix of static survival probs -- subdiagonal)
#' @param F_mat (matrix of static fertility rates -- first row)
#' @param Q (the number of kin categories 0, 1, 2, ..., Q-1, i.e., max Q-1 kin over life)
#'
#' @return vector with pmf of sisters aged s1 when Focal is y
#'
Matrix_func_age_OS  <- function(y, s1, U_mat, F_mat, Q){
  n <- nrow(U_mat)
  if(y >= s1){stop("Older sister has to be strictly older")}
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  actual_ages <- actual_ages - 1
  vec <- 0
  for(mothers_age_Focal in actual_ages){
    mothers_age_OS <- mothers_age_Focal - s1 + y
    pm <- probable_ages_of_mothering[(1+mothers_age_Focal)]
    if(mothers_age_OS %in% actual_ages){
      vec <- vec + pm*conditional_mothers_age_OS(y, s1, mothers_age_Focal, U_mat, F_mat, Q)
    }
    else{vec <- vec + pm*c(1, rep(0, Q-1))}
  }
  matrix_result <- vec
  return(as.vector(matrix_result))
}
