

#' Title Conditional (on mother's age at Focal) pdf for Younger sisters
#'
#' @param y age of Focal (at present)
#' @param s1 age of younger sister (at present)
#' @param b1 age of mother when having Focal
#' @param U_mat matrix of survival probabilities on sub-diagonal
#' @param F_mat matrix of age-specific fertilities on top-row
#' @param Q number of kin-number classes
#'
#' @return vector pdf of number / probs of younger sisters of age s1

conditional_mothers_age_YS_alt <- function(y, s1, b1, U_mat, F_mat, Q){
  #if(s1 >= y){stop("Younger sister has to be strictly younger")}
  age_mother_ys <- b1 - s1 + y
  mother_vec <- e_vector(2, Q)
  mother_vec_at_ys <- U_kin_death(b1+1, age_mother_ys-1, Q, U_mat) %*% mother_vec
  Q_mat <- Q_matrix_joe((age_mother_ys-1), Q, F_mat)
  younger_sis_at_birth <- Q_mat %*% mother_vec_at_ys
  U_prob <- U_kin_death(0, s1-1 , Q, U_mat)
  younger_sis_after_suv <- U_prob %*% younger_sis_at_birth
  return(younger_sis_after_suv)
}


#' Title
#'
#' @param y age of Focal (at present)
#' @param s1 age of younger sister (at present)
#' @param U_mat matrix of survival probabilities on sub-diagonal
#' @param F_mat matrix of age-specific fertilities on top-row
#' @param Q number of kin-number classes
#' ## Note that we sum over probable ages of mother in this case
#' @return data.fram of pdf of number / probs of younger sisters of age s1

unconditional_mothers_age_YS_alt  <- function(y, s1, U_mat, F_mat, Q, m_prob, m_age){
  alt_rho <- readRDS("data/alt_rhob1")
  rhob1 <- m_prob
  age_b1 <- m_age
  n <- nrow(U_mat)
  #if(y <= s1){stop("Younger sister has to be strictly younger")}
  probable_ages_of_mothering <- rhob1
  actual_ages <- age_b1
  dd <- 0
  full_df <- data.frame()
  for(age_mother_focal in actual_ages){
    age_mother_ys <- age_mother_focal - s1 + y
    prob_mother_age_Focal <- probable_ages_of_mothering[(age_mother_focal-min(actual_ages)+1)] ## because f_mat[1,1] corresponds to age 0 not 1
    dd <- dd + prob_mother_age_Focal
    if(age_mother_ys %in% actual_ages){
      temp_df <- data.frame(number = seq(0, Q-1),
                            prob = conditional_mothers_age_YS_alt(y, s1, age_mother_focal, U_mat, F_mat, Q),
                            age_younger_sis = rep(s1, Q),
                            age_mother_younger_sis = rep(age_mother_ys, Q),
                            reason = "mother reproduces YS",
                            prob_mf = prob_mother_age_Focal)

    }else{
      temp_df <- data.frame(number = seq(0, Q-1),
                            prob = c(1, rep(0, Q-1)),
                            age_younger_sis = rep(s1, Q),
                            age_mother_younger_sis = rep(age_mother_ys, Q),
                            reason = "mother not eligable for YS",
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


#' Title Matrix implementation of the above "unconditional_mothers_age_YS"
#'
#' @param y age of Focal (at present)
#' @param s1 age of younger sister (at present)
#' @param U_mat matrix of survival probabilities on sub-diagonal
#' @param F_mat matrix of age-specific fertilities on top-row
#' @param Q number of kin-number classes
#' ## Note that we sum over probable ages of mother in this case
#' @return vector of pdf of number / probs of younger sisters of age s1
#'
Matrix_func_age_YS_alt  <- function(y, s1, U_mat, F_mat, Q, m_prob, m_age){
  n <- nrow(U_mat)

  #if(y <= s1){stop("Younger sister has to be strictly younger")}
  probable_ages_of_mothering <- m_prob
  actual_ages <- m_age
  vec <- 0
  dd <- 0
  for(mothers_age_Focal in actual_ages){

    mothers_age_YS <- mothers_age_Focal - s1 + y
    pm <- probable_ages_of_mothering[(mothers_age_Focal-min(actual_ages)+1)]
    dd <- dd + pm

    if(mothers_age_YS %in% actual_ages){
      vec <- vec + pm*conditional_mothers_age_YS_alt(y, s1, mothers_age_Focal, U_mat, F_mat, Q)
    }
    else{vec <- vec + pm*c(1, rep(0, Q-1))}
    }

  matrix_result <- vec
  return(list(as.vector(matrix_result),dd))
}
