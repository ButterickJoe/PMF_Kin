



#' Title Conditional dist of Focal's younger aunts: conditioned on mother age b1 and gran age b2
#'
#' @param y age of Focal at present
#' @param s1 age of younger aunt at present
#' @param b1 age of mother when had Focal
#' @param b2 age of gran when had mother
#' @param U_mat matrix of mortality probs
#' @param F_mat matrix of fertility rates
#' @param Q number of possible kin
#'
#' @return vector of probs of number of kin 0,1,2,..,Q
#'
conditional_grans_age_YA <- function(y, s1, b1, b2, U_mat, F_mat, Q){
  if(s1 >= y + b1){stop("Younger aunts have to be strictly younger than mother")}
  age_gran_ya <- b1 + b2 - s1 + y
  gran_vec <- e_vector(2, Q)
  gran_vec_at_ya <- U_kin_death(b2, age_gran_ya-1, Q, U_mat) %*% gran_vec
  Q_mat <- Q_matrix(age_gran_ya, Q, F_mat)
  younger_aunt_at_birth <- Q_mat %*% gran_vec_at_ya
  U_prob <- U_kin_death(0, s1-1, Q, U_mat)
  younger_aunts_after_suv <- U_prob %*% younger_aunt_at_birth
  return(younger_aunts_after_suv)
}



#' Title unconditional dist of Focal's younger aunts:
#'
#' @param y age of Focal
#' @param s1 age of younger aunt (both at present age FOcal)
#' @param U_mat mort mat
#' @param F_mat fert mat
#' @param Q kin-numbers
#'
#' @return data frame: numbers and probs of kin at age s1
#'
unconditional_grans_age_YA <- function(y, s1, U_mat, F_mat, Q){
  n <- nrow(F_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  aunt_df <- data.frame()
  dd <- 0
  for(b1 in actual_ages){
    rho_b1 <- probable_ages_of_mothering[b1]
    for(b2 in actual_ages){
      rho_b2 <- probable_ages_of_mothering[b2]
      prob_gran_produce_YA <- rho_b1*rho_b2
      dd <- dd + prob_gran_produce_YA
      age_gran_YA <- b1 + b2 + y - s1
      if(age_gran_YA %in% actual_ages & s1 < b1 + y){
        temp_df <- data.frame(number = seq(0, Q-1),
                              prob = conditional_grans_age_YA(y, s1, b1, b2, U_mat, F_mat, Q),
                              age_younger_aunt = rep(s1, Q),
                              age_gran_younger_aunt = rep(age_gran_YA, Q),
                              reason = "gran reproduces YA",
                              prob_mf = prob_gran_produce_YA)
      }
      else{temp_df <- data.frame(number = seq(0, Q-1),
                                 prob = c(1, rep(0, Q-1)),
                                 age_younger_aunt = rep(s1, Q),
                                 age_gran_younger_aunt = rep(age_gran_YA, Q),
                                 reason = "gran reproduces YA",
                                 prob_mf = prob_gran_produce_YA)}
      aunt_df <- rbind(aunt_df, temp_df)
    }
  }
  aunt_df <- aunt_df %>%
    dplyr::mutate(prob = prob*prob_mf)
  aunt_df <- aunt_df %>% dplyr::group_by(number) %>%
    dplyr::summarise(prob = sum(prob)) %>%
    dplyr::ungroup()

  return(aunt_df)
}


#' Title matrix form for the kin-number distribution of Focal's younger aunts
#'
#' @param y age of Focal at present
#' @param s1 age of younger aunt at present
#' @param U_mat U_mat matrix of mortality probs
#' @param F_mat matrix of fertility rates
#' @param Q kin number classes
#'
#' @return vector with pmf of sisters aged s1 when Focal is y
#'
Matrix_func_age_YA  <- function(y, s1, U_mat, F_mat, Q){
  n <- nrow(F_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  vec <- 0
  dd <- 0
  for(b1 in actual_ages){
    rho_b1 <- probable_ages_of_mothering[b1]
    for(b2 in actual_ages){
      rho_b2 <- probable_ages_of_mothering[b2]
      prob_gran_produce_YA <- rho_b1*rho_b2
      dd <- dd + prob_gran_produce_YA
      age_gran_YA <- b1 + b2 + y - s1
      if(age_gran_YA %in% actual_ages & s1 < b1 + y){
        vec <- vec + prob_gran_produce_YA*conditional_grans_age_YA(y, s1, b1, b2, U_mat, F_mat, Q)
      }
      else{vec <- vec + prob_gran_produce_YA*c(1, rep(0, Q-1))}
    }
  }
  matrix_result <- vec
  return(as.vector(matrix_result))
}

#' Title matrix form the kin-number distribution of Focal's younger aunts
#'
#' @param y age of Focal at present
#' @param s1 age of younger aunt at present
#' @param U_mat U_mat matrix of mortality probs
#' @param F_mat matrix of fertility rates
#' @param Q kin number classes
#'
#' @return vector in pmf form
#'
Matrix_func_age_YA_CALC  <- function(y, s1, U_mat, F_mat, Q){
  n <- nrow(F_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  age_combs <- expand.grid(actual_ages,actual_ages)
  rho_probs <- do.call("cbind", lapply(1:nrow(age_combs), function(x){b11 <- age_combs[,1][x];b22 <- age_combs[,2][x]
  probable_ages_of_mothering[b11]*probable_ages_of_mothering[b22]})
  )
  rho_probs <- rep(1, Q) %*% rho_probs
  dist_mat <- do.call("cbind", lapply(1:nrow(age_combs),function(x){b11 <- age_combs[,1][x];b22 <- age_combs[,2][x];
  if( (b11 + y > s1) & ((b11 + b22 + y - s1) %in% actual_ages) ){
    conditional_grans_age_YA(y, s1, b11, b22, U_mat, F_mat, Q)}
  else{c(1, rep(0, Q-1))} })
  )
  matrix_result <- (rho_probs*dist_mat) %*% rep(1, ncol(dist_mat))
  return(as.vector(matrix_result))
}

