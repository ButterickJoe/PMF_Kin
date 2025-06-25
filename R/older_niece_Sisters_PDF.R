
#' Title Unconditional PDF for nieces through Focal's older sisters
#' @param y age of Focal
#' @param s2 age of older niece when Focal is y
#' @param U_mat mort mat
#' @param F_mat fert mat
#' @param Q kin-numbers
#'
#' @return data frame with numbers and probs of kin aged s2
unconditional_mothers_age_NOS  <- function(y, s2, U_mat, F_mat, Q){

  n <- nrow(U_mat)
  probable_ages_of_mothering <- unlist(lapply(1:ncol(F_mat), function(x){ mothers_age(U_mat, F_mat, x)}))
  actual_ages <- which(probable_ages_of_mothering != 0)
  U_prob <- U_kin_death(0, s2-1 , Q, U_mat)

  niece_list <- list()
  for(s1 in (y + 1) : n){
    age_OS_at_Niece <- s1 - s2
    age_foc_when_Niece_made <- y - s2
    if(age_OS_at_Niece %in% actual_ages){
      os_pdf <- unconditional_mothers_age_OS(age_foc_when_Niece_made, age_OS_at_Niece, U_mat, F_mat, Q)
      os_pdf <- os_pdf$prob
      QQ <- Q_matrix(age_OS_at_Niece, Q, F_mat)
      newborn_niece <- QQ %*% os_pdf
      nieces_f_ <- U_prob %*% newborn_niece
      niece_list[[(1+length(niece_list))]] <- nieces_f_}}
  if(length(niece_list)>0){prob <- convoluion_nth(length(niece_list), niece_list)}
  else{prob <- c(1, seq(0,Q-1))}
  df <- data.frame(number = rep(0, Q-1),
                   prob = prob)
  return(df)

}
