

### A = U + F (survivals + fertilities)

### Mixing distribution : for initial condition of mothers dist

pi_mix <- function(U, F){
  A <- F + U
  stable_rep_vec <- SD(A)
  pi <- F[1,]*stable_rep_vec / abs(sum( F[1,]*stable_rep_vec) )
  return(pi)
}


### Daughters
daughter_dist <- function(U, F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  IC <- rep(0, om ) ## define initial condition
  X[,1] <- IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde %*% X[,i] + F_tilde %*% e_vector(i, om)
  }
  return(X)
}
grand_daughter_dist <- function(U, F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  IC <- rep(0, n ) ## define initial condition
  X[,1] <- IC
  subsidy_vector <- daughter_dist(U, F)
  for(i in 1 : (n-1)){
    X[,i+1] <- U_tilde %*% X[,i] + F_tilde %*% subsidy_vector[,i]
  }
  return(X)
}
## Mothers
mother_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  for( i in 1 : n){
    IC <- IC + e_vector(i,om)*pi[i]
  }
  X[,1] <- U_tilde %*% IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde %*% X[,i]
  }
  return(X)
}

grand_mother_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  IC1 <- mother_dist(U,F)
  for( i in 1 : (n-1)){
    IC <- IC + pi[i]*IC1[,i]
  }
  X[,1] <-  U_tilde %*% IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde %*% X[,i]
  }
  return(X)
}
## Sisters...
older_sis_dist <- function(U, F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U, F)
  IC <- rep(0, om)
  IC1 <- daughter_dist(U, F)
  for( i in 1 : n){
    IC <- IC + pi[i]*IC1[,i]
  }
  X[,1] <- U_tilde %*% IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde%*%X[,i]
  }
  return(X)
}

younger_sis_dist <- function(U, F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  IC <- rep(0, om )
  X[,1] <- IC
  subsidy_vector <- mother_dist(U, F)
  for(i in 1 : (n - 1) ){
    X[,i+1] <- U_tilde %*% X[,i] + F_tilde %*% subsidy_vector[,i]
  }
  return(X)
}

# nieces
older_nieces_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  IC1 <- grand_daughter_dist(U,F)
  for(i in 1 : n){
    IC <- IC + pi[i]*IC1[,i]
  }
  X[,1] <- IC
  subsidy_vector <- older_sis_dist(U, F)
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde%*%X[,i] + F_tilde%*%subsidy_vector[,i]
  }
  return(X)
}

younger_nieces_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  X[,1] <- IC
  subsidy_vector <- younger_sis_dist(U, F)
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde%*%X[,i] + F_tilde%*%subsidy_vector[,i]
  }
  return(X)
}

## Aunts
# Aunts older than mother
older_aunts_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  IC1 <- older_sis_dist(U,F)
  for(i in 1 : n){
    IC <- IC + pi[i]*IC1[,i]
  }
  X[,1] <- U_tilde %*% IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde %*% X[,i]
  }
  return(X)
}

# Aunts younger than mother

younger_aunts_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  IC1 <-  younger_sis_dist(U,F)
  for(i in 1 : (n-1)){
    IC <- IC + pi[(i)]*( IC1[,(i+1)])
  }
  subsidy_vec <-  grand_mother_dist(U,F)
  X[,1] <-  IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde %*% X[,i] + F_tilde %*% subsidy_vec[,i]
  }
  return(X)
}

## Cousins

##
older_cousins_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  IC <- rep(0, om )
  IC1 <- older_nieces_dist(U,F)
  for(i in 1 : (n-1)){
    IC <- IC + pi[i]*IC1[,(i+1)]
  }
  subsidy_vec <- older_aunts_dist(U,F)
  X[,1] <- IC
  for(i in 1 : (n -1) ){
    X[,i+1] <- U_tilde %*% X[,i] + F_tilde %*% U_tilde %*% subsidy_vec[,i]
  }
  return(X)
}

##
younger_cousins_dist <- function(U,F){
  dims <- dim(U)
  n <- dims[1]
  m <- dims[2]
  U_tilde <- U
  F_tilde <- F
  om <- length(U_tilde[1,])
  X  <- matrix(0, nrow = om, ncol = n, byrow = TRUE)
  pi <- pi_mix(U,F)
  pi <- U %*% pi
  IC <- rep(0, om )
  IC1 <- younger_nieces_dist(U,F)
  for(i in 1 : (n-1)){
    IC <- IC + pi[(i)]*IC1[,(i+1)]
  }
  subsidy_vec <- younger_aunts_dist(U,F)
  X[,1] <-  IC
  for(i in 1 : (n-1) ){
    X[,i+1] <- U_tilde %*% X[,i] + F_tilde %*% U_tilde %*% subsidy_vec[,i]
  }
  return(X)
}







######################### Create df

create_deaths_and_cumsum_df <- function(U, F, list_dists, years, start_year){

  df_year_list <- list()
  for(j in years){
    ii <- as.numeric(j) - start_year + 1
    ## a,b,d,g,m,n
    daught <- daughter_dist(U[[ii]],F[[ii]])
    grand_daught <- grand_daughter_dist(U[[ii]],F[[ii]])
    moth <- mother_dist(U[[ii]],F[[ii]])
    grand_moth <- grand_mother_dist(U[[ii]],F[[ii]])
    old_sis <- older_sis_dist(U[[ii]],F[[ii]])
    young_sis <- younger_sis_dist(U[[ii]],F[[ii]])
    old_niece <- older_nieces_dist(U[[ii]],F[[ii]])
    young_niece <- younger_nieces_dist(U[[ii]],F[[ii]])
    older_aunt <- older_aunts_dist(U[[ii]],F[[ii]])
    young_aunt <- younger_aunts_dist(U[[ii]],F[[ii]])

    kin_list <- list(daught,grand_daught,moth,grand_moth,
                     old_sis,young_sis,old_niece,young_niece,
                     older_aunt,young_aunt)


    df_list <- list()

    for( i in 1 : length(kin_list)){
      kin_member <- list_dists[[i]]

      df <- data.frame( kin_list[[i]] )
      dims <- dim( kin_list[[i]] )
      nr <- dims[1]
      nc <- dims[2]
      col_names <- rep(0, nc)
      for(i in 1 : nc){col_names[i] <- paste0("foc_age", i)}
      colnames(df) <- col_names
      cum_sum_kin <- apply(df[1 : (nr/2), ], 2, sum)
      deaths_kin <- apply(df[ (nr/2+1) : nr, ] , 2 ,sum)
      df$cum_sum_kin <- cum_sum_kin
      df$death_kin <- deaths_kin
      df$age_kin <- seq(1 , nr, by = 1)
      df <- df%>%select(-contains("NA"))
      ### Next melt the data to plot via groups and ages

      df <- df%>%select(c("cum_sum_kin" , "death_kin" , "age_kin" ) )
      df <- melt(df, id = "age_kin")


      df$group <- rep( kin_member , nrow(df) )
      df$year <- rep(j, nrow(df))
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }

  df_year_list <- do.call("rbind", df_year_list)
  return(df_year_list)
}


## create a distribution df, with inputs a list of matrices
## from the above functions, and a list of ages to be analysed...



## Order matters: a,b,d,g,h,m,n,p,q
# are daughters, grand daughters, mothers, grand mothers, older sis, younger sis,
# nieces from older sis, nieces from younger sis, aunts older than mother,
# aunts younger than mother

create_full_distributions_df <- function(U, F, list_dists, years, start_year ){

  df_year_list <- list()
  for(j in years){
    ii <- as.numeric(j) - start_year + 1
    ## a,b,d,g,m,n
    daught <- daughter_dist(U[[ii]],F[[ii]])
    grand_daught <- grand_daughter_dist(U[[ii]],F[[ii]])
    moth <- mother_dist(U[[ii]],F[[ii]])
    grand_moth <- grand_mother_dist(U[[ii]],F[[ii]])
    old_sis <- older_sis_dist(U[[ii]],F[[ii]])
    young_sis <- younger_sis_dist(U[[ii]],F[[ii]])
    old_niece <- older_nieces_dist(U[[ii]],F[[ii]])
    young_niece <- younger_nieces_dist(U[[ii]],F[[ii]])
    older_aunt <- older_aunts_dist(U[[ii]],F[[ii]])
    young_aunt <- younger_aunts_dist(U[[ii]],F[[ii]])

    kin_list <- list(daught,grand_daught,moth,grand_moth,
                     old_sis,young_sis,old_niece,young_niece,
                     older_aunt,young_aunt)


    df_list <- list()
    for( i in 1 : length(kin_list)){

      kin_member <- list_dists[[i]]

      df <- data.frame( kin_list[[i]] )
      dims <- dim( kin_list[[i]] )
      nr <- dims[1]
      nc <- dims[2]
      col_names <- rep(0, nc)
      for(i in 1 : nc){col_names[i] <- paste0("foc_age", i)}
      colnames(df) <- col_names
      cum_sum_kin <- apply(df[1 : (nr/2), ], 2, sum)
      deaths_kin <- apply(df[ (nr/2+1) : nr, ] , 2 ,sum)
      df$cum_sum_kin <- cum_sum_kin
      df$death_kin <- deaths_kin
      df$age_kin <- seq(1 , nr, by = 1)

      df <- df%>%select(-contains("NA"))
      ### Next melt the data to plot via groups and ages
      df <- melt(df, id = "age_kin")

      df$age_focal <- as.numeric(gsub("foc_age", "", df$variable))
      df$group <- rep( kin_member , nrow(df) )
      df$year <- rep(j, nrow(df))
      df_list[[length(df_list)+1]] <- df
    }

    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }

  df_year_list <- do.call("rbind", df_year_list)
  return(df_year_list)
}



############################################################ end of functions


