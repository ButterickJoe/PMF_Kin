

# The growth rate -- the spectral radius of PM
lambda <- function(PM) {
  lead_eig <- (abs(eigen(PM, only.values = TRUE)$values))
  lead_eig <- lead_eig[which.max(lead_eig)]
  return(lead_eig)
}

SD <- function(PM) {
  spectral_stuff <- eigen(PM)
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)
}
# Find the row-eigenvector corresponding to the spectral radius (Stable reproductive values in Demographics)
RD <- function(PM) {
  spectral_stuff <- eigen(t(PM))
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)
}

###################################################### Useful matrix operations

## Constructing a unit vector with a 1 in the ith position
e_vector <- function(i, n){
  e <- rep(0, n)
  e[i] <- 1
  return(e)
}
## Creating a matrix of zeros with a 1 in the i,j-th entry
E_matrix <- function(i,j,n,m){
  E <- Matrix::Matrix(nrow = (n), ncol = (m), data = 0, sparse = TRUE)
  E[i,j] <- 1
  return(E)

}
## Creating the Vec-commutation matrix
K_perm_mat <- function(n,m){
  perm <- Matrix::Matrix(nrow = (n*m), ncol = (n*m), data = 0, sparse = TRUE)
  for(i in 1:n){
    for(j in 1:m){
      perm = perm + kronecker( E_matrix(i,j,n,m) , Matrix::t(E_matrix(i,j,n,m)) )
    }
  }
  return(perm)
}


