########### Convolution powers

#' Title
#'
#' @param n number of distributions
#' @param list_ferts list of the distributions
#'
#' @return Convolution over the distributions (Eq 3 in text)

convoluion_nth <- function(n, list_ferts){
  if(n == 1){ conv = list_ferts[[1]] }
  if(n == 2){ conv = convolve(list_ferts[[(1)]], list_ferts[[2]], conj = FALSE, type = "c") }
  if(n > 2){
    conv = convolve(list_ferts[[(n-1)]], list_ferts[[n]], conj = FALSE, type = "c")
    for(i in (n-2):(1)){
      conv = convolve(list_ferts[[i]], conv, conj = FALSE, type = "c")
    }
  }
  return(conv)
}


