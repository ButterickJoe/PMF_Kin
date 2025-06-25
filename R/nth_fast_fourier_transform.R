#' Title
#'
#' @param n number of distributions
#' @param list_ferts list of the distributions
#'
#' @return Convolution over the distributions (Eq 3 in text) by fft inverse

fourier_fft_nth <- function(n, list_ferts){
  if(n == 1){ inv_f = list_ferts[[1]] }
  if(n == 2){ inv_f = gsignal::ifft(fft(list_ferts[[1]])*fft(list_ferts[[2]])) }
  if(n > 2){
    inv_f = fft(list_ferts[[(n-1)]])*fft(list_ferts[[n]])
    for(i in (n-2):(1)){
      inv_f = fft(list_ferts[[i]])*inv_f
    }
    inv_f = gsignal::ifft(inv_f)
  }
  return(Re(inv_f))
}
