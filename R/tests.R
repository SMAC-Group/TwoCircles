#' Dixon test statistic based on spacing frequencies
#'
#' @param \code{Sk} spacing frequencies
#' @return Dixon test statistic
#' @export
#' @examples
#' data(pigeons)
#' spacings = compute_Sk(pigeons$control, pigeons$experimental)
#' dixon(spacings$Sk)
#' dixon(c(1,0,0,0,1,2))
dixon = function(Sk){
  return(sum(Sk^2))
}

#' Wilcoxon test statistic based on spacing frequencies
#'
#' @param \code{Sk} spacing frequencies
#' @return Wilcoxon test statistic
#' @export
#' @examples
#' data(pigeons)
#' spacings = compute_Sk(pigeons$control, pigeons$experimental)
#' wilcoxon(spacings$Sk)
#' wilcoxon(c(1,0,0,0,1,2))
wilcoxon = function(Sk){
  m = length(Sk)
  n = sum(Sk)
  k = 1:m
  return(sum((m - k + m*(m-1)/(2*n))*Sk))
}

#' Rao test statistic based on spacing frequencies
#'
#' @param \code{Sk} spacing frequencies
#' @return Rao test statistic
#' @export
#' @examples
#' data(pigeons)
#' spacings = compute_Sk(pigeons$control, pigeons$experimental)
#' rao(spacings$Sk)
#' rao(c(1,0,0,0,1,2))
rao = function(Sk){
  m = length(Sk)
  n = sum(Sk)
  k = 1:m
  return(sum( abs(Sk - n/m) ))
}

#' van der Waerden test statistic based on spacing frequencies
#'
#' @param \code{Sk} spacing frequencies
#' @return van der Waerden test statistic
#' @export
#' @examples
#' data(pigeons)
#' spacings = compute_Sk(pigeons$control, pigeons$experimental)
#' van.der.Waerden(spacings$Sk)
#' van.der.Waerden(c(1,0,0,0,1,2))
van.der.Waerden = function(Sk){
  m = length(Sk)
  n = sum(Sk)
  k = 1:m
  return(sum( qnorm(k/(m+1))*Sk ))
}

#' Savage test statistic based on spacing frequencies
#'
#' @param \code{Sk} spacing frequencies
#' @return Savage test statistic
#' @export
#' @examples
#' data(pigeons)
#' spacings = compute_Sk(pigeons$control, pigeons$experimental)
#' savage(spacings$Sk)
#' savage(c(1,0,0,0,1,2))
savage = function(Sk){
  m = length(Sk)
  n = sum(Sk)
  k = 1:m
  return((-1)*sum( log(1 - k/(m+1))*Sk ))
}

#' Wheeler-Watson test statistic based on spacing frequencies
#'
#' @param \code{Sk} spacing frequencies
#' @return Wheeler-Watson test statistic
#' @export
#' @examples
#' data(pigeons)
#' spacings = compute_Sk(pigeons$control, pigeons$experimental)
#' wheeler_watson(spacings$Sk)
#' wheeler_watson(c(1,0,0,0,1,2))
wheeler_watson = function(Sk){
  x = Sk
  m = length(x)
  n = sum(x)
  V1 = rep(NA,m)
  V2 = rep(NA,m)

  cte = 2*pi/(m + n)

  V1[1] = cos(cte)
  V2[1] = sin(cte)

  for (k in 2:m){
    Rk = k + sum(x[1:(k-1)])
    V1[k] = cos(cte*Rk)
    V2[k] = sin(cte*Rk)
  }

  stat = (sum(V1))^2 + (sum(V2))^2
  stat
}
