#' Monte-carlo simulation function to compute pvalues
#'
#' This function approximate by Monte-carlo simulation the distribution
#' of a given test statistic under the H0.
#' @param \code{stat} observed test statistics
#' @param \code{n} total of frequencies
#' @param \code{m} number of spacings
#' @return A list with the following structure:
#' \describe{
#'  \item{dist}{a \code{vector} containing the approximate distribution of the test statistic}
#'  \item{n}{total of frequencies}
#'  \item{m}{number of spacings}
#'  \item{test}{the considered test}
#'  }
#' @export
MC_pvalue = function(stat, n, m, test, B = 10^4, seed = 1982){
  set.seed(seed)
  res = rep(NA, B)

  if (test == "dixon"){
    for (i in 1:B){
      x = runif(m-1)
      y = runif(n)
      res[i] = dixon(compute_Sk(x,y)$Sk)
    }
  }

  if (test == "wilcox"){
    for (i in 1:B){
      x = runif(m-1)
      y = runif(n)
      res[i] = wilcoxon(compute_Sk(x,y)$Sk)
    }
  }

  if (test == "rao"){
    for (i in 1:B){
      x = runif(m-1)
      y = runif(n)
      res[i] = rao(compute_Sk(x,y)$Sk)
    }
  }

  if (test == "vdw"){
    for (i in 1:B){
      x = runif(m-1)
      y = runif(n)
      res[i] = van.der.Waerden(compute_Sk(x,y)$Sk)
    }
  }

  if (test == "ww"){
    for (i in 1:B){
      x = runif(m-1)
      y = runif(n)
      res[i] = wheeler_watson(compute_Sk(x,y)$Sk)
    }
  }

  if (test == "savage"){
    for (i in 1:B){
      x = runif(m-1)
      y = runif(n)
      res[i] = savage(compute_Sk(x,y)$Sk)
    }
  }
  list(dist = res, n = n, m = m, test = test)
}

#' Two sample Test of Homogeneity
#'
#' This function allows to perform various nonparametric test for homogeneity on two samples of circular data.
#' @param \code{x} first sample
#' @param \code{y} second sample
#' @param \code{test} considered test (\code{dixon} for Dixon test (default);
#' \code{ww} for Wheeler-Watson; \code{wilcox} for Wilcoxon test; \code{rao}
#' for Rao test; \code{vdw} for van der Waerden test; \code{savage}
#' for Savage test
#' @param \code{alpha} significance level (default = 0.05)
#' @param \code{B} number of bootstrap replications
#' @param \code{seed} seed used for random number generation
#' @param \code{type} method to compute pvalues (available methods: \code{exact} for exact computation (default) which is appropriate for small sample sizes;
#' \code{mc} for approximation based on Monte-carlo simulations)
#' @return A list with the following structure:
#' \describe{
#'  \item{cv}{a \code{list} containing the results of the exact distribution (\code{NULL} if \code{type = "mc"}), see function \code{get_critical_values} for details}
#'  \item{mc}{a \code{list} containing the results of the approximated distribution obtained by simulation (\code{NULL} if \code{type = "exact"}), see function \code{MC_pvalue} for details}
#'  \item{B}{number of bootstrap replications}
#'  \item{alpha}{significance level}
#'  \item{test}{the considered test}
#'  \item{stat}{observed test statistic}
#'  \item{spacings}{a \code{list} containing the observed spacings, see function \code{compute_Sk} for details}
#'  \item{x}{expression deparsing of the first dataset}
#'  \item{y}{expression deparsing of the second dataset}
#'  }
#' @export
#' @examples
#' # Load dataset
#' data(pigeons)
#' # Dixon test (exact pvalue)
#' circular_test(pigeons$experimental, pigeons$control)
#'
#' # Dixon test (approximated pvalue)
#' circular_test( pigeons$experimental, pigeons$control, type = "mc")
circular_test = function(x, y, test = "dixon", alpha = 0.05, B = NULL, type = "exact", seed = 1982){
  spacings = compute_Sk(x,y)

  if (sum(c("exact","mc") %in% type) == 0){
    stop("Unsupported parameter type; available options: exact and mc")
  }

  if (sum(c("dixon", "wilcox", "rao", "vdw", "savage", "ww") %in% test) == 0){
    stop("Unsupported test; available methods are dixon, ww, wilcox, rao, vdw and savage")
  }

  if (test == "dixon"){
    stat = dixon(spacings$Sk)
  }

  if (test == "wilcox"){
    stat = wilcoxon(spacings$Sk)
  }

  if (test == "rao"){
    stat = rao(spacings$Sk)
  }

  if (test == "vdw"){
    stat = van.der.Waerden(spacings$Sk)
  }

  if (test == "ww"){
    stat = wheeler_watson(spacings$Sk)
  }

  if (test == "savage"){
    stat = savage(spacings$Sk)
  }

  if (type == "exact" && (spacings$n + spacings$m) <= 30 && min(spacings$n,spacings$m) >=5){
    cv = get_critical_values(spacings$n, spacings$m, test = test, alpha = alpha)
    mc = NULL
  }else{
    if (type != "mc"){
      warning("type was changed to mc (instead of exact).")
    }
    if (is.null(B)){
      B = 10^4
    }
    cv = NULL
    mc = MC_pvalue(stat, spacings$n, spacings$m, test, B = B, seed = seed)
  }

  out = list(cv = cv, mc = mc, B = B, alpha = alpha, test = test, stat = stat,
             spacings = spacings, x = deparse(substitute(x)), y = deparse(substitute(y)))
  class(out) = "circular_test"
  out
}

#' @title Print Two sample Test of Homogeneity
#' @description
#' Prints the results of a two sample test of homogeneity
#' @method print circular_test
#' @export
#' @param x A \code{circular_test} object
#' @return Prints the results of a two sample test of homogeneity
#' @keywords internal
#' @examples
#' # Load dataset
#' data(pigeons)
#' # Dixon test (exact pvalue)
#' circular_test(pigeons$experimental, pigeons$control)
print.circular_test = function(out){
  test_list = c("dixon", "wilcox", "rao", "vdw", "savage", "ww")
  test_name = c("Dixon", "Wilcoxon", "Rao", "van der Waerden", "Savage", "Wheeler-Watson")
  cat(paste("\n      ", test_name[test_list %in% out$test], " Two Sample Test\n\n" , sep = ""))
  cat(paste("Data:  ", out$x, " and ", out$y, "\n", sep = ""))
  cat(paste("Test Statistic: ", round(out$stat,5), "\n", sep =""))

  if (is.null(out$mc)){
    cat(paste("Exact P-value: ", round(mean(out$stat <= out$cv$dist), 5), "\n", sep = ""))
    cat("Bracketing Points and Pair of Signif. Levels:\n")
    cat(paste("c1 = ",round(out$cv$brackets$c1,4) , " (p1 = " ,formatC( round(out$cv$brackets$p1,4), format='f', digits=4)  , ")\n", sep = ""))
    cat(paste("c2 = ", round(out$cv$brackets$c2,4) , " (p2 = " ,formatC( round(out$cv$brackets$p2,4), format='f', digits=4)  , ")\n", sep = ""))
  }else{
    cat(paste("Approx. P-value: ", round((1 + sum(out$stat <= out$mc$dist))/(length(out$mc$dist) + 1), 5), "\n", sep = ""))
    cat(paste("P-value stand. error: ", round(np_boot_pval(x = out$mc$dist, cut = out$stat), 7), "\n", sep = ""))
    cat(paste("based on ", out$B, " Monte-Carlo replications", sep = ""))
  }
}

np_boot_pval = function(x, cut, B = 10000){
  res = rep(NA, B)
  n = length(x)
  for (i in 1:B){
    xstar = sample(x, replace = TRUE)
    res[i] = (1 + sum(xstar <= cut))/(n + 1)
  }
  sd(res)
}

#' @title Plot null distribution associated to two sample test of Homogeneity
#' @description
#' Plots the (exact or approximated) distribution of the test statistics under H0.
#' @method plot circular_test
#' @export
#' @param x A \code{circular_test} object
#' @return Plots the (exact or approximated) distribution of the test statistics under H0.
#' @keywords internal
#' @examples
#' data(pigeons)
#' plot(circular_test(pigeons$experimental, pigeons$control))
plot.circular_test = function(x, cex.main = 1){
  out = x
  test_list = c("dixon", "wilcox", "rao", "vdw", "savage", "ww")
  test_name = c("Dixon", "Wilcoxon", "Rao", "van der Waerden", "Savage", "Wheeler-Watson")

  if (is.null(out$mc)){
    titre = paste("Exact distribution for ", test_name[test_list %in% out$test],
                  " test with n = ", out$spacings$n, " and m = ", out$spacings$m, "\n Test Statistic: ",
                  round(out$stat,4), " ; P-value: ", round(mean(out$stat <= out$cv$dist), 4), sep = "")
    plot(out$cv, main = titre, cex.main = cex.main)
  }else{
    titre = paste("Approx. distribution for ", test_name[test_list %in% out$test],
          " test with n = ", out$spacings$n, " and m = ", out$spacings$m, "\n Test Statistic: ",
          round(out$stat,4), " ; P-value: ", round(mean(out$stat <= out$mc$dist), 4), sep = "")
    plot.critical_values(out$mc, main = titre, cex.main = cex.main)
  }
  abline(v = out$stat, lty = 2, col = "red")
  #legend("topright", "Test Statistic", lty = 2, lwd = 1, col = "red", bty = "n")
}

#' Compute spacing frequencies
#'
#' @param \code{x} first sample
#' @param \code{y} second sample
#' @return A list with the following structure:
#' \describe{
#'  \item{Sk}{a \code{vector} containing the spacing frequencies}
#'  \item{n}{total of frequencies}
#'  \item{m}{number of spacings}
#' }
#' @export
#' @examples
#' data(pigeons)
#' compute_Sk(pigeons$control, pigeons$experimental)
compute_Sk = function(x,y){
  x = sort(x)
  n = length(y)
  m = length(x) + 1
  res = rep(NA,m)

  res[1] = sum(y <= x[1])
  res[m] = sum(y > x[(m-1)])
  for (i in 2:(m-1)){
    res[i] = sum((y <= x[i])*(y > x[(i-1)]))
  }
  list(Sk = res, n = n, m = m)
}

#' All permutations used to compute the exact distributions of circular tests
#'
#' @param \code{n} total of frequencies
#' @param \code{m} number of spacings
#' @return A \code{matrix} with all possible permutations
#' @export
#' @importFrom gRbase combnPrim
#' @examples
#' all_permutations(2,2)
#' all_permutations(4,3)
all_permutations = function (n = 2, m = 2)
{
  comb = combnPrim(n + m - 1, m - 1)
  comb = rbind(comb, n + m)
  if (m > 1){
    i = 2:m
    comb[i, ] = comb[i, ] - comb[i - 1, ]
  }
  t(comb - 1)
}

#' @title Compute Bracketing Values
#' @description
#' Compute bracketing values (c1, c2) corresponding to significance levels (p1, p2)
#' based on the significance level alpha (such that p1 > alpha > p2)
#' @export
#' @param x A \code{vector} of test statistics
#' @param alpha Significance level
#' @return Bracketing values (c1, c2) corresponding to significance levels (p1, p2)
#' @keywords internal
get_bracket = function(x, alpha){
  nb = length(x)
  q.alpha = round(1 + (length(x))*(1 - alpha))
  prob = (1 - sum(x < x[q.alpha])/nb)

  if (prob > alpha){
    p1 = prob
    c1 = x[q.alpha]

    for (i in 1:(nb - 1 - q.alpha)){
      if (abs(x[q.alpha] - x[(q.alpha+i)]) > 10^(-10)){
        c2 = x[(q.alpha+i)]
        break
      }
    }

    p2 = (1 - sum(x < x[(q.alpha+i)])/length(x))

  }else{
    p2 = prob
    c2 = x[q.alpha]

    for (i in 1:(q.alpha - 1)){
      if (abs(x[q.alpha] - x[(q.alpha-i)]) > 10^(-10)){
        c1 = x[(q.alpha-i)]
        break
      }
    }
    p1 = (1 - sum(x < x[(q.alpha-i)])/length(x))
  }
  return(list(c1 = c1, c2 = c2, p1 = p1, p2 = p2))
}

#' Compute critical values for circular tests
#'
#' @param \code{n} total of frequencies
#' @param \code{m} number of spacings
#' @param \code{test} considered test (\code{dixon} for Dixon test (default);
#' \code{ww} for Wheeler-Watson; \code{wilcox} for Wilcoxon test; \code{rao}
#' for Rao test; \code{vdw} for van der Waerden test; \code{savage}
#' for Savage test)
#' @param alpha Significance level
#' @return A list with the following structure:
#' \describe{
#'  \item{dist}{test distribution under H0}
#'  \item{bracket}{brackets values and associated p-values}
#'  \item{test}{considered test}
#'  \item{alpha}{considered alpha}
#'  \item{n}{total of frequencies}
#'  \item{m}{number of spacings}
#' }
#' @author Stephane Guerrier
#' @export
#' @examples
#' crit_val = get_critical_values(6,8)
#' crit_val
#' plot(crit_val)
#' crit_val = get_critical_values(12, 8, test = "ww")
#' crit_val
#' plot(crit_val)
get_critical_values = function(n, m, test = "dixon", alpha = 0.05){

  if (sum(c("dixon", "wilcox", "rao", "vdw", "savage", "ww") %in% test) == 0){
    stop("Unsupported test; available methods are dixon, ww, wilcox, rao, vdw and savage")
  }

  if ((n + m) > 30){
    stop("Unable to compute critical values (n and/or m are too large)")
  }

  if (min(n,m) < 5){
    stop("Unable to compute critical values (n and/or m are too small)")
  }


  A = all_permutations(n, m)

  if (test == "dixon"){
    x = sort(apply(A,1,dixon))
  }

  if (test == "wilcox"){
    x = sort(apply(A,1,wilcoxon))
  }

  if (test == "rao"){
    x = sort(apply(A,1,rao))
  }

  if (test == "vdw"){
    x = sort(apply(A,1,van.der.Waerden))
  }

  if (test == "ww"){
    x = sort(apply(A,1,wheeler_watson))
  }

  if (test == "savage"){
    x = sort(apply(A,1,savage))
  }

  brackets = get_bracket(x, alpha = alpha)
  out = list(dist = x, brackets = brackets, test = test, alpha = alpha, n = n, m = m)
  class(out) = "critical_values"
  out
}

#' @title Print Bracketing Values for Two Sample Test of Homogeneity
#' @description
#' Prints bracketing values and corresponding to significance levels
#' for a two sample test of homogeneity
#' @method print critical_values
#' @export
#' @param x A \code{critical_values} object
#' @return Prints bracketing values and corresponding to significance levels
#' for a two sample test of homogeneity
#' @keywords internal
#' @examples
#' # Load dataset
#' get_critical_values(6,8)
print.critical_values = function(x){
  test_list = c("dixon", "wilcox", "rao", "vdw", "savage", "ww")
  test_name = c("Dixon", "Wilcoxon", "Rao", "van der Waerden", "Savage", "Wheeler-Watson")
  cat("Bracketing values (c1, c2) corresponding to significance levels (p1, p2)\n")
  cat(paste("for ", test_name[test_list %in% x$test], " test based on the significance level ", x$alpha, "\n", sep =""))
  cat(paste("c1 = ",round(x$brackets$c1,4), " (p1 = " ,round(x$brackets$p1,4) , ")\n", sep = ""))
  cat(paste("c2 = ",round(x$brackets$c2,4), " (p2 = " ,round(x$brackets$p2,4) , ")\n", sep = ""))
}

#' @title Plot null distribution associated to two sample test of Homogeneity
#' @description
#' Plots the (exact or approximated) distribution of the test statistics under H0.
#' @method plot critical_values
#' @export
#' @param x A \code{critical_values} object
#' @return Plots the (exact or approximated) distribution of the test statistics under H0.
#' @keywords internal
#' @examples
#' plot(get_critical_values(6,8))
plot.critical_values = function(x, main = NULL, cex.main = 1){
  test_list = c("dixon", "wilcox", "rao", "vdw", "savage", "ww")
  test_name = c("Dixon", "Wilcoxon", "Rao", "van der Waerden", "Savage", "Wheeler-Watson")
  tab = table(x$dist)
  xx = as.numeric(names(tab))
  p = as.numeric(tab/sum(tab))
  if (is.null(main)){
      main =  paste("Exact distribution for ", test_name[test_list %in% x$test],
                     " test with n = ", x$n, " and m = ", x$m, sep = "")
  }

  plot(xx, p, type="h", xlim=range(xx), ylim= c(0, max(p)), lwd=2, col="blue", ylab="p(x)",
       xlab = "T", main = main, cex.main = cex.main)
  points(xx,p,pch=16,cex=1.5,col="dark red")
}
