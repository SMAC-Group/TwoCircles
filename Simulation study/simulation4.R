# Load packages
library(CircStats)
library(TwoCircles)

# Sample size
n1 = c(5,5)
n2 = c(18,24)

# Define alpha
alpha = 0.05

# Define mu
mu.length = 15
proba = seq(from = 0, to = 1, length.out = mu.length)

# Number of bootstrap replication
B = 10^4

# Progress bar initialisation
total = 2*B*length(proba)
pb = txtProgressBar(min = 0, max = total, style = 3)
counter = 0

# Results
res_rao = res_dixon = res_WW = res_wil = matrix(NA, 2, mu.length)

for (k in 1:2){
  n = n2[k]
  m = n1[k] + 1

  # Get critical values
  C.rao = get_critical_values(n = n, m = m, test = "rao", alpha = alpha)$brackets$c2
  C.dixon = get_critical_values(n = n, m = m, test = "dixon", alpha = alpha)$brackets$c2
  C.WW = get_critical_values(n = n, m = m, test = "ww", alpha = alpha)$brackets$c2
  C.Wil = get_critical_values(n = n, m = m, test = "wilcox", alpha = alpha)$brackets$c2

  # Initialisation of results matrix
  rao_results = matrix(NA, B, mu.length)
  dixon_results = matrix(NA, B, mu.length)
  WW_results = matrix(NA, B, mu.length)
  Wil_results = matrix(NA, B, mu.length)

  # Start simu
  for (i in 1:length(proba)){
    for (j in 1:B){
      # Control seed
      set.seed(j)

      # Simulation sample
      x = rnorm(n1[k], 0, 1)
      y1 = rnorm(n2[k],  0, 1)
      y2 = rnorm(n2[k],  pi, 1)
      y3 = rnorm(n2[k], -pi, 1)

      w1 = rbinom(n = n2[k], size = 1, prob = 0.5)
      w2 = rbinom(n = n2[k], size = 1, prob = proba[i])
      y = (w2 == 1)*((w1==0)*y2 + (w1==1)*y3) + (w2 == 0)*y1

      # Compute Rao Stat
      rao_results[j,i] = circular_test(x, y, test = "rao", type = "mc", B = 10, circle = FALSE)$stat >= C.rao

      # Compute Dixon Stat
      dixon_results[j,i] = circular_test(x, y, test = "dixon", type = "mc", B = 10, circle = FALSE)$stat >= C.dixon

      # Compute WW Stat
      WW_results[j,i] = circular_test(x, y, test = "ww", type = "mc", B = 10, circle = FALSE)$stat >= C.WW

      # Compute Wilcox stat
      Wil_results[j,i] = circular_test(x, y, test = "wilcox", type = "mc", B = 10, circle = FALSE)$stat >= C.Wil

      # Update progress bar
      counter = counter + 1
      setTxtProgressBar(pb, counter)
    }
  }

  # Save results
  res_rao[k,] = 100*(apply(rao_results,2,mean))
  res_dixon[k,] = 100*(apply(dixon_results,2,mean))
  res_WW[k,] = 100*(apply(WW_results,2,mean))
  res_wil[k,] = 100*(apply(Wil_results,2,mean))
}

# Close progress bar
close(pb)

res_rao
res_dixon
res_WW
res_wil

# Save results
out = list(res_rao = res_rao, res_dixon = res_dixon, res_WW = res_WW, res_wil = res_wil)
save(out, file = "Simulation study/simulation4.RData")

