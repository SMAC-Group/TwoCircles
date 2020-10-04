# Load packages
library(CircStats)
library(TwoCircles)

# Sample size (8 and 8)
n1 = c(6,6)
n2 = c(12,16)

# Define alpha
alpha = 0.05

# Define mu
mu.max = 20
mu.length = 15
delta = seq(from = 0, to = mu.max, length.out = mu.length)

# Number of bootstrap replication
B = 10^4

# Progress bar initialisation
total = 2*B*length(delta)
pb = txtProgressBar(min = 0, max = total, style = 3)
counter = 0

# Results
res_rao = res_dixon = res_WW = matrix(NA, 2, mu.length)

for (k in 1:2){
  n = n2[k]
  m = n1[k]

  # Get critical values
  C.rao = get_critical_values(n = n, m = m, test = "rao", alpha = alpha)$brackets$c2
  C.dixon = get_critical_values(n = n, m = m, test = "dixon", alpha = alpha)$brackets$c2
  C.WW = get_critical_values(n = n, m = m, test = "ww", alpha = alpha)$brackets$c2

  # Initialisation of results matrix
  rao_results = matrix(NA, B, mu.length)
  dixon_results = matrix(NA, B, mu.length)
  WW_results = matrix(NA, B, mu.length)

  # Start simu
  for (i in 1:length(delta)){
    for (j in 1:B){
      # Control seed
      set.seed(j)

      # Simulation sample
      #x = rvm(n1[k], 0, 2)
      x = runif(n1[k], 0, 2*pi)
      y = rvm(n2[k], 0, 0.000001 + delta[i])

      # Compute Rao Stat
      rao_results[j,i] = circular_test(x, y, test = "rao", type = "mc", B = 10)$stat >= C.rao

      # Compute Dixon Stat
      dixon_results[j,i] = circular_test(x, y, test = "dixon", type = "mc", B = 10)$stat >= C.dixon

      # Compute van der Waerden Stat
      WW_results[j,i] = circular_test(x, y, test = "ww", type = "mc", B = 10)$stat >= C.WW

      # Update progress bar
      counter = counter + 1
      setTxtProgressBar(pb, counter)
    }
  }

  # Save results
  res_rao[k,] = 100*(apply(rao_results,2,mean))
  res_dixon[k,] = 100*(apply(dixon_results,2,mean))
  res_WW[k,] = 100*(apply(WW_results,2,mean))
}

# Close progress bar
close(pb)

# Save results
out = list(res_rao = res_rao, res_dixon = res_dixon, res_WW = res_WW)
save(out, file = "Simulation study/simulation2.RData")
res_rao
res_WW
